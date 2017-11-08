from __future__ import division
import numpy as np
import argparse
import vcf
import itertools
import sys
import random
import multiprocessing
import pandas as pd
from scipy.stats import norm
import cPickle as pickle
from collections import *
from tools.bio import *
from cat.plots import *
from phase_lib import *
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser()
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--bam', help='input BAM')
    input_group.add_argument('--pileup', help='Path to saved parsed and converted')
    parser.add_argument('--save-pileup', help='Path to save pileup to')
    parser.add_argument('--features', required=True, help='Features file. Must have column names in the format <genome>_<paralog><#>')
    parser.add_argument('--read-pseudo', default=10, type=float)
    parser.add_argument('--deviance-mask-cutoff', default=0.0025, type=float, help='Percent of most deviant sites to mask')
    parser.add_argument('--ratio-plot')
    parser.add_argument('--deviance-plot')
    parser.add_argument('--name')
    parser.add_argument('--blacklist')
    parser.add_argument('--precomputed-C', help='precomputed C pickle. Must match inferred copy number.')
    parser.add_argument('--inferred-copy-numbers', default=['4', '2', '2', '2'], nargs=4,
                        help='Space delimited list of copy numbers in the order AB C D N')
    parser.add_argument('--consensus-fasta', default='/hive/users/cbosworth/refs/notch/notch2_aligned_consensus.fasta')
    parser.add_argument('--pileup-converter', default='/cluster/home/ifiddes/pileup2base/pileup2base.pl')
    parser.add_argument('--feature-window', help='Position range in consensus to consider', nargs=2)
    parser.add_argument('--real-genotype', help='If we know the real genotype, place it here to report that score', nargs='+')
    return parser.parse_args()


def calculate_feature_groups(inferred_copy_numbers, filtered_features):
    # create map of paralog to inferred total
    inferred_copy_numbers = {x: int(y) for x, y in zip(['AB', 'C', 'D', 'N'], inferred_copy_numbers)}
    total_copies = sum(inferred_copy_numbers.values())
    features_by_paralog = [x.split('_')[1][0] for x in filtered_features.columns]

    # for simplicity, replace all A or B with AB
    features_by_paralog = [x if x not in ['A', 'B'] else 'AB' for x in features_by_paralog]

    # split them into groups, maintaining original positions
    feature_groups = defaultdict(list)
    for i, f in enumerate(features_by_paralog):
        feature_groups[f].append(i)
    return feature_groups


def construct_C(inferred_copy_numbers, filtered_features):
    """Construct C based on previous information of counts. Only allow possibilties to be enumerated that match this information"""
    feature_groups = calculate_feature_groups(inferred_copy_numbers, filtered_features)
    # construct all possibilities for each feature group
    inferred_copy_numbers = {x: int(y) for x, y in zip(['AB', 'C', 'D', 'N'], inferred_copy_numbers)}
    possibilities = {}
    for f, positions in feature_groups.iteritems():
        inferred_copy = inferred_copy_numbers[f]
        r = np.array([np.array(x) for x in itertools.product(range(inferred_copy + 1), repeat=len(positions))
                      if sum(x) == inferred_copy])
        possibilities[f] = r

    def array_product(a1, a2):
        m1,n1 = a1.shape
        m2,n2 = a2.shape
        out = np.zeros((m1, m2, n1 + n2), dtype=int)
        out[:,:,:n1] = a1[:,None,:]
        out[:,:,n1:] = a2
        out.shape = (m1 * m2, -1)
        return out

    abc = array_product(possibilities['AB'], possibilities['C'])
    abcd = array_product(abc, possibilities['D'])
    abcdn = array_product(abcd, possibilities['N'])

    # finally, rearrange the columns to reflect the original positioning
    order = feature_groups['AB'] + feature_groups['C'] + feature_groups['D'] + feature_groups['N']
    i = np.argsort(order)
    ordered = abcdn[:,i]

    return ordered


def calculate_paratype_pseudo(inferred_copy_numbers, filtered_features):
    feature_groups = calculate_feature_groups(inferred_copy_numbers, filtered_features)
    # construct pseudocounts bsaed on the features
    ab=np.sum(filtered_features[feature_groups['AB']], axis=1)
    c=np.sum(filtered_features[feature_groups['C']], axis=1)
    d=np.sum(filtered_features[feature_groups['D']], axis=1)
    n=np.sum(filtered_features[feature_groups['N']], axis=1)

    not_ambig = ((ab>0) & (c==0) & (d==0) & (n==0)) | \
                ((ab==0) & (c>0) & (d==0) & (n==0)) | \
                ((ab==0) & (c==0) & (d>0) & (n==0)) | \
                ((ab==0) & (c==0) & (d==0) & (n>0))

    paratype_pseudo = np.array([0.01 if x else 1 for x in not_ambig])
    return paratype_pseudo


def calculate_deviance(s, filtered_data):
    expected_alt = np.multiply(s, filtered_data['coverage'])
    expected_ref = filtered_data['coverage'] - expected_alt
    actual_alt = filtered_data['alt_count']
    actual_ref = filtered_data['ref_count']
    n = filtered_data['coverage']
    deviance = (expected_alt - actual_alt) / (np.sqrt(n * s * (1 - s)))
    return deviance


def filter_deviance(deviance):
    """Remove NaN from deviance. Used for plotting"""
    return deviance.replace([np.inf, -np.inf], np.nan).dropna()


def randomize_histogram(deviance):
    """Multiply deviance vector by coin toss vector"""
    v = [random.choice([-1, 1]) for _ in range(len(deviance))]
    return np.multiply(deviance, v)


def calculate_variance(deviance, filtered_data):
    return sum(np.multiply(deviance, deviance)) / len(filtered_data['coverage'])


def reject_deviance_outliers(deviance, deviance_mask_cutoff):
    """Rejects deviance_mask_cutoff of lowest deviance values, returns the worst"""
    s = len(deviance)
    d = int(round(deviance_mask_cutoff * s))
    indices = deviance.argsort()[:d]
    return indices


def calculate_mask(S, filtered_data, deviance_mask_cutoff):
    """Constructs the masking matrix based on deviance"""
    masked = []
    for s in S:
        deviance = calculate_deviance(s, filtered_data)
        to_mask = set(reject_deviance_outliers(deviance, deviance_mask_cutoff))
        mask = [0 if x in to_mask else 1 for x in range(len(s))]
        masked.append(mask)
    return np.array(masked)


def calculate_values(filtered_features, C, paratype_pseudo, read_pseudo, deviance_mask_cutoff):
    """Calculate S and the vector R."""
    Ct = C.T

    num = np.dot(filtered_features, Ct)
    denom = np.sum(Ct[:,0])  # first column -- all columns sum to total number of copies
    #paratype_pseudo += denom
    denom = [denom] * len(paratype_pseudo) + paratype_pseudo
    S = (paratype_pseudo + num.T) / ( (2.0 * paratype_pseudo) + denom)
    #S = S.T

    #num = np.dot(filtered_features, Ct)
    #denom = np.sum(Ct, axis=0)
    #S = (1 + num) / (2.0 * 1 + denom)

    S_log = np.log(S)
    S_inv = np.log(1 - S)
    # calculate the masking matrix based on deviance
    S_mask = calculate_mask(S, filtered_data, deviance_mask_cutoff)
    # mask these matrices
    S = np.multiply(S, S_mask)
    S_log = np.multiply(S_log, S_mask)
    S_inv = np.multiply(S_inv, S_mask)

    # M is the number of alt reads, N is the number of ref reads
    M = 1.0 * filtered_data.alt_count
    N = 1.0 * filtered_data.ref_count
    R = (np.dot(M + read_pseudo, S_log.T) + np.dot(N + read_pseudo, S_inv.T)).T
    return S, R


if __name__ == '__main__':
    args = parse_args()

    _, seq = read_fasta(args.consensus_fasta, None).next()

    if args.bam is not None:
        pileup_recs = make_pileup(args.bam)
        df = convert_pileup(pileup_recs, args.pileup_converter)
        data = parse_converted_pileup(df, seq)
        if args.save_pileup is not None:
            data.to_csv(args.save_pileup, sep='\t')
    else:
        data = pd.read_csv(args.pileup, sep='\t', index_col=0)

    features = pd.read_csv(args.features, sep='\t', index_col=0)
    # find shared positions in case data is missing some
    positions = set(features.index) & set(data['loc'])
    filtered_data = data[data['loc'].isin(positions)]
    # filter features too
    filtered_features = features[features.index.isin(positions)]

    if args.blacklist is not None:
        blacklist = map(int, [x.rstrip() for x in open(args.blacklist)])
        filtered_features = filtered_features[~filtered_features.index.isin(positions)]

    # subset features by position if requested
    if args.feature_window is not None:
        tmp = filtered_features.reset_index()
        start, stop = map(int, args.feature_window)
        feature_subset = tmp[(tmp['position'] >= start) & (tmp['position'] <= stop)]
        filtered_features = feature_subset.set_index('position')

    # use a pre-computed C to save time. Must have correct genotype
    if args.precomputed_C is not None:
        C = pickle.load(open(args.precomputed_C))
    else:
        C = construct_C(args.inferred_copy_numbers, filtered_features)

    paratype_pseudo = calculate_paratype_pseudo(args.inferred_copy_numbers, filtered_features)
    S, R = calculate_values(filtered_features, C, paratype_pseudo, args.read_pseudo, args.deviance_mask_cutoff)

    R_map = {i: x for i, x in enumerate(R)}
    best_index, score = sorted(R_map.iteritems(), key=lambda x: x[1])[-1]
    best_haps = C[best_index]
    ordered = sorted(R_map.iteritems(), key=lambda x: x[1])[-10:][::-1]

    if args.deviance_plot is not None:
        with open(args.deviance_plot, 'w') as outf, PdfPages(outf) as pdf:
            if args.real_genotype is not None:
                real_genome = np.array(map(int, args.real_genotype))
                index = np.where(np.all(C == real_genome, axis=1))[0][0]
                deviance = calculate_deviance(S[index], filtered_data)
                variance = calculate_variance(deviance, filtered_data)
                fig, ax = plt.subplots()
                ax.plot(deviance.index, deviance)
                fig.suptitle('Deviance by position for {} real genotype variance = {}'.format(args.name, variance))
                ax.set_ylim(-70, 10)
                multipage_close(pdf)
                fig, ax = plt.subplots()
                g = sns.distplot(filter_deviance(deviance), ax=ax, fit=norm, kde=False)
                fig.suptitle('Deviance distribution for {} real genotype variance = {}'.format(args.name, variance))
                multipage_close(pdf)

            for i in xrange(1, 4, 1):
                index = ordered[i][0]
                deviance = calculate_deviance(S[index], filtered_data)
                variance = calculate_variance(deviance, filtered_data)
                fig, ax = plt.subplots()
                ax.plot(deviance.index, deviance)
                fig.suptitle('Deviance by position for {} hit {} variance = {}'.format(args.name, i, variance))
                ax.set_ylim(-70, 10)
                multipage_close(pdf)
                fig, ax = plt.subplots()
                g = sns.distplot(filter_deviance(deviance), ax=ax, fit=norm, kde=False)
                fig.suptitle('Deviance distribution for {} hit {} variance = {}'.format(args.name, i, variance))
                multipage_close(pdf)

    if args.ratio_plot is not None:
        with open(args.ratio_plot, 'w') as outf, PdfPages(outf) as pdf:
            g=sns.distplot(filtered_data.ratio, bins=50)
            g.set_xlim(0, 0.5)
            if args.name is not None:
                g.set_title(args.name)
            else:
                g.set_title(args.bam)
            multipage_close(pdf)


    print 'log odds: {}'.format(score)
    print ''
    print 'results: '
    for x, y in zip(filtered_features.columns, best_haps):
        if y > 0:
            print '{}: {}'.format(x, y)

    print ' '.join(features.columns)
    print 'top 10 hits:'
    for i, x in enumerate([[C[pos], pos, val] for pos, val in ordered], 1):
        print '{}: {}'.format(i, x)
