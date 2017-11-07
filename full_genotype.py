from __future__ import division
import numpy as np
import argparse
import vcf
import itertools
import sys
import multiprocessing
import pandas as pd
from collections import *
from tools.bio import *
from cat.plots import *
from phase_lib import *


def parse_args():
    parser = argparse.ArgumentParser()
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--bam', help='input BAM')
    input_group.add_argument('--pileup', help='Path to save parsed and converted pileup if desired')
    parser.add_argument('--save-pileup', help='Path to save pileup to')
    parser.add_argument('--features', required=True, help='Features file. Must have column names in the format <genome>_<paralog><#>')
    parser.add_argument('--paratype-pseudo', default=1, type=float)
    parser.add_argument('--read-pseudo', default=0, type=float)
    parser.add_argument('--ratio-plot')
    parser.add_argument('--name')
    parser.add_argument('--inferred-copy-numbers', default=['4', '2', '2', '2'], nargs=4,
                        help='Space delimited list of copy numbers in the order AB C D N')
    parser.add_argument('--consensus-fasta', default='/hive/users/cbosworth/refs/notch/notch2_aligned_consensus.fasta')
    parser.add_argument('--pileup-converter', default='/cluster/home/ifiddes/pileup2base/pileup2base.pl')
    return parser.parse_args()


def construct_C(inferred_copy_numbers, filtered_features):
    """Construct C based on previous information of counts. Only allow possibilties to be enumerated that match this information"""

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

    # construct all possibilities for each feature group
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


if __name__ == '__main__':
    args = parse_args()
    if args.bam is not None:
        pileup_recs = make_pileup(args.bam)
    else:
        pileup_recs = parse_pileup(args.pileup)
    df = convert_pileup(pileup_recs, args.pileup_converter)
    _, seq = read_fasta(args.consensus_fasta, None).next()
    data = parse_converted_pileup(df, seq)
    if args.pileup is not None:
        data.to_csv(args.pileup, sep='\t')
    features = pd.read_csv(args.features, sep='\t', index_col=0)
    # find shared positions in case data is missing some
    positions = set(features.index) & set(data['loc'])
    filtered_data = data[data['loc'].isin(positions)]
    # filter features too
    filtered_features = features[features.index.isin(positions)]

    C = construct_C(args.inferred_copy_numbers, filtered_features)
    Ct = C.T

    num = np.dot(filtered_features, Ct)
    denom = np.sum(Ct, axis=0)

    S = (args.paratype_pseudo + num) / (2.0 * args.paratype_pseudo + denom)

    S_log = np.log(S)
    S_inv = np.log(1 - S)

    # M is the number of alt reads, N is the number of ref reads
    M = 1.0 * filtered_data.alt_count
    N = 1.0 * filtered_data.ref_count
    R = np.dot(M + args.read_pseudo, S_log) + np.dot(N + args.read_pseudo, S_inv)

    if args.ratio_plot is not None:
        with open(args.ratio_plot, 'w') as outf, PdfPages(outf) as pdf:
            g=sns.distplot(filtered_data.ratio, bins=50)
            g.set_xlim(0, 0.5)
            if args.name is not None:
                g.set_title(args.name)
            else:
                g.set_title(args.bam)
            multipage_close(pdf)

    R_map = {i: x for i, x in enumerate(R)}
    best_index, score = sorted(R_map.iteritems(), key=lambda x: x[1])[-1]
    best_haps = C[best_index]

    

    #print 'log odds: {}'.format(score)
    #print 'results: '
    #for x, y in zip(filtered_features.columns, best_haps):
    #    if y > 0:
    #        print '{}: {}'.format(x, y)

    ordered = sorted(R_map.iteritems(), key=lambda x: x[1])[-10:][::-1]
    print ' '.join(features.columns)
    for i, x in enumerate([[C[pos], pos, val] for pos, val in ordered], 1):
        print '{}: {}'.format(i, x)
