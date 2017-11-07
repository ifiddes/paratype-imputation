from __future__ import division
import numpy as np
import argparse
import vcf
import itertools
import sys
import multiprocessing
import logging
import seaborn as sns
import pandas as pd
from glob import glob
from tools.fileOps import *
from tools.procOps import *
from tools.bio import *
from phase_lib import *

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bams', nargs='+', help='Phased haplotypes mapped to the consensus', required=True)
    parser.add_argument('--feature-file', help='File to write features to', required=True)
    parser.add_argument('--consensus-fasta', default='/hive/users/cbosworth/refs/notch/notch2_aligned_consensus.fasta')
    parser.add_argument('--pileup-converter', default='/cluster/home/ifiddes/pileup2base/pileup2base.pl')
    parser.add_argument('--cutoff', type=float, default=0.95, help='Percent of reads that must support alt to be considered')
    parser.add_argument('--min-cover', type=int, default=10, help='Minimum number of bases required for a feature to be considered')
    parser.add_argument('--cores', default=1, help='cores to use to speed up parsing', type=int)
    return parser.parse_args()


def construct_feature_vector(df, cutoff, min_cover):
    """For a haplotype phased genome, return a binary vector of all non-ref sites"""
    filtered = df[(df.coverage >= min_cover) & (df.ratio >= cutoff)]
    return set(filtered['loc'])


def wrapper_fn(m):
    """Wrapper for all steps to allow for parallelization"""
    # unpack, hack around map
    bam, cutoff, min_cover, seq, pileup_converter = m
    try:
        pileup_recs = make_pileup(bam)
    except Exception, e:
        raise RuntimeError('bam {} failed samtools conversion msg: {}'.format(bam, e))
    try:
        df = convert_pileup(pileup_recs, pileup_converter)
    except Exception, e:
        raise RuntimeError('bam {} failed conversion msg: {}'.format(bam, e))
    df = parse_converted_pileup(df, seq)
    feature_vector = construct_feature_vector(df, cutoff, min_cover)
    return feature_vector


if __name__ == '__main__':
    args = parse_args()

    _, seq = read_fasta(args.consensus_fasta, None).next()

    # convert each bam to a feature vector
    p = multiprocessing.Pool(processes=args.cores)
    # construct set of commands
    commands = []
    for bam in args.bams:
        commands.append([bam, args.cutoff, args.min_cover, seq, args.pileup_converter])
    feature_vectors = p.map(wrapper_fn, commands)
    p.close()
    p.join()

    positions = sorted(set.union(*feature_vectors))
    names = ['_'.join([os.path.dirname(x).split('/')[-1], os.path.basename(x).split('.')[0]]) for x in args.bams]
    feature_vector_dict = dict(zip(names, feature_vectors))

    features = []
    for p in positions:
        f = [int(p in feature_vector_dict[x]) for x in names]
        features.append([p] + f)

    features = pd.DataFrame(features, columns=['position'] + names)
    features = features.set_index('position')
    features.to_csv(args.feature_file, sep='\t')
