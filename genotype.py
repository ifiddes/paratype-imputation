from __future__ import division
import numpy as np
import argparse
import vcf
import itertools
import sys
import multiprocessing
import seaborn as sns
import pandas as pd
import cPickle as pickle
from tools.bio import *
from phase_lib import *


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', help="Maximum number of haplotypes", default=10, type=int)
    parser.add_argument('--pileup', help='input parsed pileup', required=True)
    parser.add_argument('--precomputed-s', required=True)
    parser.add_argument('--read-pseudo', default=0, type=float)
    parser.add_argument('--num-hits-to-print', default=10, type=int)
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    S_log, S_inv, features = pickle.load(open(args.precomputed_s))
    data = pd.read_csv(args.pileup, sep='\t', index_col=0)
    positions = set(features.index)

    # M is the number of alt reads, N is the number of ref reads
    M = 1.0 * filtered_data.alt_count
    N = 1.0 * filtered_data.ref_count
    R = np.dot(M + args.read_pseudo, S_log) + np.dot(N + args.read_pseudo, S_inv)

    R_map = {i: x for i, x in enumerate(R)}
    ordered = sorted(R_map.iteritems(), key=lambda x: x[1])[-10:][::-1]

    print ' '.join(features.columns)
    for i, x in enumerate([[C[pos], pos, val] for pos, val in ordered], 1):
        print '{}: {}'.format(i, x)
