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
    parser.add_argument('--features', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--paratype-pseudo', default=1, type=float)
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    features = pd.read_csv(args.features, sep='\t', index_col=0)

    k = features.shape[1]
    C = np.array([np.array(x) for x in itertools.product(range(4), repeat=k) if sum(x) <= args.n])
    Ct = C.T

    num = np.dot(features, Ct)
    denom = np.sum(Ct, axis=0)

    paratype_pseudo = 1.0
    S = (paratype_pseudo + num) / (2.0 * paratype_pseudo + denom)

    S_log = np.log(S)
    S_inv = np.log(1 - S)

    with open(args.output, 'w') as outf:
        pickle.dump([S_log, S_inv, features], outf)
