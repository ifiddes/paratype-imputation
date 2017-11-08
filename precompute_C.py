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
from collections import *
from phase_lib import *


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', help="Maximum number of haplotypes", default=10, type=int)
    parser.add_argument('--features', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--paratype-pseudo', default=0.1, type=float)
    parser.add_argument('--inferred-copy-numbers', nargs=4, default=['4', '2', '2', '2'])
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
    features = pd.read_csv(args.features, sep='\t', index_col=0)

    k = features.shape[1]
    C = construct_C(args.inferred_copy_numbers, features)
    with open(args.output, 'w') as outf:
        pickle.dump(C, outf)
