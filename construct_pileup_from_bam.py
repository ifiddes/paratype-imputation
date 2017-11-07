from __future__ import division
import numpy as np
import argparse
import vcf
import itertools
import sys
import multiprocessing
import seaborn as sns
import pandas as pd
from tools.bio import *
from phase_lib import *


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', help="Maximum number of haplotypes", default=10, type=int)
    parser.add_argument('--bam', help='input BAM')
    parser.add_argument('--output', help='Output parsed pileup')
    parser.add_argument('--consensus-fasta', default='/hive/users/cbosworth/refs/notch/notch2_aligned_consensus.fasta')
    parser.add_argument('--pileup-converter', default='/cluster/home/ifiddes/pileup2base/pileup2base.pl')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    pileup_recs = make_pileup(args.bam)
    df = convert_pileup(pileup_recs, args.pileup_converter)
    _, seq = read_fasta(args.consensus_fasta, None).next()
    data = parse_converted_pileup(df, seq)
    data.to_csv(args.output, sep='\t')
