import pandas as pd
import numpy as np
from tools.procOps import *
from tools.fileOps import *
pd.set_option('chained_assignment',None)


def convert_pileup(pileup_recs, pileup_converter):
    with TemporaryFilePath() as tmp, TemporaryFilePath() as tmp2:
        with open(tmp, 'w') as outf:
            for l in pileup_recs:
                outf.write(l + '\n')
        cmd = ['perl', pileup_converter, tmp, 0, tmp2]
        r = run_proc(cmd, stderr='/dev/null', stdout='/dev/null')
        return load_pileup(tmp2)


def load_pileup(pileup_path):
    r = [x.split() for x in open(pileup_path)]
    return pd.DataFrame(r[1:], columns=r[0])


def ref_count(s):
    return s[s.ref]


bases = {'A', 'T', 'G', 'C'}
def alt_count(s):
    return sum(s[x] for x in bases if x != s.ref)


def make_pileup(bam):
    """Construct a pileup from a bam"""
    cmd = ['samtools', 'mpileup', bam]
    return call_proc_lines(cmd)


def parse_converted_pileup(df, seq):
    df['loc'] = np.array(map(int, df['loc'])) - 1
    df['ref'] = [seq[i] for i in df['loc']]
    df = df[df.ref.isin(bases)]
    df['A'] = pd.to_numeric(df['A']) + pd.to_numeric(df['a'])
    df['C'] = pd.to_numeric(df['C']) + pd.to_numeric(df['c'])
    df['G'] = pd.to_numeric(df['G']) + pd.to_numeric(df['g'])
    df['T'] = pd.to_numeric(df['T']) + pd.to_numeric(df['t'])
    df = df[['loc', 'ref', 'A', 'C', 'G', 'T']]
    df['coverage'] = df[['A', 'T', 'G', 'C']].sum(axis=1)
    df['ref_count'] = df.apply(ref_count, axis=1)
    df['alt_count'] = df.apply(alt_count, axis=1)
    df['ratio'] = 1.0 * df.alt_count / (df.alt_count + df.ref_count)
    return df
