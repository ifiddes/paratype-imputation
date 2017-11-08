from tools.procOps import *


if __name__ == '__main__':
    # run tests
    cmd = ['python', 'full_genotype.py', '--features', 'very_reduced.features.tsv', '--name', 'H9',
           '--inferred-copy-numbers', '4 2 2 2', '--pileup', 'H9.pileup.txt', '--deviance-plot', 'H9_deviance.pdf']
    run_proc(cmd)
    cmd = ['python', 'full_genotype.py', '--features', 'very_reduced.features.tsv', '--name', 'NA12878',
           '--inferred-copy-numbers', '4 2 2 2', '--pileup', 'NA19240.pileup.txt', '--deviance-plot',  'NA12878_deviance.pdf']
    run_proc(cmd)
    cmd = ['python', 'full_genotype.py', '--features', 'very_reduced.features.tsv', '--name', 'NA19240',
           '--inferred-copy-numbers', '4 2 2 2', '--pileup', 'NA12878.pileup.txt', '--deviance-plot', 'NA19240_deviance.pdf']
    run_proc(cmd)
    cmd = ['python', 'full_genotype.py', '--features', 'very_reduced.features.tsv', '--name', 'LP6005441',
           '--inferred-copy-numbers', '4 2 1 2', '--pileup', '/hive/users/cbosworth/SimonsNormals/secondTry/LP6005441-DNA_H06.parsed.pileup',
           '--deviance-plot', 'LP6005441_deviance.pdf']
    run_proc(cmd)