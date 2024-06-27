import numpy as np
import argparse
import sys
import pandas as pd
from scipy.special import binom

def proj_sfs(n, m, k, j):
    """
    Project SFS to m alleles.
    n: total alleles
    m: projection down to m alleles
    k: obs derived allele OR obs minor allele
    j: SFS bins, usually array from 0 to m+1
    """
    if m > n:
        sys.stderr.write(f"m: {m} is larger than n: {n}. EXITING\n.")
        exit(1)
    return (binom(m, j) * binom(n-m, k-j)) /  binom(n, k)

parser = argparse.ArgumentParser(description="""Project SFS.\nAssumes a stdin
    stream of lines where each line represent a variant with two colummns
    separated by a space. First column in each line is the number of derived
    alleles for that variant. Second column in each line is the total number of
    alleles for that variant""")
parser.add_argument('--projto', type=int, help="""Number of alleles to project
    down to. Will cause error and exit if larger than number of total alleles at 
    any point.""")
parser.add_argument('--output', type=str, help="Path to output sfs.")
args = parser.parse_args()
p = args.projto
o = args.output

sfsbins = np.array([i for i in range(p+1)])
sfs = np.zeros(sfsbins.shape[0])

for line in sys.stdin:
    derived, ntot = line.rstrip().split(' ')
    sfs += proj_sfs(int(ntot), p, int(derived), sfsbins)

df = pd.DataFrame({'nderived': sfsbins, 'nsnps': sfs})
df.to_csv(args.output, sep='\t', index=False)