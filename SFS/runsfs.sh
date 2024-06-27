#!/bin/bash
# Assumes an environment with plink1.9, python and the packages: numpy, scipy, and pandas.
PROJSFS=./project_sfs.py

if [[ $# -ne 3 ]]
then
    echo -e "This script need 3 positional arguments:
        1: Path to PLINK-file prefix
        2: Number of alleles to project SFS to
        3: Path to output file"
    exit 1
fi

echo 'Note: Projecting SFS has an upper limit somewhere between 500 and 1245
individuals due to overflow of the binomial coefficient.'

# Name positional arguments
PLINKPRE=$1
PROJTO=$2
OUT=$3

# Counts
plink \
    --bfile ${PLINKPRE} \
    --freq counts \
    --keep-allele-order \
    --out ${OUT}

awk -v projectto=${PROJTO} '$5+$6>=projectto {print $5,$5+$6}' ${OUT}.frq.counts | tail -n +2 | \
    python ${PROJSFS} \
        --projto ${PROJTO} \
        --out ${OUT}

echo 'Number of variants used:'
awk -v projectto=${PROJTO} '$5+$6>=projectto {print $5,$5+$6}' ${OUT}.frq.counts | tail -n +2 | wc -l 
rm ${OUT}.log ${OUT}.nosex ${OUT}.frq.counts