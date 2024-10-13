#!/usr/bin/env bash

### Step 1 obtain input files and do LD pruning

RAW_PLINK_PREFIX=input.plink
PRUNE_SNP_LIST=input.prune.in
ls ${RAW_PLINK_PREFIX}.bed
ls ${RAW_PLINK_PREFIX}.bim
ls ${RAW_PLINK_PREFIX}.fam

plink --allow-extra-chr --bfile ${RAW_PLINK_PREFIX} --chr 1-22,xy --indep-pairwise 1000 kb 1 0.8 --maf 0.05 --geno 0.05 --out tmp.plink

### Step2 extract pruned data

PLINK_PREFIX=input.plink.maf05.geno05.ldpruned
plink --allow-extra-chr --bfile ${RAW_PLINK_PREFIX} --extract tmp.plink.prune.in --make-bed --out ${PLINK_PREFIX}

### Step2 run admixture with K=2 to obtain population allele frequency and individual ancestry proportion. As described in NGSremix README
### https://github.com/KHanghoj/NGSremix?tab=readme-ov-file#called-genotypes

admixture -j60 ${PLINK_PREFIX}.bed 2 > ${PLINK_PREFIX}.admixture.log

### Step3 run NGSremix to estimate relatedness result
### This step will generate the "ngsremix.res" file, which contains the results required for downstream full sibling or parent-offspring analysis

NGSremix -plink ${PLINK_PREFIX} -f ${PLINK_PREFIX}.2.P -q ${PLINK_PREFIX}.2.Q -P 20 -o ${PLINK_PREFIX}.kin
ls -lh ngsremix.res
