#!/bin/bash

GLVCF="/maps/projects/greenland-AUDIT/scratch/zilong/imputation/mega22/results/merge/all.merge_refpanel1_phased_impute2.vcf.gz"
GLVCF="/maps/projects/greenland-AUDIT/people/rlk420/greenland/imputation_order/results/merge/reheader/all.merge_refpanel1_phased_impute2.vcf.gz"
OUT="/maps/projects/greenland-AUDIT/people/rlk420/greenland/rfmix/results"
EUSAM="/maps/projects/greenland-AUDIT/people/rlk420/greenland/rfmix/data/EU.sample"
KGP="/maps/projects/alab/scratch/zilong/1000G.hg38/vcf"
GENMAP=/maps/projects/greenland-AUDIT/people/rlk420/greenland/rfmix/data/genetic_map_hg38.txt
RFMIX=/home/rlk420/.local/bin/rfmix

for chr in chr{22..1};do
    
outdir=$OUT/GL
mkdir -p $outdir
outvcf=$outdir/$chr.vcf.gz
## remove sites with genotype missing rate < 0.1
bcftools view -i 'F_MISSING < 0.1' -Oz -o $outvcf --threads 4 $GLVCF $chr && bcftools index -f $outvcf

outdir=$OUT/EU
mkdir -p $outdir
vcf=$KGP/CCDG_14151_B01_GRM_WGS_2020-08-05_${chr}.filtered.shapeit2-duohmm-phased.vcf.gz
euvcf=$outdir/$chr.vcf.gz

if [ ! -f $euvcf ];then
    bcftools view -S $EUSAM --threads 4 -Oz -o $euvcf $vcf && bcftools index -f $euvcf
fi

outdir="$OUT/merge"
mkdir -p $outdir
euvcf=$OUT/EU/$chr.vcf.gz
glvcf=$OUT/GL/$chr.vcf.gz
mergevcf=$outdir/$chr.vcf.gz

if [ ! -f $mergevcf.csi ];then
    tmpvcf=$outdir/tmp.$chr.vcf.gz
    bcftools merge $euvcf $glvcf --threads 4 -Oz -o $tmpvcf && bcftools index -f $tmpvcf
    bcftools view -i 'F_MISSING < 0.1 & MAF > 0.05' -Oz -o $mergevcf --threads 4 $tmpvcf && bcftools index -f $mergevcf
fi

REFMAP=$OUT/meta/ref.sample.map
ALLSAM=$OUT/meta/all.sample
REFSAM=$OUT/meta/ref.sample
ADMIXSAM=$OUT/meta/admix.sample

# bcftools query -l $mergevcf > $ALLSAM
# cut -f1 $REFMAP > $REFSAM
# # preserve order?
# rg -v -w -f $REFSAM $ALLSAM > $ADMIXSAM

outdir="$OUT/merge"
mkdir -p $outdir
refvcf=$outdir/ref.$chr.vcf.gz
admixvcf=$outdir/admix.$chr.vcf.gz
mergevcf=$outdir/$chr.vcf.gz
bcftools view -S $REFSAM -Oz -o $refvcf --threads 10 $mergevcf && bcftools index -f $refvcf 
bcftools view -S $ADMIXSAM -Oz -o $admixvcf --threads 10 $mergevcf && bcftools index -f $admixvcf

## run rfmix now
outdir=$OUT/rfmix
mkdir -p $outdir
out=$outdir/$chr
$RFMIX -f $admixvcf -r $refvcf -m $REFMAP -g $GENMAP --chromosome=$chr -o $out &> $out.llog

done
