#!/bin/bash

workdir=/maps/projects/greenland-AUDIT/people/rlk420/greenland/ld
assocdir=/maps/projects/greenland-AUDIT/people/mqr375/GWAS22/output/gwasres
datadir=/maps/projects/greenland-AUDIT/people/mqr375/GWAS22/intermediate/postgwas
phenofile=$workdir/data/phenos.txt
run=$1 ## add or rev
outprefix=gwasres_LOCO ## add or rev

# allphenos=$(for name in `grep $run $phenofile | cut -f1`;do
#                 assoc="$assocdir/$outprefix/$name/$name.gemmaloco.$run.assoc.txt"
#                 if [ ! -f $assoc ];then echo "$assoc not exist " && exit ; fi
#                 echo $assoc
#             done | paste -s -d ',')
# echo "gtime -v PCAone -B results/$run-pca/k5.svd2.residuals --clump $allphenos --clump-r2 0.001 --clump-p1 1e-6 --clump-p2 1e-6 --clump-bp 10000000 --clump-names chr,ps,p_score -n 22 -o $workdir/results/$outprefix/$run" 

# exit

if [[ $run == "rec" ]];then
    ## run PCAone on goodvar.rec data with k=5
    PCAone -b $datadir/goodvars.rec -k 5 --ld-r2 1 -n 40 -o results/pca/k5.svd2
    ## do rec first
    for name in `awk '/rec/{print $1}' $phenofile`;do
        assoc="$assocdir/$name/$name.gemma.rec.assoc.txt"
        if [ ! -s $assoc ];then echo $assoc  && exit ; fi
        outdir=$workdir/results/$name
        mkdir -p $outdir
        echo "gtime -v PCAone -B results/rec-pca/k5.svd2.residuals --clump $assoc --clump-r2 0.001 --clump-p1 1e-6 --clump-p2 1e-6 --clump-bp 10000000 --clump-names chr,ps,p_score -n 22 -o $outdir/$name.rec" 
    done | parallel -j 1 -k "sh -c {}"
fi

## do add now
if [[ $run == "add" ]];then
    ## run PCAone on goodvar.add data with k=5
    PCAone -b $datadir/goodvars.add -k 5 --ld-r2 1 -n 40 -o results/pca/k5.svd2

    for name in `awk '/add/{print $1}' $phenofile`;do
        assoc="$assocdir/$name/$name.gemma.add.assoc.txt"
        if [ ! -s $assoc ];then echo $assoc && exit ; fi
        outdir=$workdir/results/$name
        mkdir -p $outdir
        echo "gtime -v PCAone -B results/add-pca/k5.svd2.residuals --clump $assoc --clump-r2 0.001 --clump-p1 1e-6 --clump-p2 1e-6 --clump-bp 10000000 --clump-names chr,ps,p_score -n 22 -o $outdir/$name.add" 
    done | parallel -j 1 -k "sh -c {}"
fi

