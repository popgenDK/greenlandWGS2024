WORKDIR=PATHTODIR
DATA=PATHTODATA

plink \
    --bfile ${DATA} \
    --maf 0.05 \
    --geno 0.01 \
    --indep-pairwise 1000 kb 1 0.8 \
    --out ${WORKDIR}genos.maf05.geno0

plink \
    --bfile ${DATA} \
    --extract ${WORKDIR}genos.maf05.geno0.prune.in \
    --mind 0.1 \
    --make-bed \
    --out ${WORKDIR}genos.maf05.geno0.pruned

cd ${WORKDIR}
admixture -j80 ${WORKDIR}genos.maf05.geno0.pruned.bed 2 \
    > ${WORKDIR}genos.maf05.geno0.pruned.admixture.log