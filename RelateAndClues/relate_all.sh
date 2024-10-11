# Relatepaths
RELATE=PATHTORELATE/relate_v1.2.1_x86_64_static/
RELATEFORMATS=${RELATE}bin/RelateFileFormats
PREPINPUT=${RELATE}scripts/PrepareInputFiles/PrepareInputFiles.sh
RELATETHREADED=${RELATE}scripts/RelateParallel/RelateParallel.sh
RELATEALL=${RELATE}bin/Relate
ESTPOPSIZE=${RELATE}scripts/EstimatePopulationSize/EstimatePopulationSize.sh
SAMPLEBRANCH=${RELATE}scripts/SampleBranchLengths/SampleBranchLengths.sh

WORKDIR=PATHTOWORKINGDIRECTORY
DATA=${WORKDIR}data.bcf
OUTDIR=${WORKDIR}output/

# From here: http://ftp.ensembl.org/pub/current_fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh38.tar.gz
ANC=${WORKDIR}input/homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_

parallel bcftools view ${DATA} \
    -r chr{} \
    -O z \
    --threads 32 \
    -c 1 \
    -Q 0.9999 \
    -v snps \
    -o ${OUTDIR}vcf/chr{}.vcf.gz ::: {1..22}

nohup parallel $RELATEFORMATS \
    --mode ConvertFromVcf \
    --haps ${OUTDIR}hapssample/chr{}.haps \
    --sample ${OUTDIR}hapssample/chr{}.sample \
    -i ${OUTDIR}vcf/chr{} ::: {1..22} > ${OUTDIR}hapssample/convert_from_vcf.log

nohup parallel $PREPINPUT \
    --haps ${OUTDIR}hapssample/chr{}.haps \
    --sample ${OUTDIR}hapssample/chr{}.sample \
    --ancestor ${ANC}{}.highconf.fa \
    -o ${OUTDIR}prepped/chr{}_prepped \> ${OUTDIR}prepped/chr{}_prepped.log \
        ::: {1..22} > ${OUTDIR}prepped/prep_input.log

cd ${OUTDIR}ancmut

for i in {1..22}
do
    nohup nice ${RELATETHREADED} \
        --mode All \
        --seed 1 \
        --threads 100 \
        -m 2.5e-8 \
        -N 2000 \
        --haps ${OUTDIR}prepped/chr${i}_prepped.haps.gz \
        --sample ${OUTDIR}prepped/chr${i}_prepped.sample.gz \
        --map ${WORKDIR}input/genetic_map_GRCh38/genetic_map_GRCh38_chr${i}.txt \
        -o allindvs_chr${i} > allindvs_chr${i}.log
done

echo sample population group sex > ${OUTDIR}allindvs.poplabels
tail -n +3 ${OUTDIR}hapssample/chr19.sample | awk '{print $1,"GRL","GRL","NA"}' >> ${OUTDIR}allindvs.poplabels

mkdir ${OUTDIR}popsize
nohup nice ${ESTPOPSIZE} \
    -i ${OUTDIR}ancmut/allindvs \
    -m 2.5e-8 \
    --years_per_gen 28 \
    --poplabels ${OUTDIR}allindvs.poplabels \
    --threads 100 \
    --first_chr 1 \
    --last_chr 22 \
    -o ${OUTDIR}popsize/allindvs.popsize > ${OUTDIR}popsize/allindvs.popsize.log


mkdir ${OUTDIR}samplevars
parallel --link echo bash ${SAMPLEBRANCH} \
    -i ${OUTDIR}popsize/allindvs.popsize_chr{2} \
    -o ${OUTDIR}samplevars/allindvs.{1}.sample.n \
    -m 2.5e-8 \
    --coal ${OUTDIR}popsize/allindvs.popsize.coal \
    --format n \
    --num_samples 10000 \
    --first_bp {3} \
    --last_bp {3} \
    --threads 10 \
    --seed 1 \
        ::: CPT1A FADS \
        ::: 11 11 \
        ::: 68780662 61860488
