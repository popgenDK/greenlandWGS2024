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

# Prior to this, the genome has been split into 6MB regions and 150 individuals
# with inferred Inuit genetic ancestry on both haplotypes in the region randomly
# sampled.

# NCHUNKS=460
mkdir ${OUTDIR}vcf
parallel -j 10 bcftools view ${DATA} \
    -R ${OUTDIR}chunk_samples/chunks_chr{}.region \
    -S ${OUTDIR}chunk_samples/chunks_chr{}.samples \
    -c 1 \
    -Q 0.9999 \
    -v snps \
    --threads 32 \
    -O z \
    -o ${OUTDIR}vcf/chunks_chr{}.oriid.vcf.gz ::: {1..460}
parallel -j 10 bcftools reheader \
    --samples ${OUTDIR}chunk_samples/chunks_chr{}.rename \
    -o ${OUTDIR}vcf/chunks_chr{}.vcf.gz \
    ${OUTDIR}vcf/chunks_chr{}.oriid.vcf.gz ::: {1..460}
parallel rm ${OUTDIR}vcf/chunks_chr{}.oriid.vcf.gz ::: {1..460}

mkdir ${OUTDIR}hapssample
nohup parallel -j 10 $RELATEFORMATS \
    --mode ConvertFromVcf \
    --haps ${OUTDIR}hapssample/chunks_chr{}.haps \
    --sample ${OUTDIR}hapssample/chunks_chr{}.sample \
    -i ${OUTDIR}vcf/chunks_chr{} ::: {1..460} > ${OUTDIR}hapssample/convert_from_vcf.log


mkdir ${OUTDIR}prepped
for i in {1..22}
do
    nohup parallel -j 15 $PREPINPUT \
        --haps ${OUTDIR}hapssample/chunks_chr{}.haps \
        --sample ${OUTDIR}hapssample/chunks_chr{}.sample \
        --ancestor ${ANC}${i}.highconf.fa \
        -o ${OUTDIR}prepped/chunks_chr{}_prepped \
        :::: ${OUTDIR}chunk_size/chr${i}_chunks_6MB.txt > ${OUTDIR}prepped/prep_input_chr${i}.log
done

mkdir ${OUTDIR}ancmut
cd ${OUTDIR}ancmut
for i in {1..22}
do
    nohup parallel nice ${RELATEALL} \
        --mode All \
        --seed 1 \
        -m 2.5e-8 \
        -N 2000 \
        --haps ${OUTDIR}prepped/chunks_chr{}_prepped.haps.gz \
        --sample ${OUTDIR}prepped/chunks_chr{}_prepped.sample.gz \
        --map ${WORKDIR}input/genetic_map_GRCh38/genetic_map_GRCh38_chr${i}.txt \
        -o mosaic_chunks_chr{} \
        :::: ${OUTDIR}chunk_size/chr${i}_chunks_6MB.txt > chr${i}_chunks.log
    parallel rm -r mosaic_chunks_chr{} :::: ${OUTDIR}chunk_size/chr${i}_chunks_6MB.txt
done
echo sample population group sex > ${OUTDIR}mosaic.poplabels
tail -n +3 ${OUTDIR}hapssample/chunks_chr1.sample | awk '{print $1,"GRL","GRL","NA"}' >> ${OUTDIR}mosaic.poplabels

mkdir ${OUTDIR}popsize
nohup nice ${ESTPOPSIZE} \
    -i ${OUTDIR}ancmut/mosaic_chunks \
    -m 2.5e-8 \
    --years_per_gen 28 \
    --poplabels ${OUTDIR}mosaic.poplabels \
    --first_chr 1 \
    --last_chr 460 \
    -o ${OUTDIR}popsize/mosaic_chunks.popsize > ${OUTDIR}popsize/mosaic_chunks.popsize.log

# Sample 10,000 local trees for each variant of interest.
mkdir ${OUTDIR}samplevars
parallel --link bash ${SAMPLEBRANCH} \
    -i ${OUTDIR}popsize/mosaic_chunks.popsize_chr{2} \
    -o ${OUTDIR}samplevars/mosaic.{1}.sample.n \
    -m 2.5e-8 \
    --coal ${OUTDIR}popsize/mosaic_chunks.popsize.coal \
    --format n \
    --num_samples 10000 \
    --first_bp {3} \
    --last_bp {3} \
    --threads 10 \
    --seed 1 \
        ::: ATP8B1 PCCA_B LDLR TBC1D4 ACDY3 HNF1A SI \
        ::: 422 102 428 351 43 340 107 \
        ::: 57674993 136329943 11105315 75324385 24827609 120996541 165069176
