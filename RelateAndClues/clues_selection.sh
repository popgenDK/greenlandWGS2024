# Input arguments
if [[ $# -ne 8 ]]
then
echo -e "This script need 8 positional arguments:
    1: Prefix to .anc and .mut files
    2: Position of variant
    3: Derived allele frequency
    4: Path to coalescence file
    5: Number of trees to sample
    6: Thin
    7: Dominance (0=recessive, 0.5=additive, 1=dominant)
    8: Prefix to outputfiles"
exit 1
fi

# Name positional arguments
ANCMUTPRE=$1
POS=$2
DAF=$3
COAL=$4
NTREES=$5
THIN=$6
DOM=$7
OUTPRE=$8
echo 'Input parameters'
echo 'ANCMUTPRE: '${ANCMUTPRE}
echo 'POS: '${POS}
echo 'DAF: '${DAF}
echo 'COAL: '${COAL}
echo 'NTREES: '${NTREES}
echo 'THIN: '${THIN}
echo 'DOM: '${DOM}
echo 'OUTPRE: '${OUTPRE}


# PATHS
RELATE=PATHTORELATE/relate_v1.2.1_x86_64_static/
SAMPLEBRANCH=${RELATE}scripts/SampleBranchLengths/SampleBranchLengths.sh
CLUES=PATHTOCLUES/inference.py

# Environment for CLUES.
echo "Activating mamba environment for Clues"
source '[CONDAENV]/conda.sh'
mamba activate clues
cd PATHTOCLUES

# Sample trees
bash ${SAMPLEBRANCH} \
    -i ${ANCMUTPRE} \
    -o ${OUTPRE} \
    -m 2.5e-8 \
    --coal ${COAL} \
    --format b \
    --num_samples ${NTREES} \
    --first_bp ${POS} \
    --last_bp ${POS} \
    --seed 1

# Cleanup
echo 'Removing .anc, .mut, and .dist intermediary files'
rm ${OUTPRE}.anc ${OUTPRE}.mut ${OUTPRE}.dist

echo 'Running Clues (very slow)...'
python ${CLUES} \
    --times ${OUTPRE} \
    --out ${OUTPRE}.clues \
    --coal ${COAL} \
    --thin ${THIN} \
    --dom ${DOM} \
    --popFreq ${DAF} > ${OUTPRE}.clues.log

echo 'Removing Clues np arrays .npy'
rm ${OUTPRE}.clues.epochs.npy ${OUTPRE}.clues.freqs.npy ${OUTPRE}.clues.post.npy

VARNAME=$(basename "$OUTPRE")
paste <(echo ${VARNAME}) \
    <(echo ${DAF}) \
    <(echo ${DOM}) \
    <(grep 0-1000 ${OUTPRE}.clues.log | cut -f 2) \
    <(grep logLR ${OUTPRE}.clues.log | cut -d ' ' -f 2) > ${OUTPRE}.clues.tsv
rm ${OUTPRE}.timeb ${OUTPRE}.clues.log
echo 'Output format: VariantName DerivedAF DominanceCoef SelectCoefS logLR'
exit
