# Activate snakemake environment
source '[CONDAENV]/conda.sh'
mamba activate snakemake

# Go to bin dir
PIPELINEDIR=./GWASpipeline

# Build DAG
snakemake \
    --snakefile ${PIPELINEDIR}gwas_pipeline.smk \
    --configfile ${PIPELINEDIR}gwas_pipeline_config.yaml \
    --dag | dot -Tpdf > GWAS_DAG.pdf

# Run pipeline
snakemake \
    --snakefile ${PIPELINEDIR}gwas_pipeline.smk \
    --configfile ${PIPELINEDIR}gwas_pipeline_config.yaml \
    --cores 100