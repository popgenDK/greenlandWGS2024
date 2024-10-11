# Extensions.
PLINKEXTS=['.bed', '.bim', '.fam']
GEMMAGRMEXTS=['.sXX.txt', '.log.txt']
GEMMAASSEXTS=['.assoc.txt', '.log.txt']

# Number of chunks to parallelize GEMMA.
NCHUNKS=config['nchunks']
CHUNKS=[str(x).zfill(2) for x in range(NCHUNKS)]

# Extract gwas name, model and method for gwas-lines in config.
GWASCONFIGS={x['name']: dict(x) for x in config['gwas']}
GWASS=[x for x in GWASCONFIGS.keys()]
MODELS=[GWASCONFIGS[x]['model'] for x in GWASS]
METHODS=[GWASCONFIGS[x]['method'] for x in GWASS]

# Chromosome list for LOCO.
CHRS=[i for i in range(1,22)]

# Input data from config
GENOTYPES=config['genotypes']
PHENOTABLE=config['phenotable']
VARSFORGRM=config['varsforgrm'] # For GRM
GOODVARS=config['goodvars'] # For P-value and manhatten plots
INDVID=config['indvID']

# Dirs.
MODULESDIR="./modules/"
OUTDIR=config['outdir']

# Include methods.
include: f"{MODULESDIR}GEMMA.smk"
include: f"{MODULESDIR}GEMMA_LOCO.smk"

rule all:
    input:
        # Expand the gwas-lines to include the correct wild-cards.
        expand(os.path.join(OUTDIR, '{gwas}/{gwas}.{method}.{model}.assoc.txt.gz'), zip,
            gwas=GWASS, method=METHODS, model=MODELS)

rule prep_pheno:
    input:
        phenotable=PHENOTABLE
    output:
        indvs=os.path.join(OUTDIR, "{gwas}/{gwas}.indvs"),
        pheno=os.path.join(OUTDIR, "{gwas}/{gwas}.pheno"),
        covar=os.path.join(OUTDIR, "{gwas}/{gwas}.covar")
    params:
        phenotype=lambda wildcards: GWASCONFIGS[wildcards.gwas]['phenotype'],
        transformation=lambda wildcards: GWASCONFIGS[wildcards.gwas]['transformation'],
        method=lambda wildcards: GWASCONFIGS[wildcards.gwas]['method'],
        covariates=lambda wildcards: GWASCONFIGS[wildcards.gwas]['covariates'],
        indvid=INDVID
    script:
        f"{MODULESDIR}prep_pheno.R"


rule prep_geno:
    input: 
        geno=multiext(GENOTYPES, *PLINKEXTS),
        indvs=os.path.join(OUTDIR, "{gwas}/{gwas}.indvs")
    output:
        multiext(os.path.join(OUTDIR, "{gwas}/{gwas}"), *PLINKEXTS)
    params:
        outputprefix=os.path.join(OUTDIR, "{gwas}/{gwas}")
    threads: 20
    shell:"""
        plink \
            --bfile {GENOTYPES} \
            --keep {input.indvs} \
            --make-bed \
            --out {params.outputprefix}
    """


rule variantChunks:
    input:
        os.path.join(OUTDIR, "{gwas}/{gwas}.bim")
    output:
        chunks=expand(os.path.join(OUTDIR, "{{gwas}}/variant_chunks/{chunk}.txt"), chunk=CHUNKS)
    params:
        variants=os.path.join(OUTDIR, "{gwas}/{gwas}.vars"),
        chunksdir=directory(os.path.join(OUTDIR, "{gwas}/variant_chunks/"))
    shell: """
        cut -d '\t' -f 2 {input} > {params.variants}
        split -n l/{NCHUNKS} -d --additional-suffix .txt {params.variants} {params.chunksdir}
        rm {params.variants}
    """


rule qq_and_manhatten_plot:
    input:
        assoc=os.path.join(OUTDIR, "{gwas}/{gwas}.{method}.{model}.assoc.txt")
    output:
        qqplot=os.path.join(OUTDIR, "{gwas}/{gwas}.{method}.{model}.assoc.qq_manhatten_plot.png")
    params:
        goodvars=GOODVARS
    script:
        f"{MODULESDIR}qq_and_manhatten_plot.R"


rule cleanup:
    input:
        assoc=os.path.join(OUTDIR, "{gwas}/{gwas}.{method}.{model}.assoc.txt"),
        qqplot=os.path.join(OUTDIR, "{gwas}/{gwas}.{method}.{model}.assoc.qq_manhatten_plot.png")
    output:
        assocgz=os.path.join(OUTDIR, "{gwas}/{gwas}.{method}.{model}.assoc.txt.gz")
    params:
        inputprefix=os.path.join(OUTDIR, "{gwas}/{gwas}")
    shell: """
        gzip {input.assoc}
        rm {params.inputprefix}.bed
        rm {params.inputprefix}.fam
        rm {params.inputprefix}.bim
        rm {params.inputprefix}.nosex
    """
