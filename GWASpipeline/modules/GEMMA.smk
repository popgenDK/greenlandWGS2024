GEMMAADD="./gemmaADD"
GEMMAREC="./gemmaREC"

rule gemmaGRM:
    input:
        geno=multiext(os.path.join(OUTDIR, "{gwas}/{gwas}"), *PLINKEXTS),
        pheno=os.path.join(OUTDIR, "{gwas}/{gwas}.pheno"),
        covar=os.path.join(OUTDIR, "{gwas}/{gwas}.covar")
    output:
        multiext(os.path.join(OUTDIR, "{gwas}/{gwas}.gemma"), *GEMMAGRMEXTS),
    params:
        gwasprefix=os.path.join(OUTDIR, "{gwas}/")
    threads: 4
    shell:"""
    {GEMMAADD} \
        -bfile {params.gwasprefix}{wildcards.gwas} \
        -gk 2 \
        -snps {VARSFORGRM} \
        -maf 0 \
        -miss 1 \
        -r2 1 \
        -p {input.pheno} \
        -c {input.covar} \
        -outdir {params.gwasprefix} \
        -o {wildcards.gwas}.gemma
    """

rule gemmaADD_parallelized:
    input:
        geno=multiext(os.path.join(OUTDIR, "{gwas}/{gwas}"), *PLINKEXTS),
        grm=os.path.join(OUTDIR, "{gwas}/{gwas}.gemma.sXX.txt"),
        pheno=os.path.join(OUTDIR, "{gwas}/{gwas}.pheno"),
        covar=os.path.join(OUTDIR, "{gwas}/{gwas}.covar"),
        chunks=expand(os.path.join(OUTDIR, "{{gwas}}/variant_chunks/{chunk}.txt"), chunk=CHUNKS)
    output:
        expand(os.path.join(OUTDIR, "{{gwas}}/{{gwas}}.gemma.add.{chunk}{ext}"), chunk=CHUNKS, ext=GEMMAASSEXTS),
        os.path.join(OUTDIR, "{gwas}/{gwas}.gemma.add.log.txt")
    params:
        inputprefix=os.path.join(OUTDIR, "{gwas}/{gwas}"),
        outputdir=os.path.join(OUTDIR, "{gwas}/"),
        chunksdir=os.path.join(OUTDIR, "{gwas}/variant_chunks/")
    threads: 8
    shell:"""
    parallel --jobs 4 \
        {GEMMAADD} \
            -bfile {params.inputprefix} \
            -p {input.pheno} \
            -c {input.covar} \
            -lmm 3 \
            -maf 0 \
            -miss 1 \
            -r2 1 \
            -outdir {params.outputdir} \
            -o {wildcards.gwas}.gemma.add.{{/.}} \
            -snps {{}} \
            -k {input.grm} ::: {params.chunksdir}*.txt > {params.outputdir}{wildcards.gwas}.gemma.add.log.txt
    rm -r {params.chunksdir}
    """


rule gemmaREC_parallelized:
    input:
        geno=multiext(os.path.join(OUTDIR, "{gwas}/{gwas}"), *PLINKEXTS),
        grm=os.path.join(OUTDIR, "{gwas}/{gwas}.gemma.sXX.txt"),
        pheno=os.path.join(OUTDIR, "{gwas}/{gwas}.pheno"),
        covar=os.path.join(OUTDIR, "{gwas}/{gwas}.covar"),
        chunks=expand(os.path.join(OUTDIR, "{{gwas}}/variant_chunks/{chunk}.txt"), chunk=CHUNKS)
    output:
        expand(os.path.join(OUTDIR, "{{gwas}}/{{gwas}}.gemma.rec.{chunk}{ext}"), chunk=CHUNKS, ext=GEMMAASSEXTS),
        os.path.join(OUTDIR, "{gwas}/{gwas}.gemma.rec.log.txt")
    params:
        inputprefix=os.path.join(OUTDIR, "{gwas}/{gwas}"),
        outputdir=os.path.join(OUTDIR, "{gwas}/"),
        chunksdir=os.path.join(OUTDIR, "{gwas}/variant_chunks/")
    threads: 8
    shell:"""
    parallel --jobs 4 \
        {GEMMAREC} \
            -bfile {params.inputprefix} \
            -p {input.pheno} \
            -c {input.covar} \
            -lmm 3 \
            -maf 0 \
            -miss 1 \
            -r2 1 \
            -outdir {params.outputdir} \
            -o {wildcards.gwas}.gemma.rec.{{/.}} \
            -snps {{}} \
            -k {input.grm} ::: {params.chunksdir}*.txt > {params.outputdir}{wildcards.gwas}.gemma.rec.log.txt
    rm -r {params.chunksdir}
    """


rule collapseGemma:
    input:
        expand(os.path.join(OUTDIR, "{{gwas}}/{{gwas}}.gemma.{{model}}.{chunk}{ext}"), chunk=CHUNKS, ext=GEMMAASSEXTS)
    output:
        assoc=os.path.join(OUTDIR, "{gwas}/{gwas}.gemma.{model}.assoc.txt"),
        log=os.path.join(OUTDIR, "{gwas}/{gwas}.gemma.{model}.chunklogs.txt")
    params:
        chuncksprefix=os.path.join(OUTDIR, "{gwas}/{gwas}.gemma.{model}")
    shell:"""
        head -n 1 {params.chuncksprefix}.00.assoc.txt > {output.assoc}
        tail -n +2 -q {params.chuncksprefix}.*.assoc.txt >> {output.assoc}
        cat {params.chuncksprefix}.*.log.txt > {output.log}
        rm {params.chuncksprefix}.*.assoc.txt
        rm {params.chuncksprefix}.*.log.txt
    """