GEMMAADD="./gemmaADD"
GEMMAREC="./gemmaREC"


rule prep_LOCO_varfiles:
    input:
        goodSNPs=VARSFORGRM,
        allVarsBIM=GENOTYPES+'.bim'
    output:
        expand(os.path.join(OUTDIR, "snpsforGRM_LOCO/LOCO_chr{CHR}.txt"), CHR=CHRS),
        expand(os.path.join(OUTDIR, "chrvars/chr{CHR}.txt"), CHR=CHRS)
    params:
        prefixloco=os.path.join(OUTDIR, "snpsforGRM_LOCO/LOCO_chr"),
        prefixchroms=os.path.join(OUTDIR, "chrvars/chr")
    shell: """
    for i in {{1..22}}
    do
        grep -v chr${{i}}: {input.goodSNPs} > {params.prefixloco}${{i}}.txt
        grep chr${{i}}: {input.allVarsBIM} | cut -f 2 > {params.prefixchroms}${{i}}.txt
    done
    """


rule gemmaloco_GRMs:
    input:
        geno=multiext(os.path.join(OUTDIR, "{gwas}/{gwas}"), *PLINKEXTS),
        pheno=os.path.join(OUTDIR, "{gwas}/{gwas}.pheno"),
        covar=os.path.join(OUTDIR, "{gwas}/{gwas}.covar"),
        snpsforGRM_LOCO=expand(os.path.join(OUTDIR, "snpsforGRM_LOCO/LOCO_chr{CHR}.txt"), CHR=CHRS)
    output:
        expand(os.path.join(OUTDIR, "{{gwas}}/GRMs/{{gwas}}.gemmaloco.LOCO_chr{CHR}{ext}"), CHR=CHRS, ext=GEMMAGRMEXTS)
    params:
        gwasprefix=os.path.join(OUTDIR, "{gwas}/"),
        grmsdir=os.path.join(OUTDIR, "{gwas}/GRMs/"),
        locosnps=os.path.join(OUTDIR, "snpsforGRM_LOCO/LOCO_chr")
    threads: 16
    shell:"""
    parallel --jobs 4 \
        {GEMMAADD} \
            -bfile {params.gwasprefix}{wildcards.gwas} \
            -gk 2 \
            -snps {{}} \
            -maf 0 \
            -miss 1 \
            -r2 1 \
            -p {input.pheno} \
            -c {input.covar} \
            -outdir {params.grmsdir} \
            -o {wildcards.gwas}.gemmaloco.{{/.}} ::: {params.locosnps}*.txt
    """


rule gemmalocoADD_parallelized:
    input:
        geno=multiext(os.path.join(OUTDIR, "{gwas}/{gwas}"), *PLINKEXTS),
        grms=expand(os.path.join(OUTDIR, "{{gwas}}/GRMs/{{gwas}}.gemmaloco.LOCO_chr{CHR}.sXX.txt"), CHR=CHRS),
        pheno=os.path.join(OUTDIR, "{gwas}/{gwas}.pheno"),
        covar=os.path.join(OUTDIR, "{gwas}/{gwas}.covar"),
        chrvars=expand(os.path.join(OUTDIR, "chrvars/chr{CHR}.txt"), CHR=CHRS)
    output:
        expand(os.path.join(OUTDIR, "{{gwas}}/GRMs/{{gwas}}.gemmaloco.add.LOCO_chr{CHR}{ext}"), CHR=CHRS, ext=GEMMAASSEXTS),
        os.path.join(OUTDIR, "{gwas}/GRMs/{gwas}.gemmaloco.add.log.txt")
    params:
        inputprefix=os.path.join(OUTDIR, "{gwas}/{gwas}"),
        chrvarsprefix=os.path.join(OUTDIR, "chrvars/"),
        outputdir=os.path.join(OUTDIR, "{gwas}/GRMs/"),
        grmprefix=os.path.join(OUTDIR, "{gwas}/GRMs/{gwas}.gemmaloco.LOCO_")
    threads: 16
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
            -o {wildcards.gwas}.gemmaloco.add.LOCO_{{/.}} \
            -snps {{}} \
            -k {params.grmprefix}{{/.}}.sXX.txt ::: {params.chrvarsprefix}*.txt > {params.outputdir}{wildcards.gwas}.gemmaloco.add.log.txt
    """


rule gemmalocoREC_parallelized:
    input:
        geno=multiext(os.path.join(OUTDIR, "{gwas}/{gwas}"), *PLINKEXTS),
        grms=expand(os.path.join(OUTDIR, "{{gwas}}/GRMs/{{gwas}}.gemmaloco.LOCO_chr{CHR}.sXX.txt"), CHR=CHRS),
        pheno=os.path.join(OUTDIR, "{gwas}/{gwas}.pheno"),
        covar=os.path.join(OUTDIR, "{gwas}/{gwas}.covar"),
        chrvars=expand(os.path.join(OUTDIR, "chrvars/chr{CHR}.txt"), CHR=CHRS)
    output:
        expand(os.path.join(OUTDIR, "{{gwas}}/GRMs/{{gwas}}.gemmaloco.rec.LOCO_chr{CHR}{ext}"), CHR=CHRS, ext=GEMMAASSEXTS),
        os.path.join(OUTDIR, "{gwas}/GRMs/{gwas}.gemmaloco.rec.log.txt")
    params:
        inputprefix=os.path.join(OUTDIR, "{gwas}/{gwas}"),
        chrvarsprefix=os.path.join(OUTDIR, "chrvars/"),
        outputdir=os.path.join(OUTDIR, "{gwas}/GRMs/"),
        grmprefix=os.path.join(OUTDIR, "{gwas}/GRMs/{gwas}.gemmaloco.LOCO_")
    threads: 16
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
            -o {wildcards.gwas}.gemmaloco.rec.LOCO_{{/.}} \
            -snps {{}} \
            -k {params.grmprefix}{{/.}}.sXX.txt ::: {params.chrvarsprefix}*.txt > {params.outputdir}{wildcards.gwas}.gemmaloco.rec.log.txt
    """


rule collapseGemmaloco:
    input:
        expand(os.path.join(OUTDIR, "{{gwas}}/GRMs/{{gwas}}.gemmaloco.{{model}}.LOCO_chr{CHR}{ext}"), CHR=CHRS, ext=GEMMAASSEXTS)
    output:
        assoc=os.path.join(OUTDIR, "{gwas}/{gwas}.gemmaloco.{model}.assoc.txt"),
        log=os.path.join(OUTDIR, "{gwas}/{gwas}.gemmaloco.{model}.chrslogs.txt")
    params:
        chrsprefix=os.path.join(OUTDIR, "{gwas}/GRMs/{gwas}.gemmaloco.{model}.LOCO_chr"),
        grmprefix=os.path.join(OUTDIR, "{gwas}/GRMs/{gwas}.gemmaloco.LOCO_chr")
    shell:"""
        head -n 1 {params.chrsprefix}1.assoc.txt > {output.assoc}
        touch {output.log}
        for i in {{1..22}}
        do
            tail -n +2 -q {params.chrsprefix}${{i}}.assoc.txt >> {output.assoc}
            cat {params.chrsprefix}${{i}}.log.txt >> {output.log}
            gzip {params.grmprefix}${{i}}.sXX.txt
        done
        rm {params.chrsprefix}*.assoc.txt
        rm {params.chrsprefix}*.log.txt
    """