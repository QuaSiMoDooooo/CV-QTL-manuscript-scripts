#! python

# -------------
# FileName     : 04_ld_snakemake.smk
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Snakemake workflow for LD (Linkage Disequilibrium) analysis using PLINK
# -------------

rule all:
    input:
        expand("04_ld_result/chr{chr}.ld.gz", chr=range(1, 23))

rule create_directory:
    output:
        directory("04_ld_result")
    shell:
        "mkdir -p {output}"

rule run_plink:
    input:
        "03_all_vartype_merged.vcf"
    output:
        ld="04_ld_result/chr{chr}.ld.gz"
    params:
        chr="{chr}"
    shell:
        """
        plink \
        --vcf {input} \
        --double-id \
        --allow-extra-chr \
        --r2 gz \
        --ld-window 100 \
        --ld-window-kb 1000 \
        --ld-window-r2 0 \
        --make-bed \
        --chr {params.chr} \
        --out 04_ld_result/chr{params.chr}
        """