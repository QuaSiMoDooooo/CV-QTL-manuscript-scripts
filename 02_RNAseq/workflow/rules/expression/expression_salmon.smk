# -*- coding: UTF-8 -*-
#
# FileName     : salmon
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-02-08 20:39
# Last Modified: 2024-02-08 20:39
# Modified By  : EastsunW
# -------------
# Description  :
# -------------


rule expression:
    input:
        seq_1=lambda wildcards: f"{cohort.input_dir}/{cohort.list_files(wildcards.sample)[0]}",
        seq_2=lambda wildcards: f"{cohort.input_dir}/{cohort.list_files(wildcards.sample)[1]}",
    output:
        gene="results/{cohort}/expression/expression_sample/{sample}.gene.txt",
        transc="results/{cohort}/expression/expression_sample/{sample}.transcript.txt",
        tempdir=temp(
            directory("results/{cohort}/expression/expression_sample/{sample}")
        ),
        infojson="results/{cohort}/expression/expression_sample/{sample}.lib_format_counts.json",
    log:
        "logs/{cohort}/expression/{sample}.salmon.log",
    benchmark:
        "benchmarks/{cohort}/expression/{sample}.salmon.benchmark"
    threads: config["expression"]["salmon"]["threads"]
    params:
        salmon_idx=cohort.ref_config["expression"]["salmon"]["index"],
        gtf=cohort.ref_config["expression"]["salmon"]["gtf"],
        extra=config["expression"]["salmon"]["extra_params"],
    shell:
        """
        (
            salmon quant \
                {params.extra} \
                --index {params.salmon_idx} \
                --geneMap {params.gtf} \
                -1 {input.seq_1} \
                -2 {input.seq_2} \
                --threads {threads} \
                --output results/{wildcards.cohort}/expression/expression_sample/{wildcards.sample}
            mv results/{wildcards.cohort}/expression/expression_sample/{wildcards.sample}/quant.genes.sf {output.gene}
            mv results/{wildcards.cohort}/expression/expression_sample/{wildcards.sample}/quant.sf {output.transc}
            mv results/{wildcards.cohort}/expression/expression_sample/{wildcards.sample}/lib_format_counts.json {output.infojson}
        ) &> {log}
        """


rule gene_id_mapping:
    input:
        "/home/data_admin/Reference/Annotation/gencode_47.full.gtf",
    output:
        "results/GRCh38_241031/expression/gene_id_mapping.txt",
    script:
        "extract_id_map.py"


rule expression_merge:
    input:
        quantity=lambda wildcards: expand(
            f"results/{wildcards.cohort}/expression/expression_sample/{{sample}}.gene.txt",
            sample=cohort.list_samples(),
            allow_missing=False,
        ),
        idmapping=rules.gene_id_mapping.output,
    output:
        rawcount="results/{cohort}/expression/expression_cohort/{cohort}.rawcount.txt",
        rawtpm="results/{cohort}/expression/expression_cohort/{cohort}.rawtpm.txt",
        log2tpm="results/{cohort}/expression/expression_cohort/{cohort}.log2tpm.txt",
    script:
        "merge_count.R"
