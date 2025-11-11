# -*- coding: UTF-8 -*-
#
# FileName     : alignment_universal
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-01-29 18:29
# Last Modified: 2024-01-29 18:29
# Modified By  : EastsunW
# -------------
# Description  :
# -------------


rule alignment_movie:
    input:
        seq_file=lambda wildcards: f"{cohort.input_dir}/{wildcards.sample}__{wildcards.movie}.{cohort.file_suffix}",
        reference=cohort.ref_config["alignment"]["fasta"],
    output:
        mapped=temp(
            "results/{cohort}/alignment/aligned_bams/{sample}__{movie}.aligned.bam"
        ),
        index=temp(
            "results/{cohort}/alignment/aligned_bams/{sample}__{movie}.aligned.bam.csi"
        ),
    threads: config["alignment"]["threads"]
    log:
        "logs/{cohort}/alignment/align_movie/{sample}__{movie}.align_movie.log",
    benchmark:
        "benchmarks/{cohort}/alignment/align_movie/{sample}__{movie}.align_movie.benchmark"
    params:
        sample=lambda wildcards: wildcards.sample,
        extra=config["alignment"]["extra_params"],
    shell:
        """
        pbmm2 align \
            --num-threads {threads} \
            --sample {params.sample} \
            --log-file {log} \
            --bam-index CSI \
            {params.extra} \
            {input.reference} \
            {input.seq_file} \
            {output.mapped}
        """


rule alignment_merge:
    input:
        lambda wildcards: expand(
            "results/{{cohort}}/alignment/aligned_bams/{{sample}}__{movie}.aligned.bam",
            movie=cohort.list_movies(wildcards.sample),
        ),
    output:
        bam="results/{cohort}/alignment/aligned_bams/{sample}.merged.bam",
        index="results/{cohort}/alignment/aligned_bams/{sample}.merged.bam.csi",
    threads: lambda wildcards: len(cohort.list_files(sample=wildcards.sample)) * 2
    shell:
        """
        samtools merge \
            --threads {threads} \
            --write-index \
            -o {output.bam} \
            {input}
        samtools index -c {output.bam}
        """
