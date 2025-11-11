# -*- coding: UTF-8 -*-
#
# FileName     : alignment_universal
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-01-29 18:29
# Last Modified: 2024-01-29 18:29
# Modified By  : EastsunW
# -------------
# Description  : 使用hisat2将fastq文件比对到参考基因组
# -------------


rule alignment:
    input:
        seq_1=lambda wildcards: f"{cohort.input_dir}/{cohort.list_files(wildcards.sample)[0]}",
        seq_2=lambda wildcards: f"{cohort.input_dir}/{cohort.list_files(wildcards.sample)[1]}",
    output:
        raw=temp("results/{cohort}/alignment/{sample}.raw.bam"),
        sorted="results/{cohort}/alignment/{sample}.bam",
        index="results/{cohort}/alignment/{sample}.bam.bai",
        report="results/{cohort}/alignment/{sample}.report",
    threads: config["alignment"]["threads"]
    retries: 3
    log:
        "logs/{cohort}/alignment/{sample}.alignment.log",
    benchmark:
        "benchmarks/{cohort}/alignment/{sample}.alignment.benchmark"
    params:
        hisat_idx=cohort.ref_config["alignment"]["hisat_index"],
        extra=config["alignment"]["extra_params"],
    shell:
        """
        (
            hisat2 \
                -p {threads} \
                -x {params.hisat_idx} \
                -1 {input.seq_1} \
                -2 {input.seq_2} \
                --summary-file {output.report} \
                {params.extra} | \
                samtools view \
                    -@ {threads} \
                    -b -S \
                    -o {output.raw}

            samtools sort \
                {output.raw} \
                -@ {threads} \
                -o {output.sorted}

            samtools index \
                -@ {threads} \
                -o {output.index} \
                {output.sorted}
        ) &> {log}
        """
