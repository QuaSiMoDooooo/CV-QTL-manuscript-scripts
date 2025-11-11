# -*- coding: UTF-8 -*-
#
# FileName     : SV_identify_cuteSV
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-04-24 11:25
# Last Modified: 2024-04-26 21:16
# Modified By  : EastsunW
# -------------
# Description  : 使用 cuteSV 进行变异检测
# Citation     : https://github.com/tjiangHIT/cuteSV
# -------------


rule SV_cuteSV_call:
    input:
        bam=rules.alignment_merge.output.bam,
        reference=rules.alignment_movie.input.reference,
    output:
        vcf="results/{cohort}/SV/sample_cuteSV/{sample}.vcf.gz",
        vcf_index="results/{cohort}/SV/sample_cuteSV/{sample}.vcf.gz.tbi",
        temp_dir=temp(directory("results/{cohort}/SV/sample_cuteSV/{sample}")),
    log:
        "logs/{cohort}/SV/sample_cuteSV/{sample}.cuteSV_call.log",
    benchmark:
        "benchmarks/{cohort}/SV/sample_cuteSV/{sample}.cuteSV_call.benchmark"
    threads: config["SV"]["cuteSV"]["threads"]
    params:
        extra=config["SV"]["cuteSV"]["extra_params"],
    conda:
        "../../envs/SV_cutesv.yaml"
    shell:
        """
        (
            if [ ! -d {output.temp_dir} ]; then
                mkdir -p {output.temp_dir}
            fi
            cuteSV \
                --threads {threads} \
                --sample {wildcards.sample} \
                {params.extra} \
                {input.bam} \
                {input.reference} \
                results/{wildcards.cohort}/SV/sample_cuteSV/{wildcards.sample}.vcf \
                {output.temp_dir}
            bgzip results/{wildcards.cohort}/SV/sample_cuteSV/{wildcards.sample}.vcf
            tabix -p vcf {output.vcf}
        ) &> {log}
        """
