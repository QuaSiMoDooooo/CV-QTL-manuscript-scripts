# -*- coding: UTF-8 -*-
#
# FileName     : pbsv
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-01-19 19:08
# Last Modified: 2024-01-19 19:08
# Modified By  : EastsunW
# -------------
# Description  : 使用 pbsv 进行变异检测
# -------------


rule SV_pbsv_discover:
    input:
        rules.alignment_merge.output.bam,
    output:
        temp("results/{cohort}/SV/sample_pbsv/sv_discover/{sample}.svsig.gz"),
    log:
        "logs/{cohort}/pbsv/{sample}.sv_discover.log",
    params:
        tendem_rep=f'--tandem-repeats {cohort.ref_config["SV"]["tandem_repeats"]}'
        if cohort.ref_config["SV"]["tandem_repeats"]
        else "",
        extra=config["SV"]["pbsv"]["discover"]["extra_params"],
    conda:
        "../../envs/SV_pbsv.yaml"
    shell:
        """
        (
            pbsv discover \
                --hifi \
                --sample {wildcards.sample} \
                {params.tendem_rep} \
                {params.extra} \
                {input} \
                {output}
        ) &> {log}
        """


rule SV_pbsv_call:
    input:
        reference=rules.alignment_movie.input.reference,
        svsig=rules.SV_pbsv_discover.output,
    output:
        vcf="results/{cohort}/SV/sample_pbsv/sv_call/{sample}.vcf.gz",
        index="results/{cohort}/SV/sample_pbsv/sv_call/{sample}.vcf.gz.tbi",
    threads: config["SV"]["pbsv"]["call"]["threads"]
    log:
        "logs/{cohort}/pbsv/{sample}.sv_call.log",
    benchmark:
        "benchmarks/{cohort}/pbsv/{sample}.sv_call.benchmark"
    params:
        extra=config["SV"]["pbsv"]["call"]["extra_params"],
    conda:
        "../../envs/SV_pbsv.yaml"
    shell:
        """
        (
            pbsv call \
                --num-threads {threads} \
                --hifi \
                {params.extra} \
                {input.reference} \
                {input.svsig} \
                results/{wildcards.cohort}/SV/sample_pbsv/sv_call/{wildcards.sample}.vcf
            bgzip results/{wildcards.cohort}/SV/sample_pbsv/sv_call/{wildcards.sample}.vcf
            tabix -p vcf {output.vcf}
        ) &> {log}
        """
