# -*- coding: UTF-8 -*-
#
# FileName     : SV_identify_sniffles2
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-04-24 11:25
# Last Modified: 2024-04-26 09:25
# Modified By  : EastsunW
# -------------
# Description  : 使用 sniffles 进行变异检测
# Citation     : https://github.com/fritzsedlazeck/Sniffles
# -------------


rule SV_sniffles2_call:
    input:
        bam=rules.alignment_merge.output.bam,
        reference=rules.alignment_movie.input.reference,
    output:
        vcf="results/{cohort}/SV/sample_sniffles2/{sample}.vcf.gz",
        vcf_index="results/{cohort}/SV/sample_sniffles2/{sample}.vcf.gz.tbi",
    log:
        "logs/{cohort}/SV/sample_sniffles2/{sample}.sniffles2_call.log",
    benchmark:
        "benchmarks/{cohort}/SV/sample_sniffles2/{sample}.sniffles2_call.benchmark"
    threads: config["SV"]["sniffles2"]["call_threads"]
    params:
        tendem_rep=f'--tandem-repeats {cohort.ref_config["SV"]["tandem_repeats"]}'
        if cohort.ref_config["SV"]["tandem_repeats"]
        else "",
        extra=config["SV"]["sniffles2"]["call_params"],
    conda:
        "../../envs/SV_sniffles2.yaml"
    shell:
        """
        (
            sniffles \
                --allow-overwrite \
                --input {input.bam} \
                --reference {input.reference} \
                --threads {threads} \
                --sample-id {wildcards.sample} \
                {params.tendem_rep} \
                {params.extra} \
                --vcf results/{wildcards.cohort}/SV/sample_sniffles2/{wildcards.sample}.vcf
            bgzip results/{wildcards.cohort}/SV/sample_sniffles2/{wildcards.sample}.vcf
            tabix -p vcf {output.vcf}
        ) &> {log}
        """


# rule SV_sniffles2_merge:
#     input:
#         expand(
#             "results/{{cohort}}/SV/sample_sniffles2/{sample}.snf",
#             sample=cohort.list_samples(),
#             allow_missing=False,
#         ),
#     output:
#         raw_vcf=temp("results/{cohort}/SV/cohort_sniffles2/{cohort}.raw.vcf.gz"),
#         sorted_vcf="results/{cohort}/SV/cohort_sniffles2/{cohort}.sorted.vcf.gz",
#         vcf_index="results/{cohort}/SV/cohort_sniffles2/{cohort}.sorted.vcf.gz.csi",
#     log:
#         "logs/{cohort}/SV/cohort_sniffles2/{cohort}.sniffles2_merge.log",
#     benchmark:
#         "benchmarks/{cohort}/SV/cohort_sniffles2/{cohort}.sniffles2_merge.benchmark"
#     params:
#         extra=config["SV"]["sniffles2"]["merge_params"],
#     threads: config["SV"]["sniffles2"]["merge_threads"]
#     conda:
#         "SV"
#     shell:
#         """
#         (
#             sniffles \
#                 --allow-overwrite \
#                 --threads {threads} \
#                 --input {input} \
#                 --vcf results/{wildcards.cohort}/SV/cohort_sniffles2/{wildcards.cohort}.raw.vcf
#             bgzip results/{wildcards.cohort}/SV/cohort_sniffles2/{wildcards.cohort}.raw.vcf
#             bcftools sort \
#                 --write-index \
#                 -Oz5 -o {output.sorted_vcf} \
#                 -T results/{wildcards.cohort}/SV/cohort_sniffles2 \
#                 {output.raw_vcf}
#         ) &> {log}
#         """
