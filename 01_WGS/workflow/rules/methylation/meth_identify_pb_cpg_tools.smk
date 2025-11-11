# -*- coding: UTF-8 -*-
#
# FileName     : 5mc
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-02-01 20:11
# Last Modified: 2024-02-01 20:11
# Modified By  : EastsunW
# -------------
# Description  : 使用 pb_cpg_tools 进行 5mc 鉴定
# -------------


rule methylation_identify:
    input:
        bam=rules.alignment_merge.output.bam,
        bam_index=rules.alignment_merge.output.index,
    output:
        cpg_bed="results/{cohort}/methylation/methylation_sample/{sample}.combined.bed",
        cpg_bw="results/{cohort}/methylation/methylation_sample/{sample}.combined.bw",
        temp_log=temp("results/{cohort}/methylation/methylation_sample/{sample}.log"),
    log:
        "logs/{cohort}/methylation/methylation_sample/{sample}.cpg_identify.log",
    benchmark:
        "benchmarks/{cohort}/methylation/methylation_sample/{sample}.cpg_identify.benchmark"
    params:
        extra=config["methylation"]["pb_cpg_tool"]["extra_params"],
    threads: config["methylation"]["pb_cpg_tool"]["threads"]
    container:
        f"docker://quay.dockerproxy.net/pacbio/pb-cpg-tools:{config['methylation']['pb_cpg_tool']['version']}"
    shell:
        """
        (
            aligned_bam_to_cpg_scores \
                --threads {threads} \
                --bam {input.bam} \
                --output-prefix results/{wildcards.cohort}/methylation/methylation_sample/{wildcards.sample} \
                {params.extra}
        ) &> {log}
        """

rule methylation_merge:
    input:
        lambda wildcards: expand(
            "results/{{cohort}}/methylation/methylation_sample/{sample}.combined.bed",
            sample=cohort.list_samples(),
            allow_missing=False,
        ),
    output:
        "results/{cohort}/methylation/methylation_cohort/{cohort}.merged.txt",
    script:
        f"methylation_merge.R"
