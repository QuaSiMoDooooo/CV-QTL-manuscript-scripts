# -*- coding: UTF-8 -*-
#
# FileName     : QC_multiQC
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-04-15 16:25
# Last Modified: 2024-04-15 16:25
# Modified By  : EastsunW
# -------------
# Description  :
# -------------


rule multiqc:
    input:
        hisat=expand(
            "results/{{cohort}}/alignment/{sample}.report",
            sample=cohort.list_samples(),
        ),
        salmon=expand(
            "results/{{cohort}}/expression/expression_sample/{sample}.lib_format_counts.json",
            sample=cohort.list_samples(),
        ),
        rnaseqc=expand(
            "results/{{cohort}}/QC/rnaseqc/{sample}/{sample}.metrics.tsv",
            sample=cohort.list_samples(),
        ),
    output:
        "results/{cohort}/report/{cohort}.multiqc.html",
    shell:
        """
            multiqc \
                -o results/{wildcards.cohort}/report \
                -f -q \
                -n {wildcards.cohort}.multiqc.html \
                results/{wildcards.cohort}
        """
