# -*- coding: UTF-8 -*-
#
# FileName     : SV_merge_cohort
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-07-02 22:23
# Last Modified: 2024-07-02 22:23
# Modified By  : EastsunW
# -------------
# Description  :
# -------------


rule SV_merge_cohort:
    input:
        expand(
            "results/{{cohort}}/SV/sample_filtered/{sample}.by_region.vcf.gz",
            sample=cohort.list_samples(),
            allow_missing=False,
        ),
    output:
        temp("results/{cohort}/SV/cohort_merged/{cohort}.raw.vcf.gz"),
    conda:
        "../../envs/SV_process.yaml"
    threads: 20
    params:
        extend_len=100,
        overlap_rate=0.5,
    script:
        "SV_merge_CAST.py"


rule SV_merge_cohort_addAF:
    input:
        rules.SV_merge_cohort.output,
    output:
        "results/{cohort}/SV/cohort_merged/{cohort}.SV.vcf.gz",
    threads: 1
    script:
        "SV_add_AF_cohort.py"
