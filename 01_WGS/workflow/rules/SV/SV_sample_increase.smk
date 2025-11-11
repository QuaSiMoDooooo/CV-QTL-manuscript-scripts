# -*- coding: UTF-8 -*-
#
# FileName     : SV_sample_increase
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-11-28 22:27
# Last Modified: 2024-11-29 17:16
# Modified By  : EastsunW
# -------------
# Description  :
# -------------


rule SV_merge_cohort_step_1:
    input:
        f"results/{{cohort}}/SV/sample_filtered/{cohort.list_samples()[0]}.by_region.vcf.gz",
    output:
        "results/{cohort}/SV/cohort_step/data/{cohort}_1.vcf.gz",
    shell:
        "cp {input} {output}"

for n in range(2, len(cohort) + 1):
    rule:
        name: f"SV_merge_cohort_step_{n}",
        input:
            expand(
                "results/{{cohort}}/SV/sample_filtered/{sample}.by_region.vcf.gz",
                sample=cohort.list_samples()[0:n],
                allow_missing=False,
            ),
        output:
            f"results/{{cohort}}/SV/cohort_step/data/{{cohort}}_{n}.vcf.gz",
        params:
            extend_len=100,
            overlap_rate=0.5,
        script:
            "SV_merge_CAST.py"

rule SV_merge_cohort_stat:
    input:
        expand(
            "results/{{cohort}}/SV/cohort_step/data/{{cohort}}_{n}.vcf.gz",
            n=range(1, len(cohort) + 1),
        ),
    output:
        "results/{cohort}/SV/cohort_step/{cohort}_sample_increase_stat.txt",
    threads: 20
    script:
        "SV_sample_increase_stat.py"
