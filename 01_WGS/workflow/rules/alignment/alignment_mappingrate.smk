# -*- coding: UTF-8 -*-
#
# FileName     : alignment_mappingrate
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-11-02 09:04
# Last Modified: 2024-11-05 09:17
# Modified By  : EastsunW
# -------------
# Description  : 计算比对率（samtools）
# -------------


rule alignment_mappingrate:
    input:
        bam="results/{cohort}/alignment/aligned_bams/{sample}.merged.bam",
        index="results/{cohort}/alignment/aligned_bams/{sample}.merged.bam.csi",
    output:
        "results/{cohort}/alignment/stat_mappingrate/{sample}.mappingrate.summary.txt",
    threads: 10
    shell:
        """
        samtools flagstat -@ {threads} {input.bam} > {output}
        """


rule collapse_mappingrate_stat:
    input:
        expand(
            "results/{{cohort}}/alignment/stat_mappingrate/{sample}.mappingrate.summary.txt",
            sample=cohort.list_samples(),
        ),
    output:
        "results/{cohort}/alignment/{cohort}.mappingrate.summary.txt",
    run:
        import pandas as pd
        import os
        import re

        result_dict = {}
        for filepath in input:
            sample_name = os.path.basename(filepath).split(".")[0]
            with open(filepath, "r") as df:
                match = re.search(r"\((\d+\.\d+)%", df.readlines()[6])
                if match:
                    mapping_rate = float(match.group(1)) / 100
                else:
                    mapping_rate = 0
                if sample_name not in result_dict:
                    result_dict[sample_name] = {}
                result_dict[sample_name]["mapping_rate"] = mapping_rate
        result_df = pd.DataFrame.from_dict(result_dict, orient="index")
        result_df.index.name = "sample"
        result_df.to_csv(output[0], sep="\t")
