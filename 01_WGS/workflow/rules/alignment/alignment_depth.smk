# -*- coding: UTF-8 -*-
#
# FileName     : alignment_stat
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-01-14 19:41
# Last Modified: 2024-01-14 19:41
# Modified By  : EastsunW
# -------------
# Description  : 统计比对结果的测序深度和覆盖度
# -------------


rule alignment_depth:
    input:
        "results/{cohort}/alignment/aligned_bams/{sample}.merged.bam",
    output:
        dist=temp(
            "results/{cohort}/alignment/stat_depth/{sample}.mosdepth.global.dist.txt"
        ),
        per_base="results/{cohort}/alignment/stat_depth/{sample}.per-base.bed.gz",
        per_base_idx="results/{cohort}/alignment/stat_depth/{sample}.per-base.bed.gz.csi",
        summary="results/{cohort}/alignment/stat_depth/{sample}.mosdepth.summary.txt",
    log:
        "logs/{cohort}/alignment/stat_depth/{sample}.stat_depth.log",
    benchmark:
        "benchmarks/{cohort}/alignment/stat_depth/{sample}.stat_depth.benchmark"
    shell:
        """
        (mosdepth \
            --threads 4 \
            results/{wildcards.cohort}/alignment/stat_depth/{wildcards.sample} \
            {input}) &> {log}
        """


rule high_depth_region:
    input:
        rules.alignment_depth.output.per_base,
    output:
        "results/{cohort}/alignment/stat_depth/{sample}.high_depth.bed.gz",
    params:
        depth=config["SV"]["filter"]["high_depth_threshold"],
    shell:
        """
        zcat {input} | awk '$4 > {params.depth}' | gzip -c > {output}
        """


rule collapse_depth_stat:
    input:
        expand(
            "results/{{cohort}}/alignment/stat_depth/{sample}.mosdepth.summary.txt",
            sample=cohort.list_samples(),
        ),
    output:
        "results/{cohort}/alignment/{cohort}.depth.summary.txt",
    run:
        import pandas as pd
        import os
        import re

        result_dict = {}
        for filepath in input:
            sample_name = os.path.basename(filepath).split(".")[0]
            with open(filepath, "r") as df:
                df = pd.read_csv(filepath, sep="\t")
                df.set_index("chrom", inplace=True)
                if sample_name not in result_dict:
                    result_dict[sample_name] = {}
                result_dict[sample_name]["depth"] = df.loc["total", "mean"]
        result_df = pd.DataFrame.from_dict(result_dict, orient="index")
        result_df.index.name = "sample"
        result_df.to_csv(output[0], sep="\t")
