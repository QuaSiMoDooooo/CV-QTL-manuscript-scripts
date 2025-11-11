# -*- coding: UTF-8 -*-
#
# FileName     : alignment_coverage
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-04-08 18:55
# Last Modified: 2024-04-08 18:55
# Modified By  : EastsunW
# -------------
# Description  :
# -------------


rule alignment_coverage:
    input:
        bam="results/{cohort}/alignment/aligned_bams/{sample}.merged.bam",
        index="results/{cohort}/alignment/aligned_bams/{sample}.merged.bam.csi",
        seq_len=cohort.ref_config["alignment"]["index"],
    output:
        bedgraph=temp("results/{cohort}/alignment/stat_coverage/{sample}.bedgraph"),
        summary="results/{cohort}/alignment/stat_coverage/{sample}.coverage.summary.txt",
    log:
        "logs/{cohort}/alignment/stat_coverage/{sample}.stat_coverage.log",
    benchmark:
        "benchmarks/{cohort}/alignment/stat_coverage/{sample}.stat_coverage.benchmark"
    conda:
        "../../envs/alignment.yaml"
    shell:
        """
        (genomeCoverageBed \
            -bg -split \
            -ibam {input.bam} \
            > \
            {output.bedgraph} && \
        python workflow/rules/alignment/coverage_calculate.py \
            --input {output.bedgraph} \
            --chr_len {input.seq_len} \
            --output {output.summary}) &> {log}
        """


rule collapse_coverage_stat:
    input:
        expand(
            "results/{{cohort}}/alignment/stat_coverage/{sample}.coverage.summary.txt",
            sample=cohort.list_samples(),
        ),
    output:
        "results/{cohort}/alignment/{cohort}.coverage.summary.txt",
    run:
        import pandas as pd
        import os
        import re

        result_dict = {}
        for filepath in input:
            sample_name = os.path.basename(filepath).split(".")[0]
            with open(filepath, "r") as df:
                df = pd.read_csv(filepath, sep="\t")
                df.set_index("chr", inplace=True)
                if sample_name not in result_dict:
                    result_dict[sample_name] = {}
                result_dict[sample_name]["coverage"] = df.loc["Total", "coverage"]
        result_df = pd.DataFrame.from_dict(result_dict, orient="index")
        result_df.index.name = "sample"
        result_df.to_csv(output[0], sep="\t")
