# -*- coding: UTF-8 -*-
#
# FileName     : QC_rseqc
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-05-20 10:31
# Last Modified: 2024-05-20 10:31
# Modified By  : EastsunW
# -------------
# Description  : 使用RSeQC进行RNA-seq数据的质控
# Citation     : https://rseqc.sourceforge.net/
# -------------


rule geneBodyCoverage_filelist:
    input:
        expand(
            "results/{{cohort}}/alignment/{sample}.bam",
            sample=cohort.list_samples(),
        ),
    output:
        "results/{cohort}/QC/rseqc/{cohort}.BAM_files_list.txt",
    shell:
        """
        temp=({input})
        printf "%s\\n" "${{temp[@]}}" > {output}
        """


rule geneBodyCoverage:
    input:
        ancient(rules.geneBodyCoverage_filelist.output),
    output:
        "results/{cohort}/QC/rseqc/{cohort}.geneBodyCoverage.txt",
    log:
        "logs/{cohort}/QC/rseqc/{cohort}.rseqc.log",
    benchmark:
        "benchmarks/{cohort}/QC/rseqc/{cohort}.rseqc.benchmark"
    conda:
        "../../envs/QC.yaml"
    shell:
        """
        (geneBody_coverage.py \
            -r resources/reference/QC/hg38.HouseKeepingGenes.bed \
            -i {input}  \
            -o results/{wildcards.cohort}/QC/rseqc/{wildcards.cohort}
        ) &> {log}
        """
