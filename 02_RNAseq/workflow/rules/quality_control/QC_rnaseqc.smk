# -*- coding: UTF-8 -*-
#
# FileName     : QC_rnaseq2
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-04-15 16:23
# Last Modified: 2024-04-16 11:04
# Modified By  : EastsunW
# -------------
# Description  : 使用rnaseqc2进行质量控制信息生成，用于multiQC
# Citation     : 合并gtf：https://github.com/broadinstitute/gtex-pipeline/tree/master/gene_model
# Citation     : rnaseqc2：https://github.com/getzlab/rnaseqc
# -------------


rule rnaseqc2:
    input:
        "results/{cohort}/alignment/{sample}.bam",
    output:
        multiext(
            "results/{cohort}/QC/rnaseqc/{sample}/{sample}",
            ".metrics.tsv",
            ".coverage.tsv",
            ".exon_cv.tsv",
            ".exon_reads.gct",
            ".gene_tpm.gct",
            ".gene_reads.gct",
            ".gene_fragments.gct",
        ),
    log:
        "logs/{cohort}/QC/rnaseqc/{sample}.rnaseqc2.log",
    benchmark:
        "benchmarks/{cohort}/QC/rnaseqc/{sample}.rnaseqc2.benchmark"
    conda:
        "../../envs/QC.yaml"
    params:
        gtf=cohort.ref_config["quality_control"]["gtf"],
        extra=config["quality_control"]["extra_params"],
    shell:
        """
        (
            rnaseqc \
                {params.gtf} \
                {input} \
                results/{wildcards.cohort}/QC/rnaseqc/{wildcards.sample} \
                --sample {wildcards.sample} \
                {params.extra}
        ) &> {log}
        """
