# -*- coding: UTF-8 -*-
#
# FileName     : alternative_splice
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-02-27 20:24
# Last Modified: 2024-02-27 20:24
# Modified By  : EastsunW
# -------------
# Description  : 使用 leafcutter 鉴定 AS
# -------------


rule AS_bam2junc:
    input:
        rules.alignment.output.sorted,
    output:
        temp("results/{cohort}/AS/junc_files/{sample}.junc"),
    log:
        "logs/{cohort}/AS/bam2junc/{sample}.bam2junc.log",
    benchmark:
        "benchmarks/{cohort}/AS/bam2junc/{sample}.bam2junc.benchmark"
    container:
        "docker://dockerproxy.net/eastsunw/leafcutter:latest"
    params:
        extra=config["AS"]["bam2junc"]["extra_params"],
    shell:
        """
            regtools junctions extract {params.extra} {input} -o {output} &> {log}
        """


rule AS_intron_clustering:
    input:
        expand(
            "results/{{cohort}}/AS/junc_files/{sample}.junc",
            sample=cohort.list_samples(),
        ),
    output:
        "results/{cohort}/AS/intron_clustering/as_perind.counts.gz",
    log:
        "logs/{cohort}/AS/{cohort}.intron_clustering.log",
    benchmark:
        "benchmarks/{cohort}/AS/{cohort}.intron_clustering.benchmark"
    container:
        "docker://dockerproxy.net/eastsunw/leafcutter:latest"
    params:
        extra=config["AS"]["intron_clustering"]["extra_params"],
    shell:
        """
        (
            temp=({input})
            printf "%s\\n" "${{temp[@]}}" > results/{wildcards.cohort}/AS/junc_file_list.txt
            python /opt/leafcutter/clustering/leafcutter_cluster_regtools.py \
                {params.extra} \
                -j results/{wildcards.cohort}/AS/junc_file_list.txt \
                -o as \
                -r results/{wildcards.cohort}/AS/intron_clustering
        ) &> {log}
        """


rule AS_post_process:
    input:
        rules.AS_intron_clustering.output,
    output:
        "results/{cohort}/AS/intron_clustering/as_perind.counts.gz.PCs",
    log:
        "logs/{cohort}/AS/{cohort}.post_process.log",
    benchmark:
        "benchmarks/{cohort}/AS/{cohort}.post_process.benchmark"
    container:
        "docker://dockerproxy.net/eastsunw/leafcutter:latest"
    params:
        exon_file=cohort.ref_config["AS"]["exon"],
    shell:
        """
        (
            /opt/leafcutter/scripts/leafcutter_ds.R \
                -e {params.exon_file} \
                -o 
                {input}
        ) &> {log}
        """


rule AS_merge:
    input:
        rules.AS_post_process.output,
    output:
        "results/{cohort}/AS/{cohort}.AS.txt",
    shell:
        """
            head -n 1 results/{wildcards.cohort}/AS/intron_clustering/as_perind.counts.gz.phen_chr1 > {output}
            grep -v -h "^#" results/{wildcards.cohort}/AS/intron_clustering/as_perind.counts.gz.phen_chr* >> {output}
        """
