# -*- coding: UTF-8 -*-
#
# FileName     : SV_filter
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-06-21 16:48
# Last Modified: 2024-06-21 16:48
# Modified By  : EastsunW
# -------------
# Description  : 对合并后的SV进行过滤，首先保留至少两种软件的鉴定结果，接着有三种过滤条件，分别是基于深度的过滤，基于长度的过滤，基于区域的过滤
# -------------


rule SV_filter_software:
    input:
        rules.SV_merge_sample.output,
    output:
        "results/{cohort}/SV/sample_filtered/{sample}.by_software.vcf.gz",
    threads: 1
    params:
        mode="software",
        sample_name=lambda wildcards: wildcards.sample,
        software_order=config["SV"]["filter"]["software_order"],
    script:
        "SV_filter.py"


rule SV_filter_depth:
    input:
        rules.SV_filter_software.output,
    output:
        "results/{cohort}/SV/sample_filtered/{sample}.by_depth.vcf.gz",
    threads: 1
    params:
        mode="depth",
        min_depth=config["SV"]["filter"]["min_reads"],
    script:
        "SV_filter.py"


rule SV_filter_length:
    input:
        rules.SV_filter_depth.output,
    output:
        "results/{cohort}/SV/sample_filtered/{sample}.by_length.vcf.gz",
    threads: 1
    params:
        mode="length",
        length_field=config["SV"]["filter"]["length_field"],
        max_len_ins=config["SV"]["filter"]["max_len_ins"],
        max_len_del=config["SV"]["filter"]["max_len_del"],
        max_len_inv=config["SV"]["filter"]["max_len_inv"],
        max_len_dup=config["SV"]["filter"]["max_len_dup"],
    script:
        "SV_filter.py"


rule SV_filter_region:
    input:
        sv=rules.SV_filter_length.output,
        high_depth_region=rules.high_depth_region.output,
    output:
        "results/{cohort}/SV/sample_filtered/{sample}.by_region.vcf.gz",
    threads: 1
    params:
        mode="region",
        remove_contigs=config["SV"]["filter"]["remove_contig"],
        centromeres_region=cohort.ref_config["SV"]["centromeres_region"],
        gap_telomere_region=cohort.ref_config["SV"]["gap_telomere_region"],
    script:
        "SV_filter.py"
