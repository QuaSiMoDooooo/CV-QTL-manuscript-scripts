# -*- coding: UTF-8 -*-
#
# FileName     : SV_merge
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-06-21 16:47
# Last Modified: 2024-06-21 16:47
# Modified By  : EastsunW
# -------------
# Description  : 合并每个样本的三种软件的结果，分为两步，将每个样本的来自三个软件的结果（格式化后的）进行合并去冗余，得到每个样本，首先将三种软件的结果整合到一个文件中，再执行CAST聚类合并，合并后的SV会有新的ID，并且会在INFO中加入软件的信息
# -------------


rule SV_merge_sample:
    input:
        cutesv=rules.SV_format_cuteSV.output,
        svisionpro=rules.SV_format_SVisionpro.output,
        sniffles2=rules.SV_format_sniffles2.output,
        pbsv=rules.SV_format_pbsv.output,
    output:
        "results/{cohort}/SV/sample_merged/{sample}.raw.vcf.gz",
    log:
        "logs/{cohort}/SV/sample_merge/{sample}.merge.log",
    conda:
        "../../envs/SV_process.yaml"
    threads: 20
    params:
        extend_len=100,
        overlap_rate=0.5,
    script:
        "SV_merge_CAST.py"
