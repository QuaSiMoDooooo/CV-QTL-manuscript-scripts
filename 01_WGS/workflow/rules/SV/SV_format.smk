# -*- coding: UTF-8 -*-
#
# FileName     : SV_format
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-06-19 21:53
# Last Modified: 2024-06-19 21:53
# Modified By  : EastsunW
# -------------
# Description  :
# -------------


rule SV_format_cuteSV:
    input:
        rules.SV_cuteSV_call.output.vcf,
    output:
        "results/{cohort}/SV/sample_formated/cuteSV/{sample}.formated.vcf.gz",
    threads: 1
    params:
        software="cutesv",
        af_threshold=0.2,
    script:
        "SV_format.py"


rule SV_format_pbsv:
    input:
        rules.SV_pbsv_call.output.vcf,
    output:
        "results/{cohort}/SV/sample_formated/pbsv/{sample}.formated.vcf.gz",
    threads: 1
    params:
        software="pbsv",
        af_threshold=0.2,
    script:
        "SV_format.py"


rule SV_format_sniffles2:
    input:
        rules.SV_sniffles2_call.output.vcf,
    output:
        "results/{cohort}/SV/sample_formated/sniffles2/{sample}.formated.vcf.gz",
    threads: 1
    params:
        software="sniffles2",
        af_threshold=0.2,
    script:
        "SV_format.py"


rule SV_format_SVisionpro:
    input:
        rules.gzip_SVisionpro.output.simpleSV,
    output:
        "results/{cohort}/SV/sample_formated/SVision/{sample}.formated.vcf.gz",
    threads: 1
    params:
        software="svisionpro",
        af_threshold=0.2,
    script:
        "SV_format.py"
