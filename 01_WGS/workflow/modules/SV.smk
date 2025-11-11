# -*- coding: UTF-8 -*-
#
# FileName     : SV
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-04-10 15:18
# Last Modified: 2024-04-10 15:18
# Modified By  : EastsunW
# -------------
# Description  :
# -------------


include: "../rules/SV/SV_identify_pbsv.smk"
include: "../rules/SV/SV_identify_cuteSV.smk"
include: "../rules/SV/SV_identify_sniffles2.smk"
include: "../rules/SV/SV_identify_SVisionpro.smk"
# 格式化所有软件的结果
include: "../rules/SV/SV_format.smk"
# 样本SV合并、筛选
include: "../rules/SV/SV_merge_sample.smk"
include: "../rules/SV/SV_filter.smk"
# 队列SV的合并
include: "../rules/SV/SV_merge_cohort.smk"
# SV鉴定随样本增加而变化
include: "../rules/SV/SV_sample_increase.smk"
