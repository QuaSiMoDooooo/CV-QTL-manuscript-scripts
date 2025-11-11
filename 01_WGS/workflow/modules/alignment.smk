# -*- coding: UTF-8 -*-
#
# FileName     : alignment
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-04-10 15:11
# Last Modified: 2024-04-10 15:11
# Modified By  : EastsunW
# -------------
# Description  :
# -------------

conda: "../envs/alignment.yaml"


include: "../rules/alignment/alignment_pbmm2.smk"
include: "../rules/alignment/alignment_depth.smk"
include: "../rules/alignment/alignment_coverage.smk"
include: "../rules/alignment/alignment_mappingrate.smk"
