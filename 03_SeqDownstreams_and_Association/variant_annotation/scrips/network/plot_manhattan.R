#! Rscript
# -------------
# FileName     : process_GWAS
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2025-03-07 17:11
# Last Modified: 2025-03-07 21:13
# Modified By  : EastsunW
# -------------
# Description  : 处理GWAS结果，先画出曼哈顿图，然后筛选出显著的结果
# -------------

# 加载必要的包
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(doParallel))
setwd("/home/wangdy/Projects/Weibin/Downstreams/variant_annotation")


gwas_in_dir <- "data/GWAS"
