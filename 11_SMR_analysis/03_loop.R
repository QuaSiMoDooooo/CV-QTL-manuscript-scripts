#!/usr/bin/env Rscript

# -------------
# FileName     : 03_loop.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Run SMR analysis scripts
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(vroom))
suppressPackageStartupMessages(library(stringr))
setwd("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/07_QTL_GWAS_SMR")

input_dir = "/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/07_QTL_GWAS_SMR/02_besd"
dirs <- list.dirs(input_dir, full.names = TRUE, recursive = FALSE)

cis_dir = dirs[grepl("cis", dirs)]

for (dir in cis_dir) {
    dir_name <- basename(dir)
    file <- paste0(dir, "/", dir_name)
    cmd = paste0("nohup Rscript 03_smr_script.R ", file, " &")
    system(cmd)
}