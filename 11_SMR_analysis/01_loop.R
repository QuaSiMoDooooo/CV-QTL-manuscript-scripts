#!/usr/bin/env Rscript

# -------------
# FileName     : 01_loop.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Run scripts in parallel
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))

qtl_dir = "/home/wangdy/Projects/Weibin/Downstreams/phenotype_qtl/results"
qtl_list = c("eQTL","sQTL","apaQTL","meQTL")
qtl_name_list = c("expression","splicing","APA","methylation")
type_list = c("cis","trans")

for (i in 1:length(qtl_list)){
    qtl = qtl_list[i]
    qtl_name = qtl_name_list[i]
    for (j in 1:length(type_list)){
        type = type_list[j]
        command = paste("nohup Rscript 01_script.R", qtl, qtl_name, type, "&")
        print(command)
        system(command)
    }
}