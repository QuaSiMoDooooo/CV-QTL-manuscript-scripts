#!/usr/bin/env Rscript

# -------------
# FileName     : 03_smr_script.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Perform SMR analysis for QTL-GWAS associations
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(vroom))
suppressPackageStartupMessages(library(stringr))
setwd("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/07_QTL_GWAS_SMR")

qtl_file = commandArgs(trailingOnly = TRUE)[1]

qtl = basename(qtl_file)
type = str_split(qtl,"_",simplify = T)[,1]
output_dir = paste0("03_smr/", qtl)
if (!dir.exists(output_dir)){
    dir.create(output_dir, recursive = TRUE)
}

gwas_dir = "GWAS/"
gwas_file_paths = list.files(path = gwas_dir, pattern = "\\.ma$", recursive = TRUE, full.names = TRUE)

smr = function(qtl_file, type, output_dir, gwas_path){
    trait = str_split(gwas_path,"\\.",simplify = T)[,4]
    output_prefix = paste0(output_dir, "/", trait)
    output_path = paste0(output_dir, "/", trait, ".smr")
    if (!file.exists(output_path)){
        if(type == "cis"){
            cmd_cis = paste0("smr --bfile ~/data/LD_panel/mrcieu/EAS --gwas-summary ", gwas_path, " --beqtl-summary ", qtl_file, " --out ", output_prefix, " --thread-num 8")
            system(cmd_cis)
        }else if(type == "trans"){
            cmd_trans = paste0("smr --bfile ~/data/LD_panel/mrcieu/EAS --gwas-summary ", gwas_path, " --beqtl-summary ", qtl_file, " --out ", output_prefix, " --trans --trans-wind 1000 --thread-num 8")
            system(cmd_trans)
        }
    }
}

library(parallel)
mc = getOption("mc.cores", 16)
mclapply(1:length(gwas_file_paths), function(i){ smr(qtl_file = qtl_file, type = type, output_dir= output_dir, gwas_file_paths[i]) }, mc.cores = mc)