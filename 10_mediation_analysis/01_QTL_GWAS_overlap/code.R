#!/usr/bin/env Rscript

# -------------
# FileName     : code.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Identify overlapping variants between QTL and GWAS data
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

setwd("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/09_mediation_analysis/01_QTL_GWAS_overlap")

qtl_dir = "/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/09_mediation_analysis/01_cis_QTL_overlap"
qtl_list = c("eQTL","sQTL","apaQTL","meQTL")
qtl_name_list = c("expression","splicing","APA","methylation")
variant_type_list = c("SNP","SV")

gwas_dir = "/home/wangdy/Projects/Weibin/Downstreams/data_stat/results/stat/GWAS/"

for (i in 1:length(variant_type_list)){
    var = variant_type_list[i]
    if (!dir.exists(var)){
        dir.create(var)
    }
    
    qtl_path = paste0(qtl_dir,"/",var,"_all_QTL/all_cis_",var,"_QTL.tsv")
    qtl = fread(qtl_path)
    
    gwas_subdir = paste0(gwas_dir,"/",var,"/")
    gwas_paths = list.files(gwas_subdir, pattern = "*filtered.tsv", full.names = TRUE)
    gwas_paths = gwas_paths[!grepl("PHENO1", gwas_paths)]
    
    gwas = data.frame(matrix(ncol = 13, nrow = 0))
    for (j in 1:length(gwas_paths)){
        gwas_tmp = fread(gwas_paths[j])
        gwas = rbind(gwas, gwas_tmp)
    }
    
    gwas = select(gwas, Indicator, everything())
    colnames(gwas)[3] = "Variant"
    
    gwas_variants = unique(gwas$Variant)
    qtl_flt = qtl %>% filter(Variant %in% gwas_variants)
    
    fwrite(qtl_flt, file = paste0(var,"/all_cis_QTL_overlap.tsv"), sep = "\t", quote = FALSE)
    
    qtl_flt_variants = unique(qtl_flt$Variant)
    gwas_flt = gwas %>% filter(Variant %in% qtl_flt_variants)
    
    fwrite(gwas_flt, file = paste0(var,"/all_GWAS_overlap.tsv"), sep = "\t", quote = FALSE)
}