#!/usr/bin/env Rscript

# -------------
# FileName     : SNV_QTL.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Process cis-QTL data for SNP variants and filter variants with multiple molecular phenotypes
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

setwd("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/09_mediation_analysis/01_cis_QTL_overlap/SNP_all_QTL")

qtl_order = c("eQTL","sQTL","apaQTL","meQTL")
molphe_order = c("expression","splicing","APA","methylation")
qtl_dir = "/home/wtian/project/HZAU_cohort_meth/wdy_result/phenotype_qtl_results_bim3_removaNA"

file_list = list.files(path = qtl_dir, pattern = "^cis.filtered.txt.gz$", full.names = TRUE, recursive = TRUE)
file_list = file_list[grepl("SNP", file_list)]
file_list_molphe = str_replace_all(basename(dirname(dirname(file_list))), "SNP-", "")
file_list = file_list[order(match(file_list_molphe, molphe_order))]

data = fread(file_list[1], header = TRUE, sep = "\t")
data$molphe = qtl_order[1]
data$Gene = str_split(data$Phenotype,"\.", simplify = T)[,1]

for (i in 2:length(file_list)){
    qtl = qtl_order[i]
    data_temp = fread(file_list[i], header = TRUE, sep = "\t")
    data_temp$molphe = qtl
    
    if (qtl == "sQTL"){
        tmp = str_split(data_temp$Phenotype,":", simplify = T)
        data_temp$clu = paste0(tmp[,1],":",tmp[,4])
        
        tb = fread("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/08_trans_QTL_hotspot/05_transQTL_cisQTL_Genes/02_Genes_enrich/01_get_sGenes/02_clu_ensembl.tsv",sep = "\t", header = TRUE)
        colnames(tb) = c("clu","ENSEMBL")
        data_temp = left_join(data_temp, tb, by = "clu")
        data_temp = data_temp[!is.na(data_temp$ENSEMBL),]
        data_temp$Gene = data_temp$ENSEMBL
        data_temp$ENSEMBL = NULL
        data_temp$clu = NULL
    }
    if (qtl == "apaQTL"){
        data_temp$Gene = str_split(data_temp$Phenotype,"_", simplify = T)[,1]
    }
    
    data = rbind(data, data_temp)
}

data = data[order(data$Variant, data$molphe),]
data = select(data, Variant, molphe, Gene, everything())

fwrite(data, file = "all_cis_SNP_QTL.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

eQTL_var = unique(data$Variant[data$molphe == "eQTL"])
data_flt = data %>% filter(Variant %in% eQTL_var)

data_flt = data_flt %>%
  group_by(Variant) %>%
  filter(n_distinct(molphe) >= 2) %>%
  ungroup() %>%
  select(Variant, molphe, Gene, everything())

data_flt = as.data.frame(data_flt)
fwrite(data_flt, file = "all_cis_SNP_QTL_at_least_2molphe_including_eQTL.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
