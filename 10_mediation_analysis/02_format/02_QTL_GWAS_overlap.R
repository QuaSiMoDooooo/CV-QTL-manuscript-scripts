#!/usr/bin/env Rscript

# -------------
# FileName     : 02_QTL_GWAS_overlap.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Format QTL and GWAS overlap data for mediation analysis
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))

mc = getOption("mc.cores", 100)
setwd("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/09_mediation_analysis/02_format")

qtl_list = c("eQTL","sQTL","apaQTL","meQTL")
qtl_name_list = c("expression","splicing","APA","methylation")
variant_type_list = c("SNP","SV")

sample_info = fread("sample_info.txt", data.table = FALSE)
sample_info$gender = ifelse(sample_info$gender == "男", 1, 0)
sample_info = sample_info[, c(1, 2, 3)]
colnames(sample_info) = c("sample_id", "gender", "age")

expression = fread("expression.quantity.qnorm.txt", data.table = FALSE)
splicing = fread("splicing.quantity.qnorm.txt", data.table = FALSE)
APA = fread("APA.quantity.qnorm.txt", data.table = FALSE)
methylation = fread("methylation.quantity.qnorm.txt", data.table = FALSE)

input_dir = "/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/09_mediation_analysis/01_QTL_GWAS_overlap"
out_dir = "02_QTL_GWAS_overlap/"

gwas = fread("biochem_common.norm.format.txt", data.table = FALSE)

for (j in 1:length(variant_type_list)) {
    variant_type = variant_type_list[j]
    variant_path = paste0(variant_type, ".biochem_common.raw.format.t.ID.txt")
    variant = fread(variant_path, data.table = FALSE)
    var_dir = paste0(out_dir, variant_type)
    if (!dir.exists(var_dir)) {
        dir.create(var_dir, recursive = TRUE)
    }

    gwas_asso_path = paste0(input_dir, "/", variant_type, "/all_GWAS_overlap.tsv")
    gwas_asso = fread(gwas_asso_path, data.table = FALSE)
    gwas_asso = gwas_asso[gwas_asso$IndicatorType == "common", ]
    ind_list = unique(gwas_asso$Indicator)

    asso_path = paste0(input_dir, "/", variant_type, "/all_cis_QTL_overlap.tsv")
    asso = fread(asso_path, data.table = FALSE)

    for (ind in ind_list) {
        e_phe = ind
        e_df = data.frame(sample_id = colnames(gwas)[2:length(colnames(gwas))], 
                         e_values = unname(unlist(gwas[gwas$ID == e_phe, 2:length(colnames(gwas))])))
        var_list = unique(gwas_asso[gwas_asso$Indicator == ind, ]$Variant)
        
        mclapply(var_list, function(var) {
            var_df = data.frame(sample_id = colnames(variant)[2:length(colnames(variant))], 
                               var_values = unname(unlist(variant[variant$ID == var, 2:length(colnames(variant))])))
            asso_var = asso[asso$Variant == var, ]
            
            for (m in 1:length(asso_var$Phenotype)) {
                m_phe = asso_var$Phenotype[m]
                m_molphe = asso_var$molphe[m]
                out_file = paste0(out_dir, variant_type, "/", var, ":", e_phe, ":", m_phe, ":mediation.input")
                
                if (!file.exists(out_file)) {
                    if (m_molphe == "eQTL") { 
                        m_df = data.frame(sample_id = colnames(expression)[2:length(colnames(expression))], 
                                         m_values = unname(unlist(expression[expression$ID == m_phe, 2:length(colnames(expression))])))
                    }
                    if (m_molphe == "sQTL") { 
                        m_df = data.frame(sample_id = colnames(splicing)[2:length(colnames(splicing))], 
                                         m_values = unname(unlist(splicing[splicing$ID == m_phe, 2:length(colnames(splicing))])))
                    }
                    if (m_molphe == "apaQTL") { 
                        m_df = data.frame(sample_id = colnames(APA)[2:length(colnames(APA))], 
                                         m_values = unname(unlist(APA[APA$ID == m_phe, 2:length(colnames(APA))])))
                    }
                    if (m_molphe == "meQTL") { 
                        m_df = data.frame(sample_id = colnames(methylation)[2:length(colnames(methylation))], 
                                         m_values = unname(unlist(methylation[methylation$ID == m_phe, 2:length(colnames(methylation))])))
                    }
                    
                    df = left_join(var_df, sample_info, by = "sample_id") %>%
                         left_join(e_df, by = "sample_id") %>%
                         left_join(m_df, by = "sample_id")
                    
                    fwrite(df, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
                }
            }
        }, mc.cores = mc)
    }
}