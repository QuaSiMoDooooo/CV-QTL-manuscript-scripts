#!/usr/bin/env Rscript

# -------------
# FileName     : 01_cis_QTL_overlap.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Format cis-QTL overlap data for mediation analysis
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))

mc = getOption("mc.cores", 100)
setwd("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/09_mediation_analysis/02_format")

qtl_list = c("eQTL","sQTL","apaQTL","meQTL")
qtl_name_list = c("expression","splicing","APA","methylation")
variant_type_list = c("SNP","SV")

sample_info = fread("sample_info.txt",data.table = F)
sample_info$gender = ifelse(sample_info$gender == "男",1,0)
sample_info = sample_info[,c(1,2,3)]
colnames(sample_info) = c("sample_id","gender","age")

expression = fread("expression.quantity.qnorm.txt",data.table = F)
splicing = fread("splicing.quantity.qnorm.txt",data.table = F)
APA = fread("APA.quantity.qnorm.txt",data.table = F)
methylation = fread("methylation.quantity.qnorm.txt",data.table = F)

input_dir = "/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/09_mediation_analysis/01_cis_QTL_overlap"
out_dir = "01_cis_QTL_overlap/"

for (j in 1:length(variant_type_list)){
    variant_type = variant_type_list[j]
    variant_path = paste0(variant_type,".biochem_common.raw.format.t.ID.txt")
    variant = fread(variant_path,data.table = F)
    var_dir = paste0(out_dir,variant_type)
    if (!dir.exists(var_dir)){
        dir.create(var_dir, recursive = T)
    }
    asso_path=paste0(input_dir,"/",variant_type,"_all_QTL/all_cis_",variant_type,"_QTL_at_least_2molphe_including_eQTL.tsv")
    asso = fread(asso_path,data.table = F)
    var_list = unique(asso$Variant)
    
    mclapply(var_list, function(var){
        var_df = data.frame(sample_id=colnames(variant)[2:length(colnames(variant))],var_values=unname(unlist(variant[variant$ID == var,2:length(colnames(variant))])))
        asso_var = asso[asso$Variant == var,]
        asso_eQTL = asso_var[asso_var$molphe == "eQTL",]
        asso_others = asso_var[asso_var$molphe != "eQTL",]
        for (e in 1:length(asso_eQTL$Phenotype)){
            e_phe = asso_eQTL$Phenotype[e]
            e_df = data.frame(sample_id=colnames(expression)[2:length(colnames(expression))],e_values=unname(unlist(expression[expression$ID == e_phe,2:length(colnames(expression))])))
            for (m in 1:length(asso_others$Phenotype)){
                m_phe = asso_others$Phenotype[m]
                m_molphe = asso_others$molphe[m]
                out_file = paste0(out_dir,variant_type,"/",var,":",e_phe,":",m_phe,":mediation.input")
                if (!file.exists(out_file)){
                    if (m_molphe == "sQTL"){ m_df = data.frame(sample_id=colnames(splicing)[2:length(colnames(splicing))],m_values=unname(unlist(splicing[splicing$ID == m_phe,2:length(colnames(splicing))]))) }
                    if (m_molphe == "apaQTL"){ m_df = data.frame(sample_id=colnames(APA)[2:length(colnames(APA))],m_values=unname(unlist(APA[APA$ID == m_phe,2:length(colnames(APA))]))) }
                    if (m_molphe == "meQTL"){ m_df = data.frame(sample_id=colnames(methylation)[2:length(colnames(methylation))],m_values=unname(unlist(methylation[methylation$ID == m_phe,2:length(colnames(methylation))]))) }
                    df = left_join(var_df, sample_info, by = "sample_id") %>%
                        left_join(e_df, by = "sample_id") %>%
                        left_join(m_df, by = "sample_id")
                    fwrite(df,out_file,sep = "\t",quote = F,row.names = F)
                }
            }
        }
    }, mc.cores = mc)
}