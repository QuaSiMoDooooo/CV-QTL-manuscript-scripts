#! Rscript

# -------------
# FileName     : 05_calc_SV_QTL_tag_rate.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Calculate SV QTL tag rate and SNP QTL overlap analysis
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
setwd("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/05_SV_QTL_SNP_tag")

rm(list = ls())
ld = fread("01_ld_sv_snp.tsv", sep = "\t")
ld = ld[abs(ld$R2) > 0.75,]

library(data.table)
library(stringr)

compute_qtl_summary <- function(mol_type_vec = c("eQTL","sQTL","apaQTL","meQTL"),
                                reg_type = "cis",
                                base_path = "/home/wangdy/Projects/Weibin/Downstreams/phenotype_qtl/results",
                                ld_df) {
  mol_map <- list(
    eQTL = "expression",
    sQTL = "splicing",
    apaQTL = "APA",
    meQTL = "methylation"
  )
  
  sv_qtl_tag_rate = c()
  tag_snp_qtl_rate = c()
  same_level_rate = c()
  more_level_rate = c()
  less_level_rate = c()
  no_level_rate = c()
  
  for (mol_type in mol_type_vec) {
    mol_dir <- mol_map[[mol_type]]
    snp_dir <- file.path(base_path, paste0("SNP-", mol_dir), "QTL_results")
    sv_dir <- file.path(base_path, paste0("SV-", mol_dir), "QTL_results")
    
    get_qtl_file <- function(dir) {
      if (reg_type == "all") {
        f1 <- fread(file.path(dir, "cis.filtered.txt.gz"))
        f2 <- fread(file.path(dir, "trans.filtered.txt.gz"))
        rbind(f1, f2)
      } else {
        fread(file.path(dir, paste0(reg_type, ".filtered.txt.gz")))
      }
    }
    
    snp_qtl <- get_qtl_file(snp_dir)
    sv_qtl  <- get_qtl_file(sv_dir)
    
    tmp = str_split(snp_qtl$Variant, "_", simplify = T)
    tmp[,1] = sub("chr", "", tmp[,1])
    snp_qtl$SNP_B = paste(tmp[,1], tmp[,2], "SNP", sep = "_")
    
    sv_ids_uniq = unique(sv_qtl$Variant)
    tag_sv_ids = sv_ids_uniq[sv_ids_uniq %in% ld_df$sv_id]
    sv_qtl_tag_rate = c(sv_qtl_tag_rate, length(tag_sv_ids) / length(sv_ids_uniq))
    
    tag_snp_qtl_num = 0
    for (sv in tag_sv_ids) {
      tag_snp_list = ld_df$SNP_B[ld_df$sv_id == sv]
      if (any(tag_snp_list %in% snp_qtl$SNP_B)) {
        tag_snp_qtl_num = tag_snp_qtl_num + 1
      }
    }
    tag_snp_qtl_rate = c(tag_snp_qtl_rate, tag_snp_qtl_num / length(tag_sv_ids))
    
    same_level = 0; more_level = 0; less_level = 0; no_level = 0; eff_sv_num = 0
    sv_qtl_tagged = sv_qtl[Variant %in% tag_sv_ids, ]
    
    for (sv in unique(sv_qtl_tagged$Variant)) {
      sv_phe_list = sv_qtl_tagged$Phenotype[sv_qtl_tagged$Variant == sv]
      tag_snp_list = ld_df$SNP_B[ld_df$sv_id == sv]
      tag_snp_qtl_list = tag_snp_list[tag_snp_list %in% snp_qtl$SNP_B]
      
      if (length(tag_snp_qtl_list) >= 1) {
        eff_sv_num = eff_sv_num + 1
        sv_phe = unique(sv_phe_list)
        tag_snp_qtl_phe = unique(snp_qtl$Phenotype[snp_qtl$SNP_B %in% tag_snp_qtl_list])
        tag_snp_qtl_phe_num_in_sv_phe = sum(tag_snp_qtl_phe %in% sv_phe)
        
        if (tag_snp_qtl_phe_num_in_sv_phe == length(sv_phe)) {
          if (length(tag_snp_qtl_phe) == length(sv_phe)) {
            same_level = same_level + 1
          } else {
            more_level = more_level + 1
          }
        } else {
          if (tag_snp_qtl_phe_num_in_sv_phe == 0) {
            no_level = no_level + 1
          } else {
            less_level = less_level + 1
          }
        }
      }
    }
    
    same_level_rate = c(same_level_rate, same_level/eff_sv_num)
    more_level_rate = c(more_level_rate, more_level/eff_sv_num)
    less_level_rate = c(less_level_rate, less_level/eff_sv_num)
    no_level_rate = c(no_level_rate, no_level/eff_sv_num)
  }
  
  result_df = data.frame(
    qtl = mol_type_vec,
    sv_qtl_tag_rate,
    tag_snp_qtl_rate,
    same_level = same_level_rate,
    more_level = more_level_rate,
    less_level = less_level_rate,
    no_level = no_level_rate
  )
  return(result_df)
}

res_cis = compute_qtl_summary(mol_type_vec = c("eQTL","sQTL","apaQTL","meQTL"),
                              reg_type = "cis",
                              base_path = "/home/wangdy/Projects/Weibin/Downstreams/phenotype_qtl/results",
                              ld_df = ld)
fwrite(res_cis, "01_stats_cis.tsv", sep = "\t")

parallel::mclapply(c("cis","trans","all"), function(reg_type) {
  outpath = paste0("01_stats_", reg_type, ".tsv")
  if (!file.exists(outpath)) {
    res = compute_qtl_summary(mol_type_vec = c("eQTL","sQTL","apaQTL","meQTL"),
                              reg_type = reg_type,
                              base_path = "/home/wangdy/Projects/Weibin/Downstreams/phenotype_qtl/results",
                              ld_df = ld)
    fwrite(res, outpath, sep = "\t")
  }
}, mc.cores = 4)
