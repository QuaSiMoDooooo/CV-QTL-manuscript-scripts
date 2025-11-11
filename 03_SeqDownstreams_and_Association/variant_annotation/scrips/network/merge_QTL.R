#! Rscript
# -------------
# FileName     : merge_QTL
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2025-03-07 09:15
# Last Modified: 2025-03-07 17:10
# Modified By  : EastsunW
# -------------
# Description  : 合并cis和trans的QTL关联结果
# -------------

# 加载必要的包
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(doParallel))
setwd("/home/wangdy/Projects/Weibin/Downstreams/variant_annotation")

qtl_in_dir <- "data/QTL_results"
qtl_out_dir <- "results/network/QTL"

# gwas_in_dir <- "data/GWAS"
# gwas_out_dir <- "results/network/GWAS"

# cor_qtl_in_dir <- "data/QTL_data"
# cor_clinic_in_dir <- "data/Clinic"
# cor_out_dir <- "results/network/correlation"

variant_types <- c("SNP", "InDel", "MNV", "SV")
pheno_types <- c("eQTL", "apaQTL", "sQTL", "meQTL")
qtl_types <- c("cis", "trans")

# 遍历所有的QTL类型，合并cis和trans的结果
registerDoParallel(cores = 32)

foreach(variant_type = variant_types, .combine = "c") %:%
    foreach(pheno_type = pheno_types, .combine = "c") %dopar% {
        qtl_results <- list(
            cis = fread(file.path(
                qtl_in_dir,
                paste0(variant_type, "_", pheno_type),
                "cis.filtered.txt.gz"
            )),
            trans = fread(file.path(
                qtl_in_dir,
                paste0(variant_type, "_", pheno_type),
                "trans.filtered.txt.gz"
            ))
        )
        bind_rows(qtl_results) %>%
            mutate(
                variant_type = variant_type,
                pheno_type = pheno_type
            ) %>%
            fwrite(
                file.path(
                    qtl_out_dir,
                    paste0(variant_type, "_", pheno_type, ".txt")
                ),
                sep = "\t",
                quote = FALSE
            )
    }

stopImplicitCluster()
