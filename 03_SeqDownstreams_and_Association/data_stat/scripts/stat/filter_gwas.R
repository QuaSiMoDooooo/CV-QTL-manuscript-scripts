#! Rscript
# -------------
# FileName     : filter_gwas
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2025-03-04 18:53
# Last Modified: 2025-06-04 09:59
# Modified By  : EastsunW
# -------------
# Description  : 给不同指标的GWAS结果进行FDR计算并过滤
# -------------

# 加载必要的包
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))

setwd("/home/wangdy/Projects/Weibin/Downstreams/data_stat")
gwas_dir <- "data/GWAS"
gwas_output_dir <- "results/stat/GWAS"
variant_types <- c("SNP", "InDel", "MNV", "SV")

filter_gwas <- function(df, fdr_method = "fdr", fdr_threshold = 0.05, p_col = "P") {
    df %>%
        filter(TEST == "ADD") %>%
        filter_at(vars(all_of(p_col)), ~ !is.na(.)) %>%
        mutate_at(vars(all_of(p_col)), as.numeric) %>%
        mutate(FDR = p.adjust(.data[[p_col]], method = fdr_method)) %>%
        filter(FDR < fdr_threshold)
}

for (variant_type in variant_types) {
    gwas_path <- file.path(gwas_dir, paste0(variant_type, "_GWAS"))
    if (file.exists(gwas_path)) {
        merged_gwas <- data.frame()
        dir.create(file.path(gwas_output_dir, variant_type), showWarnings = FALSE)
        for (indicator_type in c("common", "male", "female")) {
            result_path <- file.path(gwas_path, paste0("biochem_", indicator_type))
            if (file.exists(result_path)) {
                # 列出目录下所有的.assoc.linear结尾的文件，并且从文件名提取出指标名，加到结果中
                gwas_files <- list.files(result_path, pattern = "\\.glm.linear$", full.names = FALSE)
                temp <- mclapply(gwas_files, function(file) {
                    indicator_name <- strsplit(file, ".glm.linear")[[1]][1]
                    gwas_df <- fread(file.path(result_path, file), showProgress = FALSE)
                    gwas_filtered <- filter_gwas(gwas_df) %>%
                        mutate(
                            Indicator = indicator_name,
                            IndicatorType = indicator_type,
                            variant_type = variant_type
                        )
                    if (nrow(gwas_filtered) > 0) {
                        fwrite(
                            gwas_filtered,
                            file.path(
                                gwas_output_dir,
                                variant_type,
                                paste0(indicator_name, ".filtered.tsv")
                            ),
                            sep = "\t",
                            quote = FALSE
                        )
                        return(gwas_filtered)
                    } else {
                        return(NULL)
                    }
                }, mc.cores = 20)
                merged_gwas <- rbind(merged_gwas, as.data.frame(do.call(rbind, temp)))
            }
        }
        fwrite(
            merged_gwas,
            file.path(
                gwas_output_dir,
                paste0(variant_type, ".merged.filtered.tsv")
            ),
            sep = "\t",
            quote = FALSE
        )
    }
}
