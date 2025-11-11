#! Rscript
# -------------
# FileName     : calculate_density
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-12-08 18:47
# Last Modified: 2024-12-17 16:04
# Modified By  : EastsunW
# -------------
# Description  :
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(parallel))
setwd("/home/wangdy/Projects/Weibin/Downstreams/data_stat")

exclude_chr <- c("chrY", "chrM")

snp_bed <- fread(
    "results/bed/SNP.all.bed",
    header = FALSE
) %>%
    filter(!V1 %in% exclude_chr)
indel_bed <- fread(
    "results/bed/InDel.all.bed",
    header = FALSE
) %>%
    filter(!V1 %in% exclude_chr)
mnv_bed <- fread(
    "results/bed/MNV.all.bed",
    header = FALSE
) %>%
    filter(!V1 %in% exclude_chr)
sv_all_bed <- fread(
    "results/bed/SV.all.bed",
    header = FALSE
) %>%
    filter(!V1 %in% exclude_chr)
tr_bed <- fread(
    "results/bed/TR.normal.bed",
    header = FALSE
) %>%
    filter(!V1 %in% exclude_chr)
repeat_bed <- fread(
    "data/repeat_hg38.bed",
    header = FALSE
) %>%
    filter(!V1 %in% exclude_chr)
gene_bed <- fread(
    "data/gene_hg38.bed",
    header = FALSE
) %>%
    filter(!V1 %in% exclude_chr)

calculate_density_thread <- function(chr, bed) {
    bed %>%
        filter(.[[1]] == chr) %>%
        genomicDensity(
            overlap = FALSE,
            count_by = "number",
            window.size = 5e5
        ) %>%
        setnames(c("chr", "start", "end", "value")) %>%
        group_by(chr) %>%
        mutate(
            min_val = min(value, na.rm = TRUE),
            max_val = max(value, na.rm = TRUE)
        ) %>%
        # 进行归一化处理
        mutate(value = (value - min_val) / (max_val - min_val)) %>%
        # 选择需要的列，这里我们保留所有列，如果只需要归一化后的value列，可以调整为select(value)
        select(-min_val, -max_val)
}

# 多线程
chr_list <- paste0("chr", c(as.character(1:22), "X"))
bed_list <- list(
    "SNP"    = snp_bed,
    "InDel"  = indel_bed,
    "MNV"    = mnv_bed,
    "SV_all" = sv_all_bed,
    "SV_INS" = sv_all_bed %>% filter(V6 == "INS"),
    "SV_DEL" = sv_all_bed %>% filter(V6 == "DEL"),
    "SV_DUP" = sv_all_bed %>% filter(V6 == "DUP"),
    "SV_INV" = sv_all_bed %>% filter(V6 == "INV"),
    "gene"   = gene_bed,
    "repeat" = repeat_bed,
    "TR"     = tr_bed
)


density_list <- lapply(
    bed_list,
    function(bed) {
        do.call(
            rbind,
            mclapply(
                X = chr_list,
                FUN = calculate_density_thread,
                bed = bed,
                mc.cores = 23
            )
        )
    }
)

lapply(names(density_list), function(name) {
    fwrite(density_list[[name]], paste0("results/stat/circos_density/", name, "_density.txt"), sep = "\t")
})
