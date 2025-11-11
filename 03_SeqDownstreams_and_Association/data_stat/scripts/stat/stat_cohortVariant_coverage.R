#! Rscript
# -------------
# FileName     : stat_cohortVariant_coverage
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-11-26 15:22
# Last Modified: 2025-03-06 21:57
# Modified By  : EastsunW
# -------------
# Description  : 统计队列的各种变异在全基因组的覆盖度，用来说明还有很多区域的功能未被解析
# -------------

# 加载必要的包
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

snp_all_stat <- fread(
    "results/bed/SNP.all.merged.bed",
    col.names = c("chr", "start", "end")
) %>%
    mutate(coverage = end - start) %>%
    select(-chr, -start, -end) %>%
    mutate(variant_type = "SNP", freq = "all") %>%
    group_by(variant_type, freq) %>%
    summarise(coverage = sum(coverage)) %>%
    mutate(percent = coverage / 3088269832)
snp_common_stat <- fread(
    "results/bed/SNP.common.merged.bed",
    col.names = c("chr", "start", "end")
) %>%
    mutate(coverage = end - start) %>%
    select(-chr, -start, -end) %>%
    mutate(variant_type = "SNP", freq = "common") %>%
    group_by(variant_type, freq) %>%
    summarise(coverage = sum(coverage)) %>%
    mutate(percent = coverage / 3088269832)
snp_rare_stat <- fread(
    "results/bed/SNP.rare.merged.bed",
    col.names = c("chr", "start", "end")
) %>%
    mutate(coverage = end - start) %>%
    select(-chr, -start, -end) %>%
    mutate(variant_type = "SNP", freq = "rare") %>%
    group_by(variant_type, freq) %>%
    summarise(coverage = sum(coverage)) %>%
    mutate(percent = coverage / 3088269832)

indel_all_stat <- fread("results/bed/InDel.all.merged.bed", col.names = c("chr", "start", "end")) %>%
    mutate(coverage = end - start) %>%
    select(-chr, -start, -end) %>%
    mutate(variant_type = "InDel", freq = "all") %>%
    group_by(variant_type, freq) %>%
    summarise(coverage = sum(coverage)) %>%
    mutate(percent = coverage / 3088269832)
indel_common_stat <- fread("results/bed/InDel.common.merged.bed", col.names = c("chr", "start", "end")) %>%
    mutate(coverage = end - start) %>%
    select(-chr, -start, -end) %>%
    mutate(variant_type = "InDel", freq = "common") %>%
    group_by(variant_type, freq) %>%
    summarise(coverage = sum(coverage)) %>%
    mutate(percent = coverage / 3088269832)
indel_rare_stat <- fread("results/bed/InDel.rare.merged.bed", col.names = c("chr", "start", "end")) %>%
    mutate(coverage = end - start) %>%
    select(-chr, -start, -end) %>%
    mutate(variant_type = "InDel", freq = "rare") %>%
    group_by(variant_type, freq) %>%
    summarise(coverage = sum(coverage)) %>%
    mutate(percent = coverage / 3088269832)

sv_all_stat <- fread("results/bed/SV.all.merged.bed", col.names = c("chr", "start", "end")) %>%
    mutate(coverage = end - start) %>%
    select(-chr, -start, -end) %>%
    mutate(variant_type = "SV", freq = "all") %>%
    group_by(variant_type, freq) %>%
    summarise(coverage = sum(coverage)) %>%
    mutate(percent = coverage / 3088269832)
sv_common_stat <- fread("results/bed/SV.common.merged.bed", col.names = c("chr", "start", "end")) %>%
    mutate(coverage = end - start) %>%
    select(-chr, -start, -end) %>%
    mutate(variant_type = "SV", freq = "common") %>%
    group_by(variant_type, freq) %>%
    summarise(coverage = sum(coverage)) %>%
    mutate(percent = coverage / 3088269832)
sv_rare_stat <- fread("results/bed/SV.rare.merged.bed", col.names = c("chr", "start", "end")) %>%
    mutate(coverage = end - start) %>%
    select(-chr, -start, -end) %>%
    mutate(variant_type = "SV", freq = "rare") %>%
    group_by(variant_type, freq) %>%
    summarise(coverage = sum(coverage)) %>%
    mutate(percent = coverage / 3088269832)

mnv_all_stat <- fread("results/bed/MNV.all.merged.bed", col.names = c("chr", "start", "end")) %>%
    mutate(coverage = end - start) %>%
    select(-chr, -start, -end) %>%
    mutate(variant_type = "MNV", freq = "all") %>%
    group_by(variant_type, freq) %>%
    summarise(coverage = sum(coverage)) %>%
    mutate(percent = coverage / 3088269832)
mnv_common_stat <- fread("results/bed/MNV.common.merged.bed", col.names = c("chr", "start", "end")) %>%
    mutate(coverage = end - start) %>%
    select(-chr, -start, -end) %>%
    mutate(variant_type = "MNV", freq = "common") %>%
    group_by(variant_type, freq) %>%
    summarise(coverage = sum(coverage)) %>%
    mutate(percent = coverage / 3088269832)
mnv_rare_stat <- fread("results/bed/MNV.rare.merged.bed", col.names = c("chr", "start", "end")) %>%
    mutate(coverage = end - start) %>%
    select(-chr, -start, -end) %>%
    mutate(variant_type = "MNV", freq = "rare") %>%
    group_by(variant_type, freq) %>%
    summarise(coverage = sum(coverage)) %>%
    mutate(percent = coverage / 3088269832)


merged_df <- rbind(
    snp_all_stat,
    snp_common_stat,
    snp_rare_stat,
    indel_all_stat,
    indel_common_stat,
    indel_rare_stat,
    mnv_all_stat,
    mnv_common_stat,
    mnv_rare_stat,
    sv_all_stat,
    sv_common_stat,
    sv_rare_stat
) %>%
    mutate(
        percent_str = paste0(signif(percent * 100, 3), "%")
    ) %>%
    select(-percent) %>%
    rename(Frequancy = freq, percent = percent_str, Variant = variant_type) %>%
    group_by(Variant)

fwrite(
    merged_df,
    "results/stat/Variant_stat/cohortVariant_coverage.txt",
    sep = "\t",
    quote = FALSE
)
