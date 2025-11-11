#! Rscript
# -------------
# FileName     : clinic_phenotype_cor
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2025-03-04 14:39
# Last Modified: 2025-03-07 17:15
# Modified By  : EastsunW
# -------------
# Description  : 将生化指标的分组和分子表型进行关联
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(doParallel))
source("https://gitee.com/eastsunw/personal_code_notebook/raw/master/plot_tools/wdy_theme.r")
setwd("/home/wangdy/Projects/Weibin/Downstreams/variant_annotation")

bio_dir <- "data/Clinic"
bio_info_path <- file.path(bio_dir, "bio_info.txt")
bio_level_path <- file.path(bio_dir, "bio_level.txt")
bio_risk_path <- file.path(bio_dir, "bio_risk.txt")

pheno_dir <- "~/Projects/Weibin/Downstreams/phenotype_qtl/results/common_data/phenotype"
pheno_expression_path <- file.path(pheno_dir, "expression.quantity.filtered.txt")
pheno_apa_path <- file.path(pheno_dir, "APA.quantity.filtered.txt")
pheno_as_path <- file.path(pheno_dir, "splicing.quantity.filtered.txt")
pheno_methylation_path <- file.path(pheno_dir, "methylation.quantity.filtered.txt")

cor_out_dir <- "results/network/correlation"
cor_expression_path <- file.path(cor_out_dir, "correlation_expression.txt")
cor_apa_path <- file.path(cor_out_dir, "correlation_apa.txt")
cor_as_path <- file.path(cor_out_dir, "correlation_as.txt")
cor_methylation_path <- file.path(cor_out_dir, "correlation_methylation.txt")

level_mapping <- c(
    "过低" = "-3",
    "低" = "-2",
    "偏低" = "-1",
    "正常" = "0",
    "偏高" = "1",
    "高" = "2",
    "过高" = "3"
)
risk_mapping <- c(
    "正常" = "0",
    "低风险" = "1",
    "中风险" = "2",
    "高风险" = "3"
)

# 导入数据
bio_info <- fread(bio_info_path, sep = "\t") %>% select(Indicator, Category)
bio_level <- fread(bio_level_path, sep = "\t") %>%
    column_to_rownames("Indicator") %>%
    filter(!if_all(.cols = everything(), .fns = is.na)) %>%
    rownames_to_column("Indicator") %>%
    pivot_longer(cols = -Indicator, names_to = "Sample", values_to = "Level") %>%
    group_by(Indicator) %>%
    filter(!Level %in% c("阴性", "阳性") & !is.na(Level)) %>%
    mutate(n_level = length(table(Level))) %>%
    filter(
        !is.na(Level),
        n_level > 1
    ) %>%
    select(-n_level) %>%
    mutate(Level = recode(Level, !!!level_mapping)) %>%
    pivot_wider(names_from = "Sample", values_from = "Level")
bio_risk <- fread(bio_risk_path, sep = "\t") %>%
    column_to_rownames("Indicator") %>%
    filter(!if_all(.cols = everything(), .fns = is.na)) %>%
    rownames_to_column("Indicator") %>%
    pivot_longer(cols = -Indicator, names_to = "Sample", values_to = "Risk") %>%
    group_by(Indicator) %>%
    mutate(n_level = length(table(Risk))) %>%
    filter(
        !is.na(Risk),
        n_level > 1
    ) %>%
    select(-n_level) %>%
    mutate(Risk = recode(Risk, !!!risk_mapping)) %>%
    pivot_wider(names_from = "Sample", values_from = "Risk")

expression <- fread(pheno_expression_path, sep = "\t")
apa <- fread(pheno_apa_path, sep = "\t")
as <- fread(pheno_as_path, sep = "\t")
methylation <- fread(pheno_methylation_path, sep = "\t")

# 计算相关性矩阵
calculate_correlation <- function(indicator, phenotype, bio_df, pheno_df) {
    indicator_df <- bio_df %>%
        filter(Indicator == indicator) %>%
        pivot_longer(
            cols = -1,
            names_to = "Sample",
            values_to = "BioValue"
        )
    phenotype_df <- pheno_df %>%
        filter(ID == phenotype) %>%
        pivot_longer(
            cols = -1,
            names_to = "Sample",
            values_to = "PhenoValue"
        )
    cor_df <- inner_join(
        x = indicator_df,
        y = phenotype_df,
        by = "Sample"
    ) %>%
        filter(!is.na(BioValue) & !is.na(PhenoValue))
    output_obj <- list(list(
        indicator = indicator,
        phenotype = phenotype,
        data = cor_df,
        cor = cor(
            x = as.numeric(cor_df$BioValue),
            y = as.numeric(cor_df$PhenoValue),
            method = "spearman"
        )
    ))
    names(output_obj) <- paste0(c(indicator, phenotype), collapse = "_")
    return(output_obj)
}

plot_correation <- function(indicator, phenotype, data) {
    ggplot(
        data = data,
        mapping = aes(
            x = BioValue,
            y = PhenoValue
        )
    ) +
        geom_boxplot() +
        geom_jitter(width = 0.2, height = 0, alpha = 0.5) +
        geom_smooth(method = "lm", se = FALSE) +
        labs(
            x = indicator,
            y = phenotype
        ) +
        theme_wdy()
}

cor_list_expression <- list()
all_indicators <- unique(bio_level$Indicator)
all_phenotypes <- unique(expression$ID)
cor_list_expression <- foreach(indicator = all_indicators, .combine = "c") %:%
    foreach(phenotype = all_phenotypes, .combine = "c") %dopar% {
        calculate_correlation(
            indicator = indicator,
            phenotype = phenotype,
            bio_df = bio_level,
            pheno_df = expression
        )
    }
# 将列表项中除了data的部分进行合并，indicator是行，phenotype是列，cor是值
cor_matrix_expression <- do.call(
    what = rbind,
    args = lapply(cor_list_expression, function(x) {
        x[c("indicator", "phenotype", "cor")]
    })
) %>%
    as.data.frame() %>%
    rename(
        Indicator = indicator,
        Phenotype = phenotype,
        Correlation = cor
    ) %>%
    pivot_wider(names_from = "phenotype", values_from = "cor") %>%
    column_to_rownames("indicator") %>%
    mutate_all(as.numeric)
fwrite(cor_matrix_expression %>% rownames_to_column("indicator"), file = cor_expression_path, sep = "\t")

# 找到相关性最高的指标以及对应的基因，返回三列，Indicator, max_cor, max_cor_phenotype
max_cor_expression <- cor_matrix_expression %>%
    rownames_to_column("Indicator") %>%
    pivot_longer(cols = -Indicator, names_to = "Phenotype", values_to = "Cor") %>%
    group_by(Indicator) %>%
    filter(Cor == max(Cor)) %>%
    ungroup() %>%
    rename(max_cor_phenotype = Phenotype) %>%
    select(Indicator, Cor, max_cor_phenotype) %>%
    arrange(desc(Cor))

# 计算并保存相关性矩阵
cor_expression <- calculate_correlation(bio_level, expression)
# cor_apa <- calculate_correlation(bio_level, apa)
# cor_as <- calculate_correlation(bio_level, as)
# cor_methylation <- calculate_correlation(bio_level, methylation)

write.table(cor_expression, file = "correlation_expression.txt", sep = "\t", quote = FALSE)
# write.table(cor_apa, file = "correlation_apa.txt", sep = "\t", quote = FALSE)
# write.table(cor_as, file = "correlation_as.txt", sep = "\t", quote = FALSE)
# write.table(cor_methylation, file = "correlation_methylation.txt", sep = "\t", quote = FALSE)
