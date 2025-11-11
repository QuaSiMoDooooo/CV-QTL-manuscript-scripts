#! Rscript
# -------------
# FileName     : stat_UCSC_qtl
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2025-01-03 22:47
# Last Modified: 2025-03-06 09:22
# Modified By  : EastsunW
# -------------
# Description  :
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
source("https://gitee.com/eastsunw/personal_code_notebook/raw/master/plot_tools/wdy_theme.r")
custom_colors <- c(
    "#FA7F6F",
    "#8ECFC9",
    "#FFBE7A",
    "#82B0D2",
    "#BEB8DC",
    "#e9e9e9",
    "#999999"
)
setwd("/home/wangdy/Projects/Weibin/Downstreams/variant_annotation")

data_set <- "UCSC"
data_dir <- paste0("results/epi_enrichment/", data_set)
out_dir <- paste0("results/epi_enrichment/", data_set, "/stat_qtl")
if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}

annotations <- c(
    "Promoter",
    "Exon",
    "Intron",
    "3UTR",
    "5UTR",
    "Enhancer",
    "CpG",
    "Repeat",
    "Intergenic"
)
variant_types <- c("SNP", "InDel", "MNV", "SV")
qtl_types <- c("eQTL", "apaQTL", "sQTL", "meQTL")

merged_enrich <- list()
for (variant_type in variant_types) {
    for (qtl_type in qtl_types) {
        merged_data <- data.frame(annotation = annotations)
        result_path <- file.path(
            data_dir,
            "qtl",
            paste0(
                data_set, "_",
                variant_type, "-", qtl_type,
                ".results.txt"
            )
        )
        if (file.exists(result_path)) {
            result <- fread(result_path, sep = "\t") %>%
                select(annotation, l2fold, qvalue) %>%
                unite(l2fold, qvalue, col = "stat", sep = "__")
            merged_data <<- left_join(
                merged_data,
                result,
                by = "annotation"
            )
            merged_enrich[[paste0(variant_type, "-", qtl_type)]] <- merged_data
        }
    }
}

plot_df <- data.frame()
temp <- sapply(
    X = names(merged_enrich),
    FUN = function(name) {
        data <- merged_enrich[[name]] %>%
            mutate(qtl_type = name) %>%
            separate(qtl_type, c("variant_type", "qtl_type"), sep = "-") %>%
            separate(stat, c("log2FC", "FDR"), sep = "__") %>%
            mutate_at(vars(log2FC, FDR), as.numeric)
        plot_df <<- rbind(plot_df, data)
    }
)
rm(temp)

plot_df_adj <- plot_df %>%
    mutate(p_sig = case_when(
        FDR >= 0.05 ~ NA,
        FDR < 0.05 ~ "*",
    )) %>%
    relocate(variant_type, qtl_type, .after = annotation) %>%
    mutate(
        annotation = factor(annotation, levels = annotations),
        variant_type = factor(variant_type, levels = variant_types),
        qtl_type = factor(qtl_type, levels = qtl_types)
    )

fwrite(
    x = plot_df_adj,
    file = file.path(out_dir, paste0("enrich_", data_set, "_qtl.txt")),
    sep = "\t",
    quote = FALSE
)

plot_by_qtl <- ggplot(
    data = plot_df_adj,
    mapping = aes(
        x = annotation,
        y = variant_type,
        fill = log2FC
    )
) +
    geom_tile() +
    scale_fill_gradientn(
        limits = range(-3, 3),
        oob = scales::squish,
        labels = c("-3", "", "", "0", "", "", "3"),
        breaks = c(-3, -2, -1, 0, 1, 2, 3),
        colors = c("#122A64", "#223290", "#3085BC", "#ffffff", "#F56B46", "#B8070D", "#670519"),
        na.value = "#ffffff"
    ) +
    geom_text(
        aes(label = p_sig),
        size = 10 / .pt,
        color = "black",
        hjust = 0.5,
        vjust = 0.5
    ) +
    facet_wrap(
        ~qtl_type,
        nrow = 1,
        strip.position = "top"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
        x = "Annotation",
        y = "Variant type",
        fill = "log2FC"
    ) +
    coord_flip() +
    wdy_theme(
        base_size = 10,
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.text.y = element_text(angle = 0),
        panel.spacing.x = unit(2, "pt"),
        panel.spacing.y = unit(2, "pt"),
        strip.background = element_rect(fill = "#ebf1f5", linewidth = NA)
    )

plot_by_variant <- ggplot(
    data = plot_df_adj,
    mapping = aes(
        x = annotation,
        y = qtl_type,
        fill = log2FC
    )
) +
    geom_tile() +
    scale_fill_gradientn(
        limits = range(-3, 3),
        oob = scales::squish,
        labels = c("-3", "", "", "0", "", "", "3"),
        breaks = c(-3, -2, -1, 0, 1, 2, 3),
        colors = c("#122A64", "#223290", "#3085BC", "#ffffff", "#F56B46", "#B8070D", "#670519"),
        na.value = "#ffffff"
    ) +
    geom_text(
        aes(label = p_sig),
        size = 10 / .pt,
        color = "black",
        hjust = 0.5,
        vjust = 0.5
    ) +
    facet_wrap(
        ~variant_type,
        nrow = 1,
        strip.position = "top"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
        x = "Annotation",
        y = "Variant type",
        fill = "log2FC"
    ) +
    coord_flip() +
    wdy_theme(
        base_size = 10,
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.text.y = element_text(angle = 0),
        panel.spacing.x = unit(2, "pt"),
        panel.spacing.y = unit(2, "pt"),
        strip.background = element_rect(fill = "#ebf1f5", linewidth = NA)
    )

ggsave(
    file = file.path(out_dir, "plot_by_qtl.pdf"),
    plot = plot_by_qtl,
    width = 5,
    height = 2.5
)
ggsave(
    file = file.path(out_dir, "plot_by_variant.pdf"),
    plot = plot_by_variant,
    width = 5,
    height = 2.5
)
