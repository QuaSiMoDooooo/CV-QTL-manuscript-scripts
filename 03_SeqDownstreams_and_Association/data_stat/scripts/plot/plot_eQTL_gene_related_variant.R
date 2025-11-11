#! Rscript
# -------------
# FileName     : plot_eQTL_gene_related_variant
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2025-02-27 21:47
# Last Modified: 2025-02-28 18:19
# Modified By  : EastsunW
# -------------
# Description  :
# -------------

# 加载必要的包
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
source("https://gitee.com/eastsunw/personal_code_notebook/raw/master/plot_tools/wdy_theme.r")
custom_colors1 <- c(
    "#EB6123",
    "#3DB1A6"
)
custom_colors2 <- c(
    "#0072B1",
    "#2F2585",
    "#FFAF0E",
    "#DB1048"
)
setwd("/home/wangdy/Projects/Weibin/Downstreams/data_stat")

variant_types <- c("SNP", "InDel", "MNV", "SV")

all_eQTL <- data.frame()

for (variant_type in variant_types) {
    cis_eQTL_path <- paste0("data/QTL_results/", variant_type, "_eQTL/cis.filtered.txt.gz")
    trans_eQTL_path <- paste0("data/QTL_results/", variant_type, "_eQTL/trans.filtered.txt.gz")
    if (file.exists(cis_eQTL_path)) {
        cis <- fread(cis_eQTL_path, sep = "\t") %>% mutate(variant_type = variant_type)
        all_eQTL <- rbind(all_eQTL, cis)
    }
    if (file.exists(trans_eQTL_path)) {
        trans <- fread(trans_eQTL_path, sep = "\t") %>% mutate(variant_type = variant_type)
        all_eQTL <- rbind(all_eQTL, trans)
    }
}

gene_associated <- all_eQTL %>%
    group_by(Gene, qtl_type, variant_type) %>%
    summarize(n = n()) %>%
    mutate_at(vars(variant_type), factor, levels = variant_types)

gene_asso_count <- gene_associated %>%
    # 将“n”列的计数转换为范围，1~5，5~10，>10
    mutate(range = case_when(
        n <= 1 ~ "single gene",
        TRUE ~ ">1 gene"
    )) %>%
    group_by(qtl_type, variant_type, range) %>%
    summarize(n = n()) %>%
    mutate_at(vars(range), factor, levels = c("single gene", ">1 gene"))

fwrite(gene_asso_count, "results/stat/QTL/eQTL_gene_related_variant_count.txt", sep = "\t")
fwrite(gene_associated, "results/stat/QTL/eQTL_gene_related_variant.txt", sep = "\t")

plot_count <- ggplot(
    data = gene_asso_count,
    aes(x = variant_type, y = n)
) +
    geom_bar(
        aes(fill = range),
        stat = "identity",
        position = "dodge"
    ) +
    scale_fill_manual(values = custom_colors1) +
    facet_wrap(~qtl_type, nrow = 1, scales = "free_y") +
    labs(
        x = "Variant Type",
        y = "Number of Genes",
        fill = "Associated gene"
    ) +
    wdy_theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "right",
        panel.spacing.x = unit(10, "pt"),
        strip.background = element_rect(color = NA, fill = "gray80")
    )

ggsave(
    "results/plot/eQTL_gene_related_variant_bar.pdf",
    plot_count,
    width = 7,
    height = 3
)


# gene_associated %>%
#     group_by(Gene) %>%
#     summarise(n = sum(n)) %>%
#     filter(n > 100 & n < 300) %>%
#     arrange(desc(n))
# all_eQTL %>%
#     filter(-log10(FDR) < 10) %>%
#     group_by(Gene) %>%
#     summarise(n = n()) %>%
#     filter(n > 100 & n < 300) %>%
#     arrange(desc(n)) %>%
#     filter(Gene == "MTARC1")


most_associated_gene <- "MTARC1"
gene_range <- c(22072e4, 22082e4)

most_associated_gene_info <- all_eQTL %>%
    filter(Gene %in% most_associated_gene) %>%
    filter(qtl_type == "cis") %>%
    separate(Variant_position, into = c("Chr", "Range"), sep = ":") %>%
    separate(Range, into = c("Start", "End"), sep = "-") %>%
    select(Variant, Chr, Start, End, FDR, variant_type) %>%
    mutate_at(vars(Start, End), as.numeric)


# 曼哈顿图
manhaten <- ggplot(
    data = most_associated_gene_info,
    mapping = aes(
        x = Start,
        y = -log10(FDR),
        color = variant_type
    )
) +
    geom_point(
        size = 2,
        alpha = 0.8,
        shape = 16
    ) +
    geom_hline(
        yintercept = -log10(5e-8),
        linetype = "dashed",
        color = "red"
    ) +
    # 将虚线下面的内容填充为灰色
    annotate(
        "rect",
        xmin = gene_range[1],
        xmax = gene_range[2],
        ymin = 0,
        ymax = -log10(5e-8),
        fill = "gray80",
        alpha = 0.5
    ) +
    scale_x_continuous(
        limits = gene_range,
        labels = scales::unit_format(unit = "Mb", scale = 1e-6, sep = "", big.mark = ","),
        n.breaks = 4
    ) +
    scale_y_continuous(
        limits = c(0, 15)
    ) +
    scale_color_manual(values = custom_colors2) +
    labs(
        x = paste0("QTL locus for ", most_associated_gene, " gene"),
        y = "-log10(FDR)",
        color = "Variant Type"
    ) +
    guides(
        color = guide_legend(override.aes = list(size = 5))
    ) +
    wdy_theme(
        legend.position = "right",
        panel.spacing = unit(5, "points"),
        strip.background = element_rect(color = NA, fill = "gray80")
    )
ggsave(
    "results/plot/eQTL_gene_related_variant_manhattan.pdf",
    manhaten,
    width = 7,
    height = 3
)
