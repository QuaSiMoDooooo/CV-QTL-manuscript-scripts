#! Rscript
# -------------
# FileName     : plot_SV_sampleSV_count
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-11-22 16:01
# Last Modified: 2024-12-14 16:05
# Modified By  : EastsunW
# -------------
# Description  : 可视化样本非冗余SV数据统计
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(gt))
suppressPackageStartupMessages(library(ggbreak))
source("https://gitee.com/eastsunw/personal_code_notebook/raw/master/plot_tools/wdy_theme.r")
custom_colors <- c(
    "#FA7F6F",
    "#8ECFC9",
    "#FFBE7A",
    "#82B0D2",
    "#BEB8DC",
    "#E7DAD2",
    "#999999"
)
setwd("/home/wangdy/Projects/Weibin/Downstreams/data_stat")

df <- fread("results/stat/Variant_stat/sampleSV_type_count.txt") %>%
    mutate(all = INS + DEL + INV + DUP)
max_y <- max(df$all)
sampleSV_barplot <- df %>%
    arrange(desc(all), desc(INS), desc(DEL), desc(INV), desc(DUP)) %>%
    select(-all) %>%
    mutate_at(vars("sample"), ~ factor(., levels = .)) %>%
    pivot_longer(-1, names_to = "SV Type", values_to = "Number of SVs") %>%
    mutate_at(vars("SV Type"), ~ factor(., levels = c("INS", "DEL", "INV", "DUP")))
max_invdup <- ceiling(sum(c(
    max(sampleSV_barplot$`Number of SVs`[sampleSV_barplot$`SV Type` == "INV"]),
    max(sampleSV_barplot$`Number of SVs`[sampleSV_barplot$`SV Type` == "DUP"])
)) / 10 + 1) * 10
sampleSV_table <- sampleSV_barplot %>%
    group_by(`SV Type`) %>%
    summarise(
        Mean = round(mean(`Number of SVs`), 0),
        Min = round(min(`Number of SVs`), 0),
        Max = round(max(`Number of SVs`), 0)
    ) %>%
    column_to_rownames("SV Type") %>%
    t() %>%
    as.data.frame() %>%
    mutate_all(~ scales::comma(.)) %>%
    rownames_to_column("Statistic")
sampleSV_table %>%
    gt(rowname_col = "Statistic") %>%
    cols_label(
        INS = html(
            fontawesome::fa("square", prefer_type = "solid", fill = custom_colors[1]),
            "INS"
        ),
        DEL = html(
            fontawesome::fa("square", prefer_type = "solid", fill = custom_colors[2]),
            "DEL"
        ),
        INV = html(
            fontawesome::fa("square", prefer_type = "solid", fill = custom_colors[3]),
            "INV"
        ),
        DUP = html(
            fontawesome::fa("square", prefer_type = "solid", fill = custom_colors[4]),
            "DUP"
        )
    ) %>%
    opt_table_outline(width = "1.5px", color = "black") %>%
    tab_options(
        column_labels.border.top.color = "black",
        column_labels.border.top.width = "1.5px",
        column_labels.border.top.style = "solid",
        table_body.border.bottom.color = "black",
        table_body.border.bottom.width = "1.5px",
        table_body.border.bottom.style = "solid",
        table.border.left.color = "black",
        table.border.left.width = "1.5px",
        table.border.left.style = "solid",
        table.border.right.color = "black",
        table.border.right.width = "1.5px",
        table.border.right.style = "solid",
        table_body.hlines.color = "gray50",
        table_body.hlines.width = "1.5px",
        table_body.hlines.style = "dotted",
        stub.border.color = "gray50",
        stub.border.width = "1.5px",
        stub.border.style = "dotted",
        table_body.border.top.color = "gray50",
        table_body.border.top.width = "1.5px",
        table_body.border.top.style = "dotted",
        column_labels.border.bottom.color = "gray50",
        column_labels.border.bottom.width = "1.5px",
        column_labels.border.bottom.style = "none"
    ) %>%
    gtsave("results/plot/sampleSV_type_stat.html", inline_css = TRUE)

pdf("results/plot/sampleSV_type_count.pdf", width = 8, height = 3, onefile = FALSE)
ggplot(
    data = sampleSV_barplot,
    mapping = aes(
        x = sample,
        y = `Number of SVs`,
        fill = `SV Type`
    )
) +
    geom_bar(
        stat = "identity",
        width = 1
    ) +
    scale_fill_manual(values = custom_colors) +
    scale_y_break(
        breaks = c(max_invdup, 2000),
        expand = expansion(0, 0),
        scales = 1,
        ticklabels = seq(2000, (ceiling(max_y / 1000) + 1) * 1000, length.out = 3)
    ) +
    scale_y_continuous(
        expand = expansion(0, 0),
        n.breaks = 3,
        limits = c(0, max_y * 1.1)
    ) +
    wdy_theme(
        base_size = 12,
        base_line_size = 1,
        legend.position = "inside",
        legend.direction = "horizontal",
        legend.justification = c(0.5, 1),
        legend.position.inside = c(0.5, 0.99),
        legend.background = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank()
    ) +
    labs(
        x = "Samples",
        fill = "SV type"
    )
dev.off()
