#! Rscript
# -------------
# FileName     : plot_sv_compare
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-12-27 11:00
# Last Modified: 2025-06-18 22:47
# Modified By  : EastsunW
# -------------
# Description  :
# -------------

# 加载必要的包
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
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
setwd("/home/wangdy/Projects/Weibin/Downstreams/variant_annotation")

data_sources <- c(
    "gnomAD",
    "HGDP",
    "HGSVC",
    "Iceland",
    "LRS15",
    "Han945"
)

# 读取未病数据的bed
weibin_ins <- fread(
    "results/SV_compare/Weibin/converted/Weibin.ins.bed",
    sep = "\t",
    header = FALSE,
    col.names = c("CHR", "START", "END", "ID")
)
weibin_del <- fread(
    "results/SV_compare/Weibin/converted/Weibin.del.bed",
    sep = "\t",
    header = FALSE,
    col.names = c("CHR", "START", "END", "ID")
)
weibin_inv <- fread(
    "results/SV_compare/Weibin/converted/Weibin.inv.bed",
    sep = "\t",
    header = FALSE,
    col.names = c("CHR", "START", "END", "ID")
)
weibin_dup <- fread(
    "results/SV_compare/Weibin/converted/Weibin.dup.bed",
    sep = "\t",
    header = FALSE,
    col.names = c("CHR", "START", "END", "ID")
)
weibin_all <- rbind(weibin_ins, weibin_del, weibin_inv, weibin_dup)

overlap_stat <- data.frame(
    ID = weibin_all$ID,
    sv_type = substr(weibin_all$ID, 4, 6)
)
for (source in data_sources) {
    source_df <- data.frame(
        ID = character(),
        overlap = character()
    ) %>%
        rename(!!source := overlap)
    for (svtype in c("ins", "del", "inv", "dup")) {
        overlap_path <- file.path("results/SV_compare", source, "result", paste0("overlap_", svtype, ".bed"))
        if (file.exists(overlap_path)) {
            overlap_sv <- fread(
                overlap_path,
                sep = "\t",
                header = FALSE,
                col.names = c(
                    "CHR1", "START1", "END1", "ID1",
                    "CHR2", "START2", "END2", "ID2"
                )
            ) %>%
                select(ID1) %>%
                mutate(overlap = 1) %>%
                distinct(ID1, .keep_all = TRUE) %>%
                rename(
                    ID = ID1,
                    !!source := overlap
                )
            source_df <<- rbind(source_df, overlap_sv)
        }
    }
    overlap_stat <<- left_join(
        x = overlap_stat,
        y = source_df,
        by = c("ID")
    )
}

overlap_stat_position <- overlap_stat %>%
    mutate_if(is.numeric, ~nafill(., fill = 0)) %>%
    left_join(
        y = weibin_all,
        by = c("ID")
    ) %>%
    mutate(length = END - START) %>%
    relocate(CHR, START, END, length, .after = ID) %>%
    mutate(
        SRS = ifelse(gnomAD + HGDP > 0, 1, 0),
        LRS = ifelse(HGSVC + LRS15 + Iceland + Han945 > 0, 1, 0),
        reported = ifelse(gnomAD + HGDP + HGSVC + Iceland + LRS15 + Han945 > 0, 1, 0)
    )

fwrite(
    x = overlap_stat_position,
    file = "results/SV_compare/stat/SV_compare_stat.txt",
    sep = "\t",
    quote = FALSE
)

plotdata_amount <- overlap_stat_position %>%
    select(ID, sv_type, gnomAD, HGDP, HGSVC, Iceland, LRS15, Han945) %>%
    pivot_longer(
        cols = -c(ID, sv_type),
        names_to = "source",
        values_to = "overlap"
    ) %>%
    mutate(
        source = factor(source, levels = data_sources),
        sv_type = factor(sv_type, levels = c("INS", "DEL", "INV", "DUP"))
    ) %>%
    mutate(type = ifelse(overlap == 0, "novel", "share")) %>%
    group_by(source, sv_type, type) %>%
    summarise(n = n()) %>%
    mutate(percentage = round(n / sum(n) * 100, 1))

plot_long <- ggplot(
    data = plotdata_amount,
    mapping = aes(
        x = source,
        y = percentage,
        fill = type
    )
) +
    geom_bar(
        position = position_stack(),
        width = 0.6,
        stat = "identity"
    ) +
    scale_y_continuous(
        limits = c(0, 100),
        expand = c(0, 0)
    ) +
    coord_flip() +
    scale_fill_manual(
        values = c(
            "share" = "#82B0D2",
            "novel" = "#E7F3FC"
        )
    ) +
    facet_wrap(
        ~sv_type,
        ncol = 1,
        strip.position = "right"
    ) +
    labs(
        x = "Data source",
        y = "Percentage (%)",
        fill = NULL
    ) +
    wdy_theme(
        base_size = 10,
        axis.line = element_blank(),
        legend.position = "top",
        panel.spacing.y = unit(5, "pt"),
        strip.background = element_rect(fill = "#ebf1f5", linewidth = NA)
    )

plot_wide <- ggplot(
    data = plotdata_amount,
    mapping = aes(
        x = sv_type,
        y = percentage,
        fill = type
    )
) +
    geom_bar(
        position = position_stack(),
        width = 0.6,
        stat = "identity"
    ) +
    scale_y_continuous(
        limits = c(0, 100),
        expand = c(0, 0)
    ) +
    scale_fill_manual(
        values = c(
            "share" = "#82B0D2",
            "novel" = "#E7F3FC"
        )
    ) +
    facet_wrap(
        ~source,
        nrow = 1,
        strip.position = "top"
    ) +
    labs(
        x = "SV type",
        y = "Percentage (%)",
        fill = NULL
    ) +
    wdy_theme(
        base_size = 10,
        axis.line = element_blank(),
        legend.position = "top",
        panel.spacing.x = unit(5, "pt"),
        strip.background = element_rect(fill = "#ebf1f5", linewidth = NA)
    )

plot_all <- ggplot(
    data = plotdata_amount %>%
        group_by(source, type) %>%
        summarise(n = sum(n)) %>%
        mutate(percentage = round(n / sum(n) * 100, 1)),
    mapping = aes(
        x = source,
        y = percentage,
        fill = type
    )
) +
    geom_bar(
        position = position_stack(),
        width = 0.6,
        stat = "identity"
    ) +
    scale_y_continuous(
        limits = c(0, 100),
        expand = c(0, 5)
    ) +
    scale_fill_manual(
        values = c(
            "share" = "#82B0D2",
            "novel" = "#E7F3FC"
        )
    ) +
    labs(
        x = "Data source",
        y = "Percentage (%)",
        fill = NULL
    ) +
    coord_flip() +
    wdy_theme(
        base_size = 12,
        axis.line = element_blank(),
        legend.position = "top",
        panel.spacing.x = unit(5, "pt"),
        strip.background = element_rect(fill = "#ebf1f5", linewidth = NA)
    )

pdf("results/SV_compare/stat/SV_compare_long.pdf", width = 2.6, height = 6)
plot(plot_long)
dev.off()

pdf("results/SV_compare/stat/SV_compare_wide.pdf", width = 8, height = 2.6)
plot(plot_wide)
dev.off()

pdf("results/SV_compare/stat/SV_compare_allsv.pdf", width = 3, height = 3)
plot(plot_all)
dev.off()

plotdata_length <- overlap_stat_position %>%
    select(ID, sv_type, length, SRS, LRS, reported) %>%
    mutate(
        SRS = ifelse(SRS == 1, "SRS", NA),
        LRS = ifelse(LRS == 1, "LRS", NA),
        reported = ifelse(reported == 0, "Novel", NA)
    ) %>%
    pivot_longer(
        cols = -c(ID, sv_type, length),
        names_to = "col",
        values_to = "platform"
    ) %>%
    select(-col) %>%
    filter(!is.na(platform)) %>%
    mutate(
        platform = factor(platform, levels = c("SRS", "LRS", "Novel")),
        sv_type = factor(sv_type, levels = c("INS", "DEL", "INV", "DUP")),
        length_str = cut(
            length,
            breaks = c(50, 500, 1000, 1500, 2000, Inf),
            labels = c(
                "50-500",
                "500-1k",
                "1k-1.5k",
                "1.5k-2k",
                "2k+"
            ),
            include.lowest = TRUE
        )
    ) %>%
    group_by(platform, length_str) %>%
    summarise(n = n()) %>%
    mutate(percentage = round(n / sum(n) * 100, 1))

plot_length <- ggplot(
    data = plotdata_length,
    mapping = aes(
        x = length_str,
        y = percentage,
        fill = platform
    )
) +
    geom_bar(
        stat = "identity",
        position = position_dodge(),
        width = 0.6
    ) +
    scale_y_continuous(
        expand = c(0, 0),
        limits = c(0, 100)
    ) +
    scale_fill_manual(values = custom_colors) +
    labs(
        x = "SV length",
        y = "Percentage (%)",
        fill = "Type"
    ) +
    wdy_theme(
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.99, 0.99),
        legend.justification = c(1, 1)
    )

pdf("results/SV_compare/stat/length_overlap.pdf", width = 3, height = 2.5)
plot(plot_length)
dev.off()
