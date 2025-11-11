# -------------
# FileName     : format_QTL_result
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-11-02 21:24
# Last Modified: 2024-11-02 21:25
# Modified By  : EastsunW
# -------------
# Description  : 给合并后的QTL结果加上一些统计量，并计算FDR值
# -------------

# 加载必要的包
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

qtl_input <- fread(snakemake@input[[1]], sep = "\t")



# -------------------- 处理近端QTL的结果-------------------- #
if (QTL[["cis"]][["neqtls"]] > 0) {
    # 计算效应值标准差和皮尔森相关系数
    cisQTL_formated <- QTL[["cis"]][["eqtls"]] %>%
        setnames(new = c(
            "variant",
            "marker",
            "t_statistic",
            "pvalue",
            "fdr",
            "beta"
        )) %>%
        mutate(
            df = QTL[["param"]][["dfFull"]],
            se = beta / t_statistic,
            r  = t_statistic / sqrt(df + t_statistic^2)
        ) %>%
        # 加上变异的位置
        left_join(
            y = variant_position %>%
                rename(
                    variant = ID,
                    variant_chr = CHR,
                    variant_start = START
                ) %>%
                select(
                    variant,
                    variant_chr,
                    variant_start
                ) %>%
                unite(
                    variant_chr, variant_start,
                    col = "variant_position",
                    sep = ":"
                ),
            by = "variant"
        ) %>%
        relocate("variant_position", .after = "variant") %>%
        # 加上表型的位置
        left_join(
            y = phenotype_position %>%
                rename(
                    marker = ID,
                    marker_chr = CHR,
                    marker_start = START,
                    marker_end = END,
                    marker_strand = STRAND
                ) %>%
                unite(
                    marker_start, marker_end,
                    col = "marker_range",
                    sep = "-"
                ) %>%
                unite(
                    marker_chr, marker_range, marker_strand,
                    col = "marker_position",
                    sep = ":"
                ),
            by = "marker"
        ) %>%
        relocate("marker_position", .after = "marker")

    fwrite(
        x    = cisQTL_formated,
        file = paste0(options[["output"]], ".cis.txt"),
        sep  = "\t"
    )
}

# -------------------- 处理远端QTL的结果-------------------- #
if (QTL[["trans"]][["neqtls"]] > 0) {
    # 计算效应值标准差和皮尔森相关系数
    transQTL_formated <- QTL[["trans"]][["eqtls"]] %>%
        setnames(new = c(
            "variant",
            "marker",
            "t_statistic",
            "pvalue",
            "fdr",
            "beta"
        )) %>%
        mutate(
            df = QTL[["param"]][["dfFull"]],
            se = beta / t_statistic,
            r  = t_statistic / sqrt(df + t_statistic^2)
        ) %>%
        # 加上变异的位置
        left_join(
            y = variant_position %>%
                rename(
                    variant = ID,
                    variant_chr = CHR,
                    variant_start = START
                ) %>%
                select(
                    variant,
                    variant_chr,
                    variant_start
                ) %>%
                unite(
                    variant_chr, variant_start,
                    col = "variant_position",
                    sep = ":"
                ),
            by = "variant"
        ) %>%
        relocate("variant_position", .after = "variant") %>%
        # 加上表型的位置
        left_join(
            y = phenotype_position %>%
                rename(
                    marker = ID,
                    marker_chr = CHR,
                    marker_start = START,
                    marker_end = END,
                    marker_strand = STRAND
                ) %>%
                unite(
                    marker_start, marker_end,
                    col = "marker_range",
                    sep = "-"
                ) %>%
                unite(
                    marker_chr, marker_range, marker_strand,
                    col = "marker_position",
                    sep = ":"
                ),
            by = "marker"
        ) %>%
        relocate("marker_position", .after = "marker")

    fwrite(
        x    = transQTL_formated,
        file = paste0(options[["output"]], ".trans.txt"),
        sep  = "\t"
    )
}
