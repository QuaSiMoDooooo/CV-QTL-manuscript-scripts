#! Rscript
# -------------
# FileName     : format_QTL_results
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2025-01-20 16:00
# Last Modified: 2025-01-20 16:52
# Modified By  : EastsunW
# -------------
# Description  : 将QTL的结果进行处理
# -------------

# 加载必要的包
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))

get_options <- function() {
    usage <- "usage: Rscript %prog [options] <arguments>"
    description <- "给QTL结果加上表型和变异的位置和统计量"
    option_list <- list(
        make_option(
            c("--variant_position"),
            type = "character",
            dest = "variant_position",
            metavar = "<filepath>",
            action = "store",
            help = "变异的位置文件路径, 前三列必须是ID、染色体、位置（start）"
        ),
        make_option(
            c("--phenotype_position"),
            type    = "character",
            dest    = "phenotype_position",
            metavar = "<filepath>",
            action  = "store",
            help    = "表型定量的位置文件, 包含名字, 染色体, 开始位置, 结束位置四列"
        ),
        make_option(
            c("--qtl"),
            type = "character",
            dest = "qtl",
            metavar = "<filepath>",
            action = "store",
            help = "matrixEQTL的原始结果"
        ),
        make_option(
            c("--window"),
            type    = "numeric",
            dest    = "window",
            metavar = "<path>",
            default = 1e6,
            action  = "store",
            help    = "顺式QTL的窗口长度"
        ),
        make_option(
            c("--o_cis"),
            type    = "character",
            dest    = "o_cis",
            metavar = "<path>",
            action  = "store",
            help    = "cis输出文件的路径"
        ),
        make_option(
            c("--o_trans"),
            type    = "character",
            dest    = "o_trans",
            metavar = "<path>",
            action  = "store",
            help    = "trans输出文件的路径"
        )
    )
    opt_obj <- OptionParser(
        usage           = usage,
        description     = description,
        option_list     = option_list,
        add_help_option = TRUE
    )
    opts <- optparse::parse_args(opt_obj)
    if (length(commandArgs(TRUE)) == 0) {
        optparse::print_help(opt_obj)
        if (interactive()) {
            stop("help requested")
        } else {
            quit(status = 0)
        }
    } else {
        return(opts)
    }
}

get_distance <- function(vl, vr, pl, pr) {
    # vl:variant_left, vr:variant_right
    # pl:pheontype_left, pr:phenotype_right
    o_len <- abs(max(vl, vr, pl, pr) - min(vl, vr, pl, pr)) -
        abs(max(vl, vr) - max(pl, pr)) -
        abs(min(vl, vr) - min(pl, pr))
    if (o_len >= 0) {
        return(0)
    } else {
        if (vr < pl) {
            return(vr - pl)
        } else if (vl > pr) {
            return(vl - pr)
        }
    }
}


chr_len <- c(
    "chr1"  = 248956422, "chr2"  = 242193529,
    "chr3"  = 198295559, "chr4"  = 190214555,
    "chr5"  = 181538259, "chr6"  = 170805979,
    "chr7"  = 159345973, "chr8"  = 145138636,
    "chr9"  = 138394717, "chr10" = 133797422,
    "chr11" = 135086622, "chr12" = 133275309,
    "chr13" = 114364328, "chr14" = 107043718,
    "chr15" = 101991189, "chr16" = 90338345,
    "chr17" = 83257441,  "chr18" = 80373285,
    "chr19" = 58617616,  "chr20" = 64444167,
    "chr21" = 46709983,  "chr22" = 50818468,
    "chrX"  = 156040895, "chrY"  = 57227415
)

options <- get_options()
if (any(unlist(options) == "")) {
    optparse::print_help(opt_obj)
    stop(
        paste0(
            "These argument are not specified:\n\t --",
            paste0(
                names(unlist(options)[unlist(options) == ""]),
                collapse = "\n\t --"
            )
        )
    )
}

variant_position <- fread(
    options[["variant_position"]],
    sep = "\t",
    fill = TRUE
)
phenotype_position <- fread(
    options[["phenotype_position"]],
    sep = "\t",
    fill = TRUE
)

qtl_all_processed <- fread(
    options[["qtl"]],
    sep = "\t",
    fill = TRUE
) %>%
    rename(
        Variant = SNP,
        Phenotype = gene,
        "T_stat" = "t-stat",
        "P" = "p-value",
        "Beta" = "beta"
    ) %>%
    mutate(
        Df = 124,
        Se = Beta / T_stat,
        R_pearson = T_stat / sqrt(Df + T_stat^2)
    ) %>%
    left_join(
        y = variant_position,
        by = c("Variant" = "ID")
    ) %>%
    rename(
        variant_chr = CHR,
        variant_start = START,
        variant_end = END,
        Ref = REF,
        Alt = ALT
    ) %>%
    left_join(
        y = phenotype_position,
        by = c("Phenotype" = "ID")
    ) %>%
    rename(
        phenotype_chr = CHR,
        phenotype_start = START,
        phenotype_end = END,
        phenotype_strand = STRAND,
        Gene = GENE,
    ) %>%
    mutate(distance = apply(., 1, function(row) {
        if (row["variant_chr"] != row["phenotype_chr"]) {
            return(Inf)
        } else {
            return(get_distance(
                as.numeric(row["variant_start"]),
                as.numeric(row["variant_end"]),
                as.numeric(row["phenotype_start"]),
                as.numeric(row["phenotype_end"])
            ))
        }
    })) %>%
    mutate(qtl_type = ifelse(
        abs(distance) <= options[["window"]],
        "cis",
        "trans"
    )) %>%
    unite(
        variant_start, variant_end,
        col = "variant_pos",
        sep = "-"
    ) %>%
    unite(
        variant_chr, variant_pos,
        col = "Variant_position",
        sep = ":"
    ) %>%
    relocate(Variant_position, Ref, Alt, .after = Variant) %>%
    unite(
        phenotype_start, phenotype_end,
        col = "phenotype_pos",
        sep = "-"
    ) %>%
    unite(
        phenotype_chr, phenotype_pos, phenotype_strand,
        col = "Phenotype_position",
        sep = ":"
    ) %>%
    relocate(Phenotype_position, Gene, .after = Phenotype) %>%
    group_by(qtl_type) %>%
    mutate(FDR = p.adjust(P, method = "fdr")) %>%
    relocate(FDR, .after = P) %>%
    arrange(FDR, P)

fwrite(
    qtl_all_processed %>% filter(qtl_type == "cis"),
    options[["o_cis"]],
    sep = "\t",
    quote = FALSE
)
fwrite(
    qtl_all_processed %>% filter(qtl_type == "trans"),
    options[["o_trans"]],
    sep = "\t",
    quote = FALSE
)
