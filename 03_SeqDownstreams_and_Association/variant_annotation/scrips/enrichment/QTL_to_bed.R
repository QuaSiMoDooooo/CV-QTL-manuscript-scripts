#! Rscript
# -------------
# FileName     : QTL_to_bed
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-12-26 15:54
# Last Modified: 2025-03-05 20:46
# Modified By  : EastsunW
# -------------
# Description  : 将QTL结果中的变异转换成bed格式，用于基因组区域富集分析
# -------------

# 加载必要的包
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

library("optparse")

get_options <- function() {
    usage <- "usage: Rscript %prog [options] <arguments>"
    description <- "Convert QTL results to bed format"
    option_list <- list(
        make_option(
            c("-q", "--qtl"),
            type    = "character",
            dest    = "qtl",
            metavar = "<qtltype>",
            default = FALSE,
            action  = "store",
            help    = "QTL type eg. eQTL, sQTL, apaQTL, meQTL"
        ),
        make_option(
            c("-v", "--variant"),
            type    = "character",
            dest    = "variant",
            metavar = "<path>",
            default = FALSE,
            action  = "store",
            help    = "Variant type eg. SNP, InDel, MNV, SV"
        )
    )
    opt_obj <- OptionParser(
        usage           = usage,
        description     = description,
        add_help_option = TRUE,
        option_list     = option_list
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


# 读取QTL结果
base_dir <- "/home/wangdy/Projects/Weibin/Downstreams/variant_annotation/data/QTL_results"
input_file <- file.path(
    base_dir,
    paste0(options$variant, "_", options$qtl)
)
qtl_data <- list(
    cis = fread(file.path(input_file, "cis.filtered.txt.gz")),
    trans = fread(file.path(input_file, "trans.filtered.txt.gz"))
)
bed_files <- lapply(
    qtl_data,
    function(x) {
        x %>%
            select(Variant, Variant_position) %>%
            rename(ID = Variant) %>%
            separate(Variant_position, into = c("chr", "start", "end"), sep = ":|-") %>%
            relocate(ID, .after = end) %>%
            mutate(
                start = as.numeric(start) - 1,
                end = as.numeric(end)
            )
    }
)
merged_bed <- bind_rows(bed_files)
fwrite(
    merged_bed,
    file.path("data/Variants/QTL", paste0(options$variant, "_", options$qtl, ".bed")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
)
