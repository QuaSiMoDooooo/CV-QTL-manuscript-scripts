#! Rscript
# -------------
# FileName     : split_file
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-11-01 22:21
# Last Modified: 2024-11-01 22:52
# Modified By  : EastsunW
# -------------
# Description  :
# -------------

# 加载必要的包
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
library("optparse")

get_options <- function() {
    usage <- "usage: Rscript %prog [options] <arguments>"
    description <- "description"
    option_list <- list(
        make_option(
            c("-i", "--input"),
            type    = "character",
            dest    = "input",
            metavar = "<path>",
            default = FALSE,
            action  = "store",
            help    = "Input file path"
        ),
        make_option(
            c("-n", "--n_split"),
            type    = "numeric",
            dest    = "n_split",
            metavar = "<int>",
            default = FALSE,
            action  = "store",
            help    = "Total file to be splited"
        ),
        make_option(
            c("--rank"),
            type    = "numeric",
            dest    = "rank",
            metavar = "<int>",
            default = FALSE,
            action  = "store",
            help    = "The n.st file of splited file"
        ),
        make_option(
            c("-o", "--output"),
            type    = "character",
            dest    = "output",
            metavar = "<path>",
            default = "%slice%.txt",
            action  = "store",
            help    = "Output file path, %slice% will be replaced by slice number"
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


get_split_index <- function(n_split, total) {
    per_part <- ceiling(total / n_split)
    part_content <- c(0, rep(per_part, n_split - 1), total - (n_split - 1) * per_part)
    result <- list()
    for (i in seq_len(n_split)) {
        result <- c(
            result,
            list(c(
                1 + (i - 1) * part_content[i],
                (i - 1) * part_content[i] + part_content[i + 1]
            ))
        )
    }
    return(result)
}


df <- fread(
    options[["input"]],
    sep = "\t"
)
split_numbers <- get_split_index(
    options[["n_split"]], nrow(df)
)
fwrite(
    x = df[split_numbers[[options[["rank"]]]][1]:split_numbers[[options[["rank"]]]][2], ],
    file = gsub("%slice%", sprintf("%02d", options[["rank"]]), options[["output"]]),
    sep = "\t",
    quote = FALSE,
    col.names = TRUE
)
