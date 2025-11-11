#! /home/wangdy/miniforge3/envs/peer/bin/Rscript
# -------------
# FileName     : prepare_covariants
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-09-12 14:12
# Last Modified: 2024-10-14 15:33
# Modified By  : EastsunW
# -------------
# Description  : 准备表型的协变量-PEER，输入的数据每行一个指标，每列一个样本
# -------------
set.seed(2024)
library("optparse")
suppressMessages(library("peer"))
suppressPackageStartupMessages(library(data.table))

get_options <- function() {
    usage       <- "usage: Rscript %prog [options] <arguments>"
    description <- "接受一个数据框，生成指定数量的peer因子，每行一个指标（如基因），每列一个样本，第一行和第一列分别是行名和列名"
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
            c("-o", "--output"),
            type    = "character",
            dest    = "output",
            metavar = "<path>",
            default = FALSE,
            action  = "store",
            help    = "Output file path"
        ),
        make_option(
            c("-n", "--n_peer"),
            type    = "numeric",
            dest    = "n_peer",
            metavar = "<integer>",
            default = FALSE,
            action  = "store",
            help    = "Number of peer factors to generate"
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

# 根据样本数判断PEER的数量
get_peer_num <- function(sample) {
    if (sample > 400) {
        return(100)
    } else if (round(0.25 * sample) + 8 < sample - 3) {
        return(round(0.25 * sample))
    } else if (sample - 11 <= 0) {
        return(0)
    } else {
        return(sample - 3 - 8)
    }
}

# 读取数据
data_raw <- fread(
    options[["input"]],
    sep = "\t",
    na.strings = c("NA", "NaN", "NULL", "")
)
samples <- colnames(data_raw)[-1]
indicators <- data_raw[[1]]
rownames(data_raw) <- indicators
# 进行peer分析时，需要列为指标，行为样本，数据框需要是matrix格式，而非数据框
data_matrix <- t(data_raw[, ..samples])
colnames(data_matrix) <- indicators
if (options[["n_peer"]]) {
    peer_num <- options[["n_peer"]]
} else {
    peer_num <- get_peer_num(length(samples))
}

model <- PEER()
PEER_setPhenoMean(model, data_matrix)
PEER_setNk(model, peer_num)
PEER_setNmax_iterations(model, 100)

print(paste("PEER因子数：", peer_num))
invisible(PEER_update(model))

# peer完成后的结果每列是一个因子，每行是一个指标，由于QTL中需要列是样本，因此在这里转置
factors <- as.data.frame(t(PEER_getX(model)))
na_row <- apply(factors, 1, function(x) {
    all(is.na(x))
    }
)
print(dim(factors))
colnames(factors) <- samples

peer_data <- cbind(
    data.frame(cov = paste0("PEER", c(1:peer_num))),
    factors
)
fwrite(
    x = peer_data,
    file = options[["output"]],
    sep = "\t"
)
