#! /usr/bin/Rscript
# -------------
# FileName     : impute_data
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-09-19 09:25
# Last Modified: 2024-09-19 16:22
# Modified By  : EastsunW
# -------------
# Description  : 输入一个数据框，插补其中的缺失值
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(mice))

library("optparse")
method_help <- "Methods for imputation:
                    Method                 Data type    Description
                    ===========================================================================
                    pmm                    any          Predictive mean matching
                    midastouch             any          Weighted predictive mean matching
                    sample                 any          Random sample from observed values
                    cart                   any          Classification and regression trees
                    rf                     any          Random forest imputations
                    mean                   numeric      Unconditional mean imputation
                    norm                   numeric      Bayesian linear regression
                    norm.nob               numeric      Linear regression ignoring model error
                    norm.boot              numeric      Linear regression using bootstrap
                    norm.predict           numeric      Linear regression, predicted values
                    lasso.norm             numeric      Lasso linear regression
                    lasso.select.norm      numeric      Lasso select + linear regression
                    quadratic              numeric      Imputation of quadratic terms
                    ri                     numeric      Random indicator for nonignorable data
                    logreg                 binary       Logistic regression
                    logreg.boot            binary       Logistic regression with bootstrap
                    lasso.logreg           binary       Lasso logistic regression
                    lasso.select.logreg    binary       Lasso select + logistic regression
                    polr                   ordered      Proportional odds model
                    polyreg                unordered    Polytomous logistic regression
                    lda                    unordered    Linear discriminant analysis
                    2l.norm                numeric      Level-1 normal heteroscedastic
                    2l.lmer                numeric      Level-1 normal homoscedastic, lmer
                    2l.pan                 numeric      Level-1 normal homoscedastic, pan
                    2l.bin                 binary       Level-1 logistic, glmer
                    2lonly.mean            numeric      Level-2 class mean
                    2lonly.norm            numeric      Level-2 class normal
                    2lonly.pmm             any          Level-2 class predictive mean matching"
get_options <- function() {
    usage       <- "usage: Rscript %prog [options] <arguments>"
    description <- "插补你输入的数据框中的缺失值，输出一个插补后的新数据框，每行一个样本，每列一个变量，第一列是样本名，第一行是变量名"
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
            default = "imputed.txt",
            action  = "store",
            help    = "Output file path"
        ),
        make_option(
            c("-m", "--method"),
            type    = "character",
            dest    = "method",
            metavar = "<method>",
            default = "pmm",
            action  = "store",
            help    = method_help
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

raw_data <- fread(
    options[["input"]],
    sep = "\t",
    na.strings = c("", "NA")
) %>%
    column_to_rownames(colnames(.)[1])

col_name <- colnames(raw_data)
colnames(raw_data) <- paste0("V", seq_len(ncol(raw_data)))

print("Imputing data...")
imputed_result <- mice(
    raw_data,
    m = 5,
    maxit = 50,
    method = options[["method"]],
    seed = 2024,
    printFlag = FALSE
)

imputed_df <- complete(imputed_result)
colnames(imputed_df) <- col_name
imputed_df <- imputed_df %>%
    rownames_to_column("Sample")

fwrite(
    imputed_df,
    options[["output"]],
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
