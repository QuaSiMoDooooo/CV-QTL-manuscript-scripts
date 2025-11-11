#! Rscript
# -------------
# FileName     : matrixEQTL
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-10-22 20:57
# Last Modified: 2025-01-20 22:06
# Modified By  : EastsunW
# -------------
# Description  : 进行QTL分析, 注意不会区分顺式和反式
# -------------

# 加载必要的包
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(MatrixEQTL))


get_options <- function() {
    usage <- "usage: Rscript %prog [options] <arguments>"
    description <- "自动执行QTL分析，需要匹配的基因型和表型数据，变异和位置也需要匹配"
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
            c("--variant_genotype"),
            type = "character",
            dest = "variant_genotype",
            metavar = "<filepath>",
            action = "store",
            help = "变异的基因型矩阵文件路径, 每行一个变异, 每列一个样本，第一列是变异ID，第一行是列名"
        ),
        make_option(
            c("--phenotype_quantity"),
            type    = "character",
            dest    = "phenotype_quantity",
            metavar = "<filepath>",
            action  = "store",
            help    = "表型的定量矩阵文件, 每行一个表型marker, 每列一个样本，第一列是markerID，第一行是列名"
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
            c("--covariate"),
            type = "character",
            dest = "covariate",
            metavar = "<filepath>",
            action = "store",
            help = "协变量文件的文件, 每行一个协变量, 每列一个样本，第一列是协变量的ID，第一行是列名"
        ),
        make_option(
            c("--pvalue"),
            type    = "double",
            dest    = "pvalue",
            metavar = "<path>",
            default = 1,
            action  = "store",
            help    = "QTL的pvalue阈值, 默认为1不筛选"
        ),
        make_option(
            c("-w", "--window"),
            type    = "numeric",
            dest    = "window",
            metavar = "<path>",
            default = 1e6,
            action  = "store",
            help    = "顺式QTL的窗口长度"
        ),
        make_option(
            c("--output"),
            type    = "character",
            dest    = "output",
            metavar = "<path>",
            action  = "store",
            help    = "原始QTL结果的输出路径前缀"
        ),
        make_option(
            c("--threads"),
            type    = "numeric",
            dest    = "threads",
            metavar = "<num>",
            action  = "store",
            help    = "原始QTL结果的输出路径前缀"
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

# -------------------- 读取基因型 -------------------- #
suppressMessages({
    variant_genotype <- SlicedData$new()
    variant_genotype$fileDelimiter <- "\t"
    variant_genotype$fileOmitCharacters <- "NA"
    variant_genotype$fileSkipRows <- 1
    variant_genotype$fileSkipColumns <- 1
    variant_genotype$fileSliceSize <- 5000
    variant_genotype$LoadFile(options[["variant_genotype"]])
})
variant_position <- as.data.frame(fread(
    options[["variant_position"]],
    sep = "\t",
    fill = TRUE
))

# -------------------- 读取表型 -------------------- #
## 读取原始表型数据
phenotype_quantity_raw <- fread(
    options[["phenotype_quantity"]],
    sep = "\t",
    fill = TRUE
)
phenotype_splited <- split(
    x = phenotype_quantity_raw,
    f = seq(0, nrow(phenotype_quantity_raw) - 1) %/% (nrow(phenotype_quantity_raw) %/% options[["threads"]] + 1) + 1
)
phenotype_position <- as.data.frame(fread(
    options[["phenotype_position"]],
    sep = "\t",
    fill = TRUE
))

# -------------------- 读取协变量 -------------------- #
suppressMessages({
    covariates <- SlicedData$new()
    covariates$fileDelimiter <- "\t"
    covariates$fileOmitCharacters <- "NA"
    covariates$fileSkipRows <- 1
    covariates$fileSkipColumns <- 1
    covariates$LoadFile(options[["covariate"]])
})

phenotype_splited_sliced <- parallel::mclapply(
    X = phenotype_splited,
    FUN = function(slice) {
        slice_matrix <- slice %>%
            column_to_rownames("ID") %>%
            as.matrix()
        # 将矩阵转换为SlicedData对象
        phenotype_quantity <- SlicedData$new()
        phenotype_quantity$CreateFromMatrix(slice_matrix)
        phenotype_quantity$ResliceCombined(5000)
        return(phenotype_quantity)
    },
    mc.cores = options[["threads"]]
)

# -------------------- 进行QTL分析 -------------------- #
temp <- parallel::mclapply(
    X = names(phenotype_splited_sliced),
    FUN = function(slice_name) {
        suppressMessages({
            qtl_obj <- Matrix_eQTL_main(
                useModel              = modelLINEAR,
                snps                  = variant_genotype,
                gene                  = phenotype_splited_sliced[[slice_name]],
                cvrt                  = covariates,
                snpspos               = variant_position[, c(1:3)],
                genepos               = phenotype_position[, c(1:4)],
                cisDist               = options[["window"]],
                verbose               = FALSE,
                min.pv.by.genesnp     = FALSE,
                noFDRsaveMemory       = TRUE,
                errorCovariance       = numeric(),
                pvOutputThreshold.cis = 1,  # cis全部输出
                pvOutputThreshold     = 0,  # trans全部不输出
                output_file_name.cis  = paste0(options[["output"]], "_", slice_name)  # cis
            )
        })
    },
    mc.cores = options[["threads"]]
)
