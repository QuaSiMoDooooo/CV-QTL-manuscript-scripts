#! /usr/bin/Rscript
# -------------
# FileName     : format_biochemistry
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-09-04 15:49
# Last Modified: 2024-09-20 19:54
# Modified By  : EastsunW
# -------------
# Description  : 将数值型和非数值型的变量以及性别特有的指标分开，先对整个数值型的变量进行缺失值填充，然后对数值型变量进行z-score标准化，接着从中分离出男性和女性特有的指标。
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(mice))

non_numeric_indicators <- c(
    "FOBT",
    "Anti-ANA antibody",
    "Nuclear fluorescence model",
    "Cytoplasm fluorescence model",
    "Cell mitosis fluorescence model",
    "Urine Glyphosate",
    "Urine mycotoxins",
    "Lung CT",
    "Head CT",
    "Pancreas ultrasound",
    "Breast ultrasound",
    "Urinary tract ultrasound",
    "Thyroid ultrasound",
    "Cardiac ultrasound",
    "Carotid artery ultrasound",
    "Lower limb artery ultrasound",
    "Gynecological ultrasound",
    "Prostate ultrasound"
)

male_only_indicators <- c(
    "fPSA",
    "tPSA",
    "Testosterone",
    "Free testosterone"
)
female_only_indicators <- c(
    "β-hCG",
    "HE4",
    "CA15-3",
    "FSH",
    "AMH",
    "Progesterone"
)

data_raw <- fread(
    snakemake@input[[1]],
    sep = "\t",
    na.strings = ""
) %>%
    select_at(vars(c("Indicator", starts_with("HN-"))))

data_nonnum <- data_raw %>%
    filter(Indicator %in% non_numeric_indicators)

# 将数值型变量转换为数值型
data_num <- data_raw %>%
    filter(!Indicator %in% non_numeric_indicators) %>%
    mutate_at(vars(starts_with("HN-")), ~ gsub("^(<|>)", "", .)) %>%
    mutate_at(vars(starts_with("HN-")), as.numeric)

data_num_common <- data_num %>%
    filter(!Indicator %in% male_only_indicators) %>%
    filter(!Indicator %in% female_only_indicators)
data_num_male <- data_num %>%
    filter(Indicator %in% male_only_indicators) %>%
    select_at(vars(c(1:7, starts_with("HN-1"))))
data_num_female <- data_num %>%
    filter(Indicator %in% female_only_indicators) %>%
    select_at(vars(c(1:7, starts_with("HN-2"))))

fwrite(
    data_nonnum,
    snakemake@output[["non_numeric"]],
    sep = "\t"
)
fwrite(
    data_num_common,
    snakemake@output[["numeric_common"]],
    sep = "\t"
)
fwrite(
    data_num_male,
    snakemake@output[["numeric_male"]],
    sep = "\t"
)
fwrite(
    data_num_female,
    snakemake@output[["numeric_female"]],
    sep = "\t"
)
