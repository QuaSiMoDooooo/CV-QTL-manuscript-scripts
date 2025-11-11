#!/usr/bin/env Rscript

# -------------
# FileName     : 01_script.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Prepare QTL data for SMR analysis by creating ESD and flist files
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))

args <- commandArgs(trailingOnly = TRUE)
qtl <- args[1]
qtl_name <- args[2]
type <- args[3]

cur_path <- "/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/07_QTL_GWAS_SMR"

rsid <- fread("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/07_QTL_GWAS_SMR/simplify_rsid_ref_bed/all_qtl_filtered.sorted.uniq.rsid.redundant.bed")
judge_vec = !duplicated(paste(rsid$V1,rsid$V2,sep = "_"))
rsid_ref = data.frame(chr=rsid$V1, pos=rsid$V2, rsid=rsid$V4)[judge_vec,]

frq <- fread("/home/wtian/project/HZAU_cohort_meth/wdy_result/SNP_plink_genotype_process_bim3/SNP.biochem_common.frq")
tmp = str_split(frq$SNP,"_",simplify = T)
frq$chr_pos = paste(tmp[,1],tmp[,2],sep = "_")
rsid_ref$chr_pos = paste(rsid_ref$chr,rsid_ref$pos,sep = "_")
frq_ref = frq[frq$chr_pos %in% rsid_ref$chr_pos,]
rsid_ref$chr_pos = NULL
tmp = select(frq_ref,chr_pos,MAF)
tmp$chr = str_split(tmp$chr_pos,"_",simplify = T)[,1]
tmp$pos = str_split(tmp$chr_pos,"_",simplify = T)[,2]
tmp$chr_pos = NULL
frq_ref = select(tmp,chr,pos,MAF)
frq_ref = frq_ref[frq_ref$MAF > 0,]

rsid_ref$chr_pos = paste(rsid_ref$chr,rsid_ref$pos,sep = "_")
frq_ref$chr_pos = paste(frq_ref$chr,frq_ref$pos,sep = "_")
rsid_ref$chr = NULL
rsid_ref$pos = NULL
frq_ref$chr = NULL
frq_ref$pos = NULL

qtl_dir <- "/home/wangdy/Projects/Weibin/Downstreams/phenotype_qtl/results"
file_path <- paste0(qtl_dir, "/SNP-", qtl_name, "/QTL_results/", type, ".filtered.txt.gz")

data <- fread(file_path, showProgress = FALSE, sep = "\t")
tmp <- str_split(data$Variant, "_", simplify = TRUE)
data$chr_pos <- paste(tmp[,1], tmp[,2], sep = "_")
data$chr <- tmp[,1]
data$pos <- tmp[,2]

data <- left_join(data, rsid_ref, by = "chr_pos")
data <- data[!is.na(data$rsid), ]

data <- left_join(data, frq_ref, by = "chr_pos")
data <- data[!is.na(data$MAF), ]

data_slt <- data %>% select(chr, pos, rsid, Phenotype, Alt, Ref, MAF, Beta, P, Se)
colnames(data_slt) <- c("Chr", "BP", "SNP", "phe", "A1", "A2", "Freq", "Beta", "p", "se")
data_slt <- data_slt %>% select(phe, Chr, SNP, BP, A1, A2, Freq, Beta, se, p) %>% arrange(phe, Chr, BP)
data_slt$Chr <- sub("chr", "", data_slt$Chr)

esd_dir = "01_esd"
if (!dir.exists(esd_dir)){
    dir.create(esd_dir)
}
type_qtl_dir = paste0(esd_dir,"/",type,"_",qtl)
if (!dir.exists(type_qtl_dir)){
    dir.create(type_qtl_dir)
}
for (i in unique(data_slt$phe)) {
    data_slt_tmp <- data_slt[data_slt$phe == i, ]
    file_name <- paste0(type_qtl_dir, "/", i, ".esd")
    write.table(data_slt_tmp[, -c("phe")], file_name, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
}

flist_dir <- "01_flist"
if (!dir.exists(flist_dir)){
    dir.create(flist_dir)
}
flist_name <- paste0(flist_dir, "/", type, "_", qtl_name, ".flist")
if (!file.exists(flist_name)) {
    if (type == "cis") {
        tmp <- data[!duplicated(data$Phenotype), ]
        tmp$ProbeBp <- as.numeric(tmp$pos) - as.numeric(tmp$distance)
        tmp$phechr <- str_split(tmp$Phenotype_position, ":", simplify = TRUE)[, 1]
        tmp$phechr <- sub("chr", "", tmp$phechr)
        flist <- data.frame(ProbeID = tmp$Phenotype, Chr = tmp$phechr, ProbeBp = tmp$ProbeBp, Gene = tmp$Gene, GeneticDistance = 0, Orientation = NA)
        esd_path_list <- paste0(cur_path, "/", esd_dir, "/", type, "_", qtl, "/", flist$ProbeID, ".esd")
        flist$PathOfEsd <- esd_path_list
        flist <- flist %>% select(Chr, ProbeID, GeneticDistance, ProbeBp, Gene, Orientation, PathOfEsd)
        if (all(is.na(flist$Gene))) flist$Gene <- NA
        write.table(flist, flist_name, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
        print(paste0("flist file generated: ", flist_name))
    } else if (type == "trans") {
        data1 = data[data$distance != "Inf", ]
        data2 = data[data$distance == "Inf", ]
        tmp1 = data1[!duplicated(data1$Phenotype), ]
        tmp1$ProbeBp <- as.numeric(tmp1$pos) - as.numeric(tmp1$distance)
        tmp1$phechr <- str_split(tmp1$Phenotype_position, ":", simplify = TRUE)[, 1]
        tmp1$phechr <- sub("chr", "", tmp1$phechr)
        data2 = data2[!data2$Phenotype %in% tmp1$Phenotype, ]
        tmp2 = data2[!duplicated(data2$Phenotype), ]
        tmp2$ProbeBp <- 999
        tmp2$phechr <- str_split(tmp2$Phenotype_position, ":", simplify = TRUE)[, 1]
        tmp2$phechr <- sub("chr", "", tmp2$phechr)
        tmp = rbind(tmp1, tmp2)
        flist <- data.frame(ProbeID = tmp$Phenotype, Chr = tmp$phechr, ProbeBp = tmp$ProbeBp, Gene = tmp$Gene, GeneticDistance = 0, Orientation = NA)
        esd_path_list <- paste0(cur_path, "/", esd_dir, "/", type, "_", qtl, "/", flist$ProbeID, ".esd")
        flist$PathOfEsd <- esd_path_list
        flist <- flist %>% select(Chr, ProbeID, GeneticDistance, ProbeBp, Gene, Orientation, PathOfEsd)
        if (all(is.na(flist$Gene))) flist$Gene <- NA
        write.table(flist, flist_name, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
        print(paste0("flist file generated: ", flist_name))
    }
}