#! Rscript
# -------------
# FileName     : gnomAD
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : compare 2-MNVs and InDels in gnomAD with CV-datasets
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(GenomicRanges))
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
options(stringsAsFactors = FALSE)
# options(scipen = 999)
mc <- getOption('mc.cores', 36)
getwd()

# MNV
wb = fread("GRCh38_241031.mnv.withID.txt")
tmp = str_split(wb$pos,",",simplify = T)
wb_bed = data.frame(chr = wb$`#chr`,
                    snv1 = tmp[,1],
                    snv2 = tmp[,2])
table(wb_bed$snv2 > wb_bed$snv1)
gnomad = fread("gnomad.mnv2.DisAll.bed")
table(gnomad$snv2 > gnomad$snv1)
wb_bed$mnv_site = paste(wb_bed$chr,wb_bed$snv1,wb_bed$snv2,sep = "_")
gnomad$mnv_site = paste(gnomad$chr,gnomad$snv1,gnomad$snv2,sep = "_")
table(wb_bed$mnv_site %in% gnomad$mnv_site)

# InDel
wb = fread("InDel.all.bed")
tmp = str_split(wb$V4,"_",simplify = T)
table(length(tmp[,3])==length(tmp[,4]))
table(unlist(lapply(1:length(tmp[,3]),function(i){
  return(nchar(tmp[i,3])==nchar(tmp[i,4]))
})))
wb$type = ifelse(nchar(tmp[,3])>nchar(tmp[,4]),"DEL","INS")
wb$V5 = NULL
colnames(wb) = c("chr","start","end","id","type")
head(wb)
wb_df = data.frame(
  chr = wb$chr,
  start = wb$start,
  end = wb$end,
  type = wb$type
)
exome = fread("gnomad.exomes.v4.1.sites.indel.info")
genome = fread("gnomad.genomes.v4.1.sites.indel.info")
gnomad = rbind(exome,genome)
gnomad$v3_nchar = nchar(gnomad$V3)
gnomad$v4_nchar = nchar(gnomad$V4)
gnomad$V3 = NULL
gnomad$V4 = NULL
gnomad$end = gnomad$V2 + max(gnomad$v3_nchar,gnomad$v4_nchar) - 1
gnomad$type = ifelse(gnomad$v3_nchar>gnomad$v4_nchar,"DEL","INS")
gnomad_df = data.frame(
  chr = gnomad$V1,
  start = gnomad$V2,
  end = gnomad$end,
  type = gnomad$type
)
head(gnomad_df)
gnomad_df = gnomad_df %>% distinct()
calc_overlap_ratio <- function(wb_df, gnomad_df, sv_type, mc.cores = 64) {
  wb_sub <- wb_df %>% filter(type == sv_type)
  gnomad_sub <- gnomad_df %>% filter(type == sv_type)
  if (nrow(wb_sub) == 0 | nrow(gnomad_sub) == 0) return(NA)
  wb_gr <- GRanges(seqnames = wb_sub$chr,
                   ranges = IRanges(start = wb_sub$start, end = wb_sub$end))
  gnomad_gr <- GRanges(seqnames = gnomad_sub$chr,
                       ranges = IRanges(start = gnomad_sub$start, end = gnomad_sub$end))
  hits <- findOverlaps(wb_gr, gnomad_gr)
  if (length(hits) == 0) return(0)
  overlap_valid_vec <- mclapply(seq_along(hits), function(i) {
    wb_idx <- queryHits(hits)[i]
    gn_idx <- subjectHits(hits)[i]
    ov_start <- max(start(wb_gr[wb_idx]), start(gnomad_gr[gn_idx]))
    ov_end <- min(end(wb_gr[wb_idx]), end(gnomad_gr[gn_idx]))
    ov_len <- max(0, ov_end - ov_start + 1)
    wb_len <- width(wb_gr[wb_idx])
    gn_len <- width(gnomad_gr[gn_idx])
    if (ov_len >= 0.5 * wb_len || ov_len >= 0.5 * gn_len) {
      return(wb_idx)
    } else {
      return(NA_integer_)
    }
  }, mc.cores = mc.cores)
  overlap_valid_idx <- unique(na.omit(unlist(overlap_valid_vec)))
  ratio <- length(overlap_valid_idx) / length(wb_gr)
  return(ratio)
}
ratio_INS <- calc_overlap_ratio(wb_df, gnomad_df, "INS", mc.cores = 64)
ratio_DEL <- calc_overlap_ratio(wb_df, gnomad_df, "DEL", mc.cores = 64)
cat("INS overlap ratio:", ratio_INS, "\n")
cat("DEL overlap ratio:", ratio_DEL, "\n")