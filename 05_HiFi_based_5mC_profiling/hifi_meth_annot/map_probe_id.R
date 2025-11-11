#! Rscript
# -------------
# FileName     : map_probe_id
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : map probe id to the HiFi-based 5mC sites
# -------------

suppressPackageStartupMessages(library(ChAMP))
suppressPackageStartupMessages(library(easylift))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
options(scipen=999)

# get probe info
chain <- "/home/wtian/R/x86_64-pc-linux-gnu-library/4.3/easylift/extdata/hg19ToHg38.over.chain.gz"
# 450k
data(hm450.manifest.hg19)
df450 <- hm450.manifest.hg19
seqname_vec <- df450$CpG_chrm
start_vec <- df450$CpG_beg+1
end_vec <- df450$CpG_end
probe_id <- rownames(df450)
gr <- GRanges(seqname=seqname_vec,ranges = IRanges(start=start_vec,end = end_vec),probe_id=probe_id)
genome(gr) <- "hg19"
df450_hg38 <- data.frame(easylift(gr, "hg38", chain))
colnames(df450_hg38) <- c("chr","hg38_start","hg38_end","width","strand","probe_id")
df450$probe_id <- rownames(df450)
hm450.manifest.hg38 <- inner_join(df450_hg38, df450, by="probe_id")
hm450.manifest.hg38 <- hm450.manifest.hg38 %>% select(probe_id, everything())
write.table(hm450.manifest.hg38, "hm450.manifest.hg38.tsv", sep="\t", quote=F, row.names=F,col.names = T)
# 850k
data(EPIC.manifest.hg19)
df850 <- EPIC.manifest.hg19
seqname_vec <- df850$CpG_chrm
start_vec <- df850$CpG_beg+1
end_vec <- df850$CpG_end
probe_id <- rownames(df850)
gr <- GRanges(seqname=seqname_vec,ranges = IRanges(start=start_vec,end = end_vec),probe_id=probe_id)
genome(gr) <- "hg19"
df850_hg38 <- data.frame(easylift(gr, "hg38", chain))
colnames(df850_hg38) <- c("chr","hg38_start","hg38_end","width","strand","probe_id")
df850$probe_id <- rownames(df850)
EPIC.manifest.hg38 <- inner_join(df850_hg38, df850, by="probe_id")
EPIC.manifest.hg38 <- EPIC.manifest.hg38 %>% select(probe_id, everything())
write.table(EPIC.manifest.hg38, "EPIC.manifest.hg38.tsv", sep="\t", quote=F, row.names=F,col.names = T)
# merge
merge_df <- rbind(df450, df850)
merge_df <- merge_df[!duplicated(merge_df$probeID),]
merge_df = merge_df[,c("probeID","CpG_chrm","CpG_beg")]
merge_df$chr = str_split(merge_df$CpG_chrm,"chr",simplify=T)[,2]
merge_df$pos= merge_df$CpG_beg+1
merge_df = merge_df[,c("chr","pos","probeID")]
write.table(merge_df, "merge_manifest.hg19.tsv", sep="\t", quote=F, row.names=F,col.names = T)
merge_df <- rbind(hm450.manifest.hg38, EPIC.manifest.hg38)
merge_df <- merge_df[!duplicated(merge_df$probe_id),]
write.table(merge_df, "merge_manifest.hg38.tsv", sep="\t", quote=F, row.names=F,col.names = T)






