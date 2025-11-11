#!/usr/bin/env Rscript

# -----------------------------------------------------------------
# FileName     : SV_QTL_tag_by_GWASCatalog_lead_SNVs.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -----------------------------------------------------------------
# Description  : Identify SV QTLs that tag GWAS Catalog lead SNVs through LD analysis
# -----------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(vroom))
suppressPackageStartupMessages(library(pryr))
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
options(stringsAsFactors = FALSE)
mc <- getOption('mc.cores', 36)

gwas = fread("/home/wtian/data/GWAS_catalog/gwas_catalog_v1.0-associations_e114_r2025-07-10.tsv")
gwas = gwas %>% select(PUBMEDID,`DISEASE/TRAIT`, REGION, CHR_ID, CHR_POS,SNPS, `REPORTED GENE(S)`,`MAPPED_GENE`, `OR or BETA`,`P-VALUE`)
colnames(gwas) <- c("pubmed_id", "trait", "region", "chr", "pos", "snp", "reported_gene", "mapped_gene", "effect", "pval")
gwas$pval = as.numeric(gwas$pval)

gwas = gwas %>% filter(str_detect(snp, "rs|:"))
gwas = gwas %>% filter(!(str_detect(snp, "rs") & str_detect(snp, ":")))

rsid_snps1 = gwas$snp[str_detect(gwas$snp, "rs")]
chrpos_snps2 = gwas$snp[str_detect(gwas$snp, ":")]

rsid_snps1_bed = fread("rsid_snps1.bed", header = F)

gwas_rsid_1 = gwas %>% filter(str_detect(snp, "rs"))
gwas_rsid_1 = gwas_rsid_1 %>% filter(gwas_rsid_1$snp %in% rsid_snps1_bed$V4)
gwas_chrpos_2 = gwas %>% filter(str_detect(snp, ":"))
tmp = str_split(gwas_chrpos_2$snp, ":\\s*", simplify = T)
gwas_chrpos_2$snp_format = paste0(tmp[,1], "_", tmp[,2])

rsid_snps1_bed = rsid_snps1_bed[!duplicated(rsid_snps1_bed$V4),]
colnames(rsid_snps1_bed) = c("chr", "start", "end", "snp")
gwas_rsid_1 = left_join(gwas_rsid_1, rsid_snps1_bed, by = "snp")
gwas_rsid_1 = rename(gwas_rsid_1, chr.x  = "chr")
gwas_rsid_1$snp_format = paste0(gwas_rsid_1$chr, "_", gwas_rsid_1$start)
gwas_rsid_1$chr.y = NULL
gwas_rsid_1$start = NULL
gwas_rsid_1$end = NULL

gwas_flt_snp_ldformat = rbind(gwas_rsid_1, gwas_chrpos_2)
gwas_flt_snp_ldformat$snp_format = paste0(str_replace_all(gwas_flt_snp_ldformat$snp_format, "chr", ""), "_SNP")

fwrite(gwas_flt_snp_ldformat, file = "gwas_flt_snp_ldformat.tsv", row.names = F, col.names = T, quote = F, sep = "\t")

ld = data.frame(vroom("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/02_complexVar_SNP_LD/05_ld_process/01_ld_subsets/SNP_SV.ld", delim = ","))
ld = ld %>% filter(SNP_A %in% gwas_flt_snp_ldformat$snp_format & R2 > 0.75)
ld = ld %>% select(SNP_A, SNP_B, R2)
colnames(ld) = c("snp_format", "sv_format", "r2")

gwas_flt_snp_ld_sv = gwas_flt_snp_ldformat %>% filter(snp_format %in% ld$snp_format)
gwas_flt_snp_sv_ld = left_join(gwas_flt_snp_ld_sv, ld, by = "snp_format", relationship = "many-to-many")

# Add SV annotations
sv_info = fread("GRCh38_241031.SV.annovar.output.variant_function")
sv_info_df = data.frame(
    gene_region = sv_info$V1, gene = sv_info$V2, sv_type = str_replace_all(sv_info$V15, "<|>", ""),
    chr = str_replace_all(sv_info$V3, "chr", ""), start = sv_info$V4, end = sv_info$V5, sv_id = sv_info$V13
)
sv_info_df$sv_format = paste0(sv_info_df$chr, "_", sv_info_df$start,  "_", sv_info_df$sv_type)
dup_sv_format = sv_info_df$sv_format[duplicated(sv_info_df$sv_format)]
sv_info_df = sv_info_df[!sv_info_df$sv_format %in% dup_sv_format,]
gwas_flt_snp_ld_sv_info = left_join(gwas_flt_snp_sv_ld, sv_info_df, by = "sv_format")

fwrite(gwas_flt_snp_ld_sv_info, file = "gwas_flt_snp_ld_sv_info.tsv", row.names = F, col.names = T, quote = F, sep = "\t")

gwas_flt_snp_ld_sv_info = fread("gwas_flt_snp_ld_sv_info.tsv")
write.table(c(table(gwas_flt_snp_ld_sv_info$gene_region)), file = "gwas_flt_snp_ld_sv_info_gene_region.txt", row.names = T, col.names = F, quote = F, sep = "\t")

sv_qtl = fread("merged_cis_sv_qtl.txt.gz")
gwas_flt_snp_ld_sv_info_qtl = gwas_flt_snp_ld_sv_info %>% filter(sv_id %in% sv_qtl$Variant)
fwrite(gwas_flt_snp_ld_sv_info_qtl, file = "gwas_flt_snp_ld_sv_info_qtl.tsv", row.names = F, col.names = T, quote = F, sep = "\t")

gwas_flt_snp_ld_sv_info_qtl = fread("gwas_flt_snp_ld_sv_info_qtl.tsv")
write.table(c(table(gwas_flt_snp_ld_sv_info_qtl$gene_region)), file = "gwas_flt_snp_ld_sv_info_qtl_gene_region.txt", row.names = T, col.names = F, quote = F, sep = "\t")

sv_qtl_gwas_catalog = sv_qtl %>% filter(Variant %in% gwas_flt_snp_ld_sv_info_qtl$sv_id)

# Filter for regions of interest
interest_regions = c("exonic", "splicing")
gwas_flt_snp_ld_sv_info_qtl_interest = gwas_flt_snp_ld_sv_info_qtl[gwas_flt_snp_ld_sv_info_qtl$gene_region %in% interest_regions,]
gwas_flt_snp_ld_sv_info_qtl_interest = gwas_flt_snp_ld_sv_info_qtl_interest[gwas_flt_snp_ld_sv_info_qtl_interest$gene == gwas_flt_snp_ld_sv_info_qtl_interest$reported_gene,]

# Count SV occurrences across traits
sv_trait_count <- gwas_flt_snp_ld_sv_info_qtl_interest[, .(n_trait = uniqueN(trait)), by = sv_id]
sv_multi_trait <- sv_trait_count[n_trait > 1]

snv_qtl = data.frame(vroom("merged_cis_snv_qtl.txt.gz", delim = "\t"))
snv_qtl$chr_pos = paste0(str_split(snv_qtl$Variant,"_",simplify = T)[,1],"_",str_split(snv_qtl$Variant,"_",simplify = T)[,2])