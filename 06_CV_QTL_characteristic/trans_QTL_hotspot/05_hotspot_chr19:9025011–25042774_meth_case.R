#!/usr/bin/env Rscript

# -------------
# FileName     : 05_hotspot_chr19:9025011–25042774_meth_case.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Analyze methylation QTLs in chromosome 19 hotspot region (9025011-25042774) and identify promoter overlaps
# -------------

setwd("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/01_quick_stats_for_paper/trans_hotspot_loci")
library(vroom)
library(dplyr)
library(stringr)

qtl_dir <- "~/project/HZAU_cohort_meth/wdy_result/phenotype_qtl_results_bim3_removaNA"
variants <- c("SNP")

hot_chr <- "19"
hot_start <- 9025011
hot_end   <- 25042774

all_records <- data.frame()

for (var in variants) {
  qtl_file <- file.path(qtl_dir, paste0(var, "-methylation"), "QTL_results", "trans.filtered.txt.gz")
  if (!file.exists(qtl_file)) next
  
  message("Processing: ", qtl_file)
  
  qtl <- vroom(qtl_file, col_select = c("Variant_position","Phenotype","Phenotype_position")) %>%
    as.data.frame()
  
  coords <- str_split_fixed(qtl$Variant_position, "[:-]", 3)
  qtl$chr <- str_remove(coords[,1], "chr")
  qtl$pos <- as.integer(coords[,2])
  
  in_hot <- qtl[qtl$chr == hot_chr & qtl$pos >= hot_start & qtl$pos <= hot_end, 
                c("Phenotype","Phenotype_position")]
  
  if (nrow(in_hot) > 0) {
    all_records <- rbind(all_records, in_hot)
  }
}

all_records <- all_records[!duplicated(all_records$Phenotype), ]

out_file <- "methylation_chr19_hotspot_phenotypes_with_positions.tsv"
write.table(all_records, out_file, sep = "\t", row.names = FALSE, quote = FALSE)

message("Done! Found ", nrow(all_records), 
        " unique methylation phenotypes with positions. Saved to: ", out_file)

library(parallel)
library(GenomicRanges)

gene_file <- "~/data/annot/GRCh38/GENCODE/release-46/hg38/gencode.v46.annotation.genesymbol_start_strand.csv"
genes <- vroom(gene_file, delim = "\t") %>% as.data.frame()

promoters <- mclapply(1:nrow(genes), function(i) {
  g <- genes[i, ]
  start_pos <- as.integer(g$start)
  if (g$strand == "+") {
    tss <- start_pos
    prom_start <- max(1, tss - 7500)
    prom_end   <- tss + 2500
  } else {
    tss <- start_pos
    prom_start <- max(1, tss - 2500)
    prom_end   <- tss + 7500
  }
  GRanges(seqnames = g$chr,
          ranges = IRanges(prom_start, prom_end),
          strand = g$strand,
          gene_symbol = g$gene_symbol)
}, mc.cores = 60)

promoters <- do.call(c, promoters)
promoters_reduced <- reduce(promoters, with.revmap = TRUE)

meth_file <- "methylation_chr19_hotspot_phenotypes_with_positions.tsv"
meth <- vroom(meth_file) %>% as.data.frame()

coords <- str_split_fixed(meth$Phenotype_position, "[:-]", 4)
meth_gr <- GRanges(
  seqnames = coords[,1],
  ranges = IRanges(as.integer(coords[,2]), as.integer(coords[,3])),
  strand = coords[,4],
  Phenotype = meth$Phenotype
)

hits <- findOverlaps(meth_gr, promoters_reduced)

result <- data.frame(
  Phenotype = meth_gr$Phenotype[queryHits(hits)],
  Meth_Position = paste0(seqnames(meth_gr)[queryHits(hits)], ":",
                         start(meth_gr)[queryHits(hits)], "-",
                         end(meth_gr)[queryHits(hits)], ":",
                         strand(meth_gr)[queryHits(hits)]),
  Promoter_Region = paste0(seqnames(promoters_reduced)[subjectHits(hits)], ":",
                           start(promoters_reduced)[subjectHits(hits)], "-",
                           end(promoters_reduced)[subjectHits(hits)])
)

out_file <- "methylation_chr19_promoter_overlap.tsv"
write.table(result, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
message("Done! Found ", nrow(result), " methylation sites in promoter regions. Saved to: ", out_file)

Promoter_Region = unique(sort(result$Promoter_Region))
tmp = str_split_fixed(Promoter_Region, "[:-]", 3)
tmp = data.frame(chr = str_remove(tmp[,1], "chr"), start = as.integer(tmp[,2]), end = as.integer(tmp[,3]))
write.table(tmp, "promoter_region.bed", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)