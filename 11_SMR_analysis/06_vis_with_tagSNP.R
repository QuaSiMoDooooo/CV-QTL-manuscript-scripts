#!/usr/bin/env Rscript

# -------------
# FileName     : 06_vis_with_tagSNP.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Visualize SMR results with tag SNP information
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(vroom))

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
options(stringsAsFactors = FALSE)
mc <- getOption('mc.cores', 36)

# Read LD information
ld = fread("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/02_complexVar_SNP_LD/05_ld_process/01_ld_subsets/SNP_SV.ld", sep=",")
ld = ld[abs(ld$R2) > 0.75,]
ld$tagSNP_chrpos = paste0("chr",ld$CHR_A,"_",ld$BP_A)
ld_snp = data.frame(
    chr = paste0("chr",ld$CHR_A), pos = ld$BP_A, end = ld$BP_A + 1
)
fwrite(ld_snp, file = "06_ld_snp.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Read intersect results
ld_rsid = fread("06_ld_snp_rsid.bed")
colnames(ld_rsid) = c("chr", "start", "end", "rsid")
tag_rsid = unique(ld_rsid$rsid)
tag_snp = tag_rsid

# Load SMR results
smr = fread("05_smr_all_cis_flt.tsv")
top_snp = unique(smr$topSNP)
smr$tag_snp = ifelse(smr$topSNP %in% tag_snp, "tag", "non")

data = smr %>% select(trait, mol_phe, probeID, topSNP, topSNP_chr, topSNP_bp, tag_snp, p_SMR)
data$topSNP_chr = factor(data$topSNP_chr, levels = seq(1, 22))

# Build chromosome information table
chr_info <- data.frame(
  chr = paste0("chr", 1:22),
  chr_length = c(
    248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
    159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
    114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
    58617616, 64444167, 46709983, 50818468
  )
)

chr_info <- chr_info %>%
  mutate(chr_num = 1:22) %>%
  mutate(chr_start = lag(cumsum(chr_length), default = 0))

# Add normalized genomic positions
data <- data %>%
  mutate(chr = paste0("chr", topSNP_chr)) %>%
  left_join(chr_info, by = "chr") %>%
  mutate(
    genomic_pos = chr_start + topSNP_bp,
    shape_type = ifelse(tag_snp == "non", "circle", "triangle"),
    neg_log10_p = -log10(p_SMR)
  )
data$topSNP_chr = factor(data$topSNP_chr, levels = seq(1, 22))

# Set color mapping
molphe_colors <- c(
  "APA" = "#9467BD",
  "expression" = "#FF7F0E",
  "splicing" = "#D82728",
  "methylation" = "#B77554"
)

# Set shape mapping
shape_values <- c("circle" = 16, "triangle" = 17)

# Extract chromosome centers for x-axis labels
chr_labels <- chr_info %>%
  mutate(
    center = chr_start + chr_length / 2,
    chr_label = gsub("chr", "", chr)
  )

data$trait = factor(data$trait, levels = sort(unique(data$trait), decreasing = TRUE))

# Create plot
p <- ggplot(data, aes(x = genomic_pos, y = trait)) +
  geom_point(
    aes(size = neg_log10_p, color = mol_phe, shape = shape_type),
    position = position_jitter(width = 2e6, height = 0.3),
    alpha = 0.5
  ) +
  scale_shape_manual(
  values = shape_values,
  name = "Tag SNP",
  labels = c("circle" = "no", "triangle" = "yes")
  ) +
  scale_color_manual(values = molphe_colors) +
  scale_x_continuous(
    breaks = chr_labels$center,
    labels = chr_labels$chr_label,
    expand = c(0.01, 0)
  ) +
  theme_bw() +
  labs(
    x = "Genomic Position (by Chromosome)",
    y = "Trait",
    size = expression(-log[10](SMR~P)),
    shape = "Tag SNP",
    color = "Molecular Phenotype"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

ggsave("06_vis_with_tagSNP_1.pdf", plot = p, width = 10, height = 10)