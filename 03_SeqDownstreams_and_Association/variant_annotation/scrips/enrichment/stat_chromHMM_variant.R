#! Rscript
# -------------
# FileName     : stat_chromHMM
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-12-31 22:35
# Last Modified: 2025-01-06 23:04
# Modified By  : EastsunW
# -------------
# Description  : 整合、可视化ChromHMM的富集结果
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
source("https://gitee.com/eastsunw/personal_code_notebook/raw/master/plot_tools/wdy_theme.r")
custom_colors <- c(
    "#FA7F6F",
    "#8ECFC9",
    "#FFBE7A",
    "#82B0D2",
    "#BEB8DC",
    "#e9e9e9",
    "#999999"
)
setwd("/home/wangdy/Projects/Weibin/Downstreams/variant_annotation")

variant_types <- c("SNP", "InDel", "MNV", "SV")
annotations <- c(
    "Promoter", "Enhancer", "Transcriped", "ZNF_repeat",
    "Heterchromatin", "Polycomb", "Quiescent"
)

# 细胞系对应组织
cell_line_annotation <- c(
    "E112" = "Thymus--Thymus",
    "E093" = "Thymus--Fetal Thymus",
    "E078" = "Sm.Muscle--Duodenum Smooth Muscle",
    "E076" = "Sm.Muscle--Colon Smooth Muscle",
    "E103" = "Sm.Muscle--Rectal Smooth Muscle",
    "E111" = "Sm.Muscle--Stomach Smooth Muscle",
    "E099" = "Other--Placenta Amnion",
    "E086" = "Other--Fetal Kidney",
    "E088" = "Other--Fetal Lung",
    "E097" = "Other--Ovary",
    "E087" = "Other--Pancreatic lslets",
    "E080" = "Other--Fetal Adrenal Gland",
    "E091" = "Other--Placenta",
    "E066" = "Other--Liver",
    "E098" = "Other--Pancreas",
    "E096" = "Other--Lung",
    "E113" = "Other--Spleen",
    "E054" = "Neurosph--Ganglion Eminence derived primary cultured neurospheres",
    "E053" = "Neurosph--Cortex derived primary cultured neurospheres",
    "E052" = "Myosat--Muscle Satellite Cultured Cells",
    "E100" = "Muscle--Psoas Muscle",
    "E108" = "Muscle--Skeletal Muscle Female",
    "E107" = "Muscle--Skeletal Muscle Male",
    "E089" = "Muscle--Fetal Muscle Trunk",
    "E090" = "Muscle--Fetal Muscle Leg",
    "E026" = "Mesench--Bone Marrow Derived Cultured Mesenchymal Stem Cells",
    "E049" = "Mesench--Mesenchymal Stem Cell Derived Chondrocyte Cultured Cells",
    "E025" = "Mesench--Adipose Derived Mesenchymal Stem Cell Cultured Cells",
    "E023" = "Mesench--Mesenchymal Stem Cell Derived Adipocyte Cultured Cells",
    "E020" = "iPSC--iPS-20b cells",
    "E019" = "iPSC--iPS-18 Cells",
    "E018" = "iPSC--iPS-15b Cells",
    "E021" = "iPSC--iPS DF 6.9 Cells",
    "E022" = "iPSC--iPS DF 19.11 Cells",
    "E017" = "IMR90--IMR90 fetal lung fibroblasts Cell Line",
    "E029" = "HSC&B-cell--Primary monocytes from peripheral blood",
    "E031" = "HSC&B-cell--Primary B cells from cord blood",
    "E035" = "HSC&B-cell--Primary hematopoietic stem cells",
    "E051" = "HSC&B-cell--Primary hematopoietic stem cells G-cSF-mobilized Male",
    "E050" = "HSC&B-cell--Primary hematopoietic stem cells G-CSF-mobilized Female",
    "E036" = "HSC&B-cell--Primary hematopoietic stem cells short term culture",
    "E032" = "HSC&B-cell--Primary B cells from peripheral blood",
    "E046" = "HSC&B-cell--Primary Natural Killer cells from peripheral blood",
    "E030" = "HSC&B-cell--Primary neutrophils from peripheral blood",
    "E083" = "Heart--Fetal Heart",
    "E104" = "Heart--Right Atrium",
    "E095" = "Heart--Left Ventricle",
    "E105" = "Heart--Right Ventricle",
    "E065" = "Heart--Aorta",
    "E002" = "ESC--ES-WA7 Cells",
    "E008" = "ESC--H9 Cells",
    "E001" = "ESC--ES-3 Cells",
    "E015" = "ESC--HUES6 Cells",
    "E014" = "ESC--HUES48 Cells",
    "E016" = "ESC--HUES64 Cells",
    "E003" = "ESC--H1 Cells",
    "E024" = "ES-deriv--ES-UCSF4 Cells",
    "E007" = "ES-deriv--H1 Derived Neuronal Progenitor Cultured Cells",
    "E009" = "ES-deriv--H9 Derived Neuronal Progenitor cultured cells",
    "E010" = "ES-deriv--H9 Derived Neuron Cultured Cells",
    "E013" = "ES-deriv--hESC Derived CD56+ Mesoderm Cultured Cells",
    "E012" = "ES-deriv--hESC Derived CD56+ Ectoderm Cultured Cells",
    "E011" = "ES-deriv--hESC Derived CD184+ Endodemm Cultured Cells",
    "E004" = "ES-deriv--H1 BMP4 Derived Mesendoderm Cultured cells",
    "E005" = "ES-deriv--H1 BMP4 Derived Trophoblast Cultured Cells",
    "E006" = "ES-deriv--H1 Derived Mesenchymal Stem Cells",
    "E055" = "Epithelial--Foreskin Fibroblast Primary Cells skin01",
    "E056" = "Epithelial--Foreskin Fibroblast Primary Cells skin02",
    "E059" = "Epithelial--Foreskin Melanocyte Primary cells skin01",
    "E061" = "Epithelial--Foreskin Melanocyte Primary Cells skin03",
    "E057" = "Epithelial--Foreskin Keratinocyte Primary cells skin02",
    "E058" = "Epithelial--Foreskin Keratinocyte Primary Cells skin03",
    "E028" = "Epithelial--Breast variant Human Mammary Epithelial Cells (vHMEC)",
    "E027" = "Epithelial--Breast Myoepithelial Primary Cells",
    "E114" = "ENCODE--A549 EtOH 0.02pct Lung carcinoma Cell Line",
    "E115" = "ENCODE--Dnd41 TCell Leukemia cell Line",
    "E116" = "ENCODE--GM12878 Lymphoblastoid Cells",
    "E117" = "ENCODE--HeLa-S3 Cervical carcinoma cell Line",
    "E118" = "ENCODE--HepG2 Hepatocellular Carcinoma Cell Line",
    "E119" = "ENCODE--HMEC Mammary Epithelial Primary Cells",
    "E120" = "ENCODE--HSMM Skeletal Muscle Myoblasts Cells",
    "E121" = "ENCODE--HSMM cell derived Skeletal Muscle Myotubes Cells",
    "E122" = "ENCODE--HUVEC Umbilical Vein Endothelial Primary cells",
    "E123" = "ENCODE--K562 Leukemia cells",
    "E124" = "ENCODE--Monocytes-CD14+ RO01746 Primary cells",
    "E125" = "ENCODE--NH-A Astrocytes Primary Cells",
    "E126" = "ENCODE--NHDF-Ad Adult Dermal Fibroblast Primary Cells",
    "E127" = "ENCODE--NHEK-Epidermal Keratinocyte Primary Cells",
    "E128" = "ENCODE--NHLF Lung Fibroblast Primary Cells",
    "E129" = "ENCODE--Osteoblast Primary cells",
    "E092" = "Digestive--Fetal stomach",
    "E085" = "Digestive--Fetal Intestine Small",
    "E084" = "Digestive--Fetal Intestine Large",
    "E109" = "Digestive--Small Intestine",
    "E106" = "Digestive--Sigmoid Colon",
    "E075" = "Digestive--Colonic Mucosa",
    "E101" = "Digestive--Rectal Mucosa Donor 29",
    "E102" = "Digestive--Rectal Mucosa Donor 31",
    "E110" = "Digestive--Stomach Mucosa",
    "E077" = "Digestive--Duodenum Mucosa",
    "E079" = "Digestive--Esophagus",
    "E094" = "Digestive--Gastric",
    "E071" = "Brain--Brain Hippocampus Middle",
    "E074" = "Brain--Brain Substantia Nigra",
    "E068" = "Brain--Brain Anterior Caudate",
    "E069" = "Brain--Brain Cingulate Gyrus",
    "E072" = "Brain--Brain Inferior Temporal Lobe",
    "E067" = "Brain--Brain Angular Gyrus",
    "E073" = "Brain--Brain Dorsolateral Prefrontal Cortex",
    "E070" = "Brain--Brain Germinal Matrix",
    "E082" = "Brain--Fetal Brain Female",
    "E081" = "Brain--Fetal Brain Male",
    "E062" = "Blood--Primary mononuclear cells from peripheral blood",
    "E063" = "Adipose--Adipose Nuclei",
    "E034" = "Blood--Primary T cells from peripheral blood",
    "E045" = "Blood--Primary T cells effector/memory enriched from peripheral blood",
    "E033" = "Blood--Primary T cells from cord blood",
    "E044" = "Blood--Primary T regulatory cells from peripheral blood",
    "E043" = "Blood--Primary T helper cells from peripheral blood",
    "E039" = "Blood--Primary T helper naive cells from peripheral blood",
    "E041" = "Blood--Primary T helper cells PMA-l stimulated",
    "E042" = "Blood--Primary T helper 17 cells PMA-l stimulated",
    "E040" = "Blood--Primary T helper memory cells from peripheral blood 1",
    "E037" = "Blood--Primary T helper memory cells from peripheral blood 2",
    "E048" = "Blood--Primary T CD8+ memory cells from peripheral blood",
    "E038" = "Blood--Primary T helper naive cells from peripheral blood",
    "E047" = "Blood--PrimaryT CD8+ naive cells from peripheral blood"
)

cell_line_annotation_df <- data.frame(
    cell_line = names(cell_line_annotation),
    description = as.character(cell_line_annotation)
) %>%
    separate(
        col = "description",
        into = c("category", "description"),
        sep = "--"
    ) %>%
    arrange(cell_line) %>%
    unite(cell_line, description, col = "description", sep = ": ", remove = FALSE)

data_dir <- "results/epi_enrichment/chromHMM"
out_dir <- "results/epi_enrichment/chromHMM/stat_variant"
if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}

merged_cell_line_enrich <- list()
for (variant_type in variant_types) {
    if (variant_type == "SV") {
        freq_types <- c("all", "common", "rare", "INS", "DEL", "DUP", "INV")
    } else {
        freq_types <- c("all", "common", "rare")
    }
    for (freq_type in freq_types) {
        merged_data <- data.frame(
            annotation = annotations
        )
        for (n in 1:129) {
            cell_line <- paste0("E", str_pad(n, 3, pad = "0"))
            result_path <- file.path(
                data_dir,
                "variant",
                paste0(
                    cell_line, "_",
                    variant_type, "_",
                    freq_type,
                    ".results.txt"
                )
            )
            if (file.exists(result_path)) {
                result <- fread(result_path, sep = "\t") %>%
                    select(annotation, l2fold, qvalue) %>%
                    unite(l2fold, qvalue, col = "stat", sep = "__") %>%
                    rename(!!cell_line := stat)
                merged_data <<- left_join(
                    merged_data,
                    result,
                    by = "annotation"
                )
                merged_cell_line_enrich[[paste0(variant_type, "_", freq_type)]] <- merged_data
            }
        }
    }
}

# 处理成ggplot的数据格式，第一列是变异类型（列表的name），第二列是annotation，第三列是cell_line，第四列是p_sig
plot_df <- data.frame()
temp <- sapply(
    X = names(merged_cell_line_enrich),
    FUN = function(name) {
        data <- merged_cell_line_enrich[[name]] %>%
            mutate(variant_type = name) %>%
            separate(variant_type, c("variant_type", "freq_type"), sep = "_") %>%
            pivot_longer(
                cols = -c(annotation, variant_type, freq_type),
                names_to = "cell_line",
                values_to = "stat"
            ) %>%
            separate(stat, c("log2FC", "FDR"), sep = "__") %>%
            mutate_at(vars(log2FC, FDR), as.numeric)
        plot_df <<- rbind(plot_df, data)
    }
)
rm(temp)
plot_df_adj <- plot_df %>%
    left_join(
        cell_line_annotation_df,
        by = c("cell_line" = "cell_line")
    ) %>%
    mutate(p_sig = case_when(
        FDR >= 0.05 ~ NA,
        FDR < 0.05 ~ "*",
    )) %>%
    relocate(cell_line, category, description, .before = annotation) %>%
    mutate(
        annotation = factor(annotation, levels = annotations),
        variant_type = factor(variant_type, levels = variant_types),
        freq_type = factor(freq_type, levels = freq_types)
    )

fwrite(
    x = plot_df_adj,
    file = file.path(out_dir, "enrich_chromHMM_variant.txt"),
    sep = "\t",
    quote = FALSE
)

plot_by_annotation_all <- ggplot(
    data = plot_df_adj %>% filter(freq_type == "all"),
    mapping = aes(
        x = cell_line,
        y = variant_type,
        fill = log2FC
    )
) +
    geom_tile() +
    scale_fill_gradientn(
        limits = range(-3, 3),
        oob = scales::squish,
        labels = c("-3", "", "", "0", "", "", "3"),
        breaks = c(-3, -2, -1, 0, 1, 2, 3),
        colors = c("#122A64", "#223290", "#3085BC", "#ffffff", "#F56B46", "#B8070D", "#670519"),
        na.value = "#ffffff"
    ) +
    facet_wrap(
        ~annotation,
        ncol = 1,
        strip.position = "right"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
        x = "Cell line",
        y = "Variant type",
        fill = "log2FC"
    ) +
    wdy_theme(
        base_size = 10,
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.y = element_text(angle = 0),
        panel.spacing.x = unit(2, "pt"),
        panel.spacing.y = unit(2, "pt"),
        strip.background = element_rect(fill = "#ebf1f5", linewidth = NA)
    )

plot_by_variant_all <- ggplot(
    data = plot_df_adj %>% filter(freq_type == "all"),
    mapping = aes(
        x = cell_line,
        y = annotation,
        fill = log2FC
    )
) +
    geom_tile() +
    scale_fill_gradientn(
        limits = range(-3, 3),
        oob = scales::squish,
        labels = c("-3", "", "", "0", "", "", "3"),
        breaks = c(-3, -2, -1, 0, 1, 2, 3),
        colors = c("#122A64", "#223290", "#3085BC", "#ffffff", "#F56B46", "#B8070D", "#670519"),
        na.value = "#ffffff"
    ) +
    facet_wrap(
        ~variant_type,
        ncol = 1,
        strip.position = "right"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
        x = "Cell line",
        y = "Annotation",
        fill = "log2FC"
    ) +
    wdy_theme(
        base_size = 10,
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.y = element_text(angle = 0),
        panel.spacing.x = unit(2, "pt"),
        panel.spacing.y = unit(2, "pt"),
        strip.background = element_rect(fill = "#ebf1f5", linewidth = NA)
    )

plot_by_annotation_freq <- ggplot(
    data = plot_df_adj %>% filter(freq_type %in% c("common", "rare")),
    mapping = aes(
        x = cell_line,
        y = variant_type,
        fill = log2FC
    )
) +
    geom_tile() +
    scale_fill_gradientn(
        limits = range(-3, 3),
        oob = scales::squish,
        labels = c("-3", "", "", "0", "", "", "3"),
        breaks = c(-3, -2, -1, 0, 1, 2, 3),
        colors = c("#122A64", "#223290", "#3085BC", "#ffffff", "#F56B46", "#B8070D", "#670519"),
        na.value = "#ffffff"
    ) +
    facet_grid(
        rows = vars(annotation),
        cols = vars(freq_type)
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
        x = "Cell line",
        y = "Variant type",
        fill = "log2FC"
    ) +
    wdy_theme(
        base_size = 10,
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.y = element_text(angle = 0),
        panel.spacing.x = unit(2, "pt"),
        panel.spacing.y = unit(2, "pt"),
        strip.background = element_rect(fill = "#ebf1f5", linewidth = NA)
    )

plot_by_SV_by_annotation <- ggplot(
    data = plot_df_adj %>% filter(freq_type %in% c("INS", "DEL", "DUP", "INV")),
    mapping = aes(
        x = cell_line,
        y = freq_type,
        fill = log2FC
    )
) +
    geom_tile() +
    scale_fill_gradientn(
        limits = range(-3, 3),
        oob = scales::squish,
        labels = c("-3", "", "", "0", "", "", "3"),
        breaks = c(-3, -2, -1, 0, 1, 2, 3),
        colors = c("#122A64", "#223290", "#3085BC", "#ffffff", "#F56B46", "#B8070D", "#670519"),
        na.value = "#ffffff"
    ) +
    facet_wrap(
        ~annotation,
        ncol = 1,
        strip.position = "right"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
        x = "Cell line",
        y = "SV type",
        fill = "log2FC"
    ) +
    wdy_theme(
        base_size = 10,
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.y = element_text(angle = 0),
        panel.spacing.x = unit(2, "pt"),
        panel.spacing.y = unit(2, "pt"),
        strip.background = element_rect(fill = "#ebf1f5", linewidth = NA)
    )

plot_by_SV_by_svtype <- ggplot(
    data = plot_df_adj %>% filter(freq_type %in% c("INS", "DEL", "DUP", "INV")),
    mapping = aes(
        x = cell_line,
        y = annotation,
        fill = log2FC
    )
) +
    geom_tile() +
    scale_fill_gradientn(
        limits = range(-3, 3),
        oob = scales::squish,
        labels = c("-3", "", "", "0", "", "", "3"),
        breaks = c(-3, -2, -1, 0, 1, 2, 3),
        colors = c("#122A64", "#223290", "#3085BC", "#ffffff", "#F56B46", "#B8070D", "#670519"),
        na.value = "#ffffff"
    ) +
    facet_wrap(
        ~freq_type,
        ncol = 1,
        strip.position = "right"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
        x = "Cell line",
        y = "Annotation",
        fill = "log2FC"
    ) +
    wdy_theme(
        base_size = 10,
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.y = element_text(angle = 0),
        panel.spacing.x = unit(2, "pt"),
        panel.spacing.y = unit(2, "pt"),
        strip.background = element_rect(fill = "#ebf1f5", linewidth = NA)
    )

plot_by_annotation_all_blood <- ggplot(
    data = plot_df_adj %>% filter(freq_type == "all", category == "Blood"),
    mapping = aes(
        x = cell_line,
        y = variant_type,
        fill = log2FC
    )
) +
    geom_tile() +
    scale_fill_gradientn(
        limits = range(-3, 3),
        oob = scales::squish,
        labels = c("-3", "", "", "0", "", "", "3"),
        breaks = c(-3, -2, -1, 0, 1, 2, 3),
        colors = c("#122A64", "#223290", "#3085BC", "#ffffff", "#F56B46", "#B8070D", "#670519"),
        na.value = "#ffffff"
    ) +
    geom_text(
        aes(label = p_sig),
        size = 10 / .pt,
        color = "black",
        hjust = 0.5,
        vjust = 0.5
    ) +
    facet_wrap(
        ~annotation,
        nrow = 1,
        strip.position = "top"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
        x = "Cell line",
        y = "Variant type",
        fill = "log2FC"
    ) +
    coord_flip() +
    wdy_theme(
        base_size = 10,
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.text.y = element_text(angle = 0),
        panel.spacing.x = unit(2, "pt"),
        panel.spacing.y = unit(2, "pt"),
        strip.background = element_rect(fill = "#ebf1f5", linewidth = NA)
    )

plot_by_variant_all_blood <- ggplot(
    data = plot_df_adj %>% filter(freq_type == "all", category == "Blood"),
    mapping = aes(
        x = cell_line,
        y = annotation,
        fill = log2FC
    )
) +
    geom_tile() +
    scale_fill_gradientn(
        limits = range(-3, 3),
        oob = scales::squish,
        labels = c("-3", "", "", "0", "", "", "3"),
        breaks = c(-3, -2, -1, 0, 1, 2, 3),
        colors = c("#122A64", "#223290", "#3085BC", "#ffffff", "#F56B46", "#B8070D", "#670519"),
        na.value = "#ffffff"
    ) +
    geom_text(
        aes(label = p_sig),
        size = 10 / .pt,
        color = "black",
        hjust = 0.5,
        vjust = 0.5
    ) +
    facet_wrap(
        ~variant_type,
        nrow = 1,
        strip.position = "top"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
        x = "Cell line",
        y = "Annotation",
        fill = "log2FC"
    ) +
    coord_flip() +
    wdy_theme(
        base_size = 10,
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.text.y = element_text(angle = 0),
        panel.spacing.x = unit(2, "pt"),
        panel.spacing.y = unit(2, "pt"),
        strip.background = element_rect(fill = "#ebf1f5", linewidth = NA)
    )

plot_by_annotation_freq_blood <- ggplot(
    data = plot_df_adj %>% filter(freq_type %in% c("common", "rare"), category == "Blood"),
    mapping = aes(
        x = cell_line,
        y = variant_type,
        fill = log2FC
    )
) +
    geom_tile() +
    scale_fill_gradientn(
        limits = range(-3, 3),
        oob = scales::squish,
        labels = c("-3", "", "", "0", "", "", "3"),
        breaks = c(-3, -2, -1, 0, 1, 2, 3),
        colors = c("#122A64", "#223290", "#3085BC", "#ffffff", "#F56B46", "#B8070D", "#670519"),
        na.value = "#ffffff"
    ) +
    geom_text(
        aes(label = p_sig),
        size = 10 / .pt,
        color = "black",
        hjust = 0.5,
        vjust = 0.5
    ) +
    facet_grid(
        rows = vars(freq_type),
        cols = vars(annotation)
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
        x = "Cell line",
        y = "Variant type",
        fill = "log2FC"
    ) +
    coord_flip() +
    wdy_theme(
        base_size = 10,
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.text.y = element_text(angle = 0),
        panel.spacing.x = unit(2, "pt"),
        panel.spacing.y = unit(2, "pt"),
        strip.background = element_rect(fill = "#ebf1f5", linewidth = NA)
    )

plot_by_SV_by_annotation_blood <- ggplot(
    data = plot_df_adj %>% filter(freq_type %in% c("INS", "DEL", "DUP", "INV"), category == "Blood"),
    mapping = aes(
        x = cell_line,
        y = freq_type,
        fill = log2FC
    )
) +
    geom_tile() +
    scale_fill_gradientn(
        limits = range(-3, 3),
        oob = scales::squish,
        labels = c("-3", "", "", "0", "", "", "3"),
        breaks = c(-3, -2, -1, 0, 1, 2, 3),
        colors = c("#122A64", "#223290", "#3085BC", "#ffffff", "#F56B46", "#B8070D", "#670519"),
        na.value = "#ffffff"
    ) +
    geom_text(
        aes(label = p_sig),
        size = 10 / .pt,
        color = "black",
        hjust = 0.5,
        vjust = 0.5
    ) +
    facet_wrap(
        ~annotation,
        nrow = 1,
        strip.position = "top"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
        x = "Cell line",
        y = "SV type",
        fill = "log2FC"
    ) +
    coord_flip() +
    wdy_theme(
        base_size = 10,
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.text.y = element_text(angle = 0),
        panel.spacing.x = unit(2, "pt"),
        panel.spacing.y = unit(2, "pt"),
        strip.background = element_rect(fill = "#ebf1f5", linewidth = NA)
    )

plot_by_SV_by_svtype_blood <- ggplot(
    data = plot_df_adj %>% filter(freq_type %in% c("INS", "DEL", "DUP", "INV"), category == "Blood"),
    mapping = aes(
        x = cell_line,
        y = annotation,
        fill = log2FC
    )
) +
    geom_tile() +
    scale_fill_gradientn(
        limits = range(-3, 3),
        oob = scales::squish,
        labels = c("-3", "", "", "0", "", "", "3"),
        breaks = c(-3, -2, -1, 0, 1, 2, 3),
        colors = c("#122A64", "#223290", "#3085BC", "#ffffff", "#F56B46", "#B8070D", "#670519"),
        na.value = "#ffffff"
    ) +
    geom_text(
        aes(label = p_sig),
        size = 10 / .pt,
        color = "black",
        hjust = 0.5,
        vjust = 0.5
    ) +
    facet_wrap(
        ~freq_type,
        nrow = 1,
        strip.position = "top"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
        x = "Cell line",
        y = "Annotation",
        fill = "log2FC"
    ) +
    coord_flip() +
    wdy_theme(
        base_size = 10,
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.text.y = element_text(angle = 0),
        panel.spacing.x = unit(2, "pt"),
        panel.spacing.y = unit(2, "pt"),
        strip.background = element_rect(fill = "#ebf1f5", linewidth = NA)
    )

ggsave(
    file = file.path(out_dir, "plot_by_annotation_all.pdf"),
    plot = plot_by_annotation_all,
    width = 8,
    height = 4
)
ggsave(
    file = file.path(out_dir, "plot_by_variant_all.pdf"),
    plot = plot_by_variant_all,
    width = 8,
    height = 4
)
ggsave(
    file = file.path(out_dir, "plot_by_annotation_freq.pdf"),
    plot = plot_by_annotation_freq,
    width = 10,
    height = 5
)
ggsave(
    file = file.path(out_dir, "plot_by_SV_by_annotation.pdf"),
    plot = plot_by_SV_by_annotation,
    width = 8,
    height = 5
)
ggsave(
    file = file.path(out_dir, "plot_by_SV_by_svtype.pdf"),
    plot = plot_by_SV_by_svtype,
    width = 8,
    height = 5
)

ggsave(
    file = file.path(out_dir, "plot_by_annotation_all_blood.pdf"),
    plot = plot_by_annotation_all_blood,
    width = 9,
    height = 3
)
ggsave(
    file = file.path(out_dir, "plot_by_variant_all_blood.pdf"),
    plot = plot_by_variant_all_blood,
    width = 6,
    height = 3.5
)
ggsave(
    file = file.path(out_dir, "plot_by_annotation_freq_blood.pdf"),
    plot = plot_by_annotation_freq_blood,
    width = 10,
    height = 5
)
ggsave(
    file = file.path(out_dir, "plot_by_SV_by_annotation_blood.pdf"),
    plot = plot_by_SV_by_annotation_blood,
    width = 9,
    height = 3
)
ggsave(
    file = file.path(out_dir, "plot_by_SV_by_svtype_blood.pdf"),
    plot = plot_by_SV_by_svtype_blood,
    width = 8,
    height = 4
)
