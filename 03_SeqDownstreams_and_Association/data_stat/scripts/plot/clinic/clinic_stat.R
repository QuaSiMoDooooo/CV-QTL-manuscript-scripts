#! Rscript
# -------------
# FileName     : clinic_stat
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2025-02-28 22:44
# Last Modified: 2025-03-04 12:25
# Modified By  : EastsunW
# -------------
# Description  : 统计各种分类指标的各种范围和风险的数量
# -------------

# 加载必要的包
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(viridis))
# 绘图主题
source("https://gitee.com/eastsunw/personal_code_notebook/raw/master/plot_tools/wdy_theme.r")
setwd("/home/wangdy/Projects/Weibin/Downstreams/data_stat")
plot_dir <- "results/plot/clinic"
stat_dir <- "results/stat/clinic"

level_mapping <- c(
    "过低" = "Too Low",
    "低" = "Low",
    "偏低" = "Slightly Low",
    "正常" = "Normal",
    "偏高" = "Slightly High",
    "高" = "High",
    "过高" = "Too High"
)

risk_mapping <- c(
    "正常" = "Normal",
    "低风险" = "Low Risk",
    "中风险" = "Medium Risk",
    "高风险" = "High Risk"
)

category_mapping <- c(
    "测量类" = "Measurement",
    "营养类" = "Nutrition",
    "病理类" = "Pathology",
    "环境毒素" = "Toxin",
    "激素类" = "Hormone",
    "影像类" = "Imaging",
    "问卷类" = "Questionnaire"
)


# 读取指标的信息，用来获取指标的类型
bio_info_path <- "data/Clinic/bio_info.txt"
bio_level_path <- "data/Clinic/bio_level.txt"
bio_risk_path <- "data/Clinic/bio_risk.txt"

bio_info <- fread(bio_info_path, sep = "\t") %>% select(Indicator, Category)

bio_level <- fread(bio_level_path, sep = "\t") %>%
    column_to_rownames("Indicator") %>%
    filter(!if_all(.cols = everything(), .fns = is.na)) %>%
    rownames_to_column("Indicator")

bio_level_stat <- bio_level %>%
    pivot_longer(cols = -Indicator, names_to = "Sample", values_to = "Level") %>%
    filter(!Level %in% c("阴性", "阳性") & !is.na(Level)) %>%
    group_by(Indicator, Level) %>%
    summarise(
        n = n()
    ) %>%
    group_by(Indicator) %>%
    mutate(percent = n / sum(n)) %>%
    ungroup() %>%
    mutate(
        Level = recode(Level, !!!level_mapping)
    ) %>%
    mutate_at(
        vars(Level),
        factor,
        levels = rev(level_mapping)
    ) %>%
    left_join(bio_info, by = "Indicator") %>%
    mutate(
        Category = recode(Category, !!!category_mapping)
    ) %>%
    filter(!Category %in% c("Questionnaire"))

fwrite(bio_level_stat, file.path(stat_dir, "bio_level_stat.txt"), sep = "\t")

level_plot <- function(data, indicator_type) {
    ggplot(
        data = data,
        mapping = aes(
            x = percent,
            y = Indicator,
            fill = Level
        )
    ) +
        geom_bar(stat = "identity", position = "stack") +
        scale_x_continuous(expand = c(0, 0)) +
        scale_fill_manual(values = c(
            # 缺乏用蓝色系，过多用红色系，正常用绿色系，色盲友好
            "Too Low" = "#440154",
            "Low" = "#3B528B",
            "Slightly Low" = "#21908C",
            "Normal" = "#5DC863",
            "Slightly High" = "#FDE725",
            "High" = "#F37041",
            "Too High" = "#D73027"
        )) +
        labs(
            x = "Number of sample",
            y = paste0(indicator_type, " Indicator"),
            fill = "Indicator Level"
        ) +
        wdy_theme(
            base_size = 10
        )
}

level_Nutrition_plot <- level_plot(bio_level_stat %>% filter(Category == "Nutrition"), "Nutrition")
level_Pathology_plot <- level_plot(bio_level_stat %>% filter(Category == "Pathology"), "Pathology")
level_Toxin_plot <- level_plot(bio_level_stat %>% filter(Category == "Toxin"), "Toxin")
level_Hormone_plot <- level_plot(bio_level_stat %>% filter(Category == "Hormone"), "Hormone")

ggsave(
    file.path(plot_dir, "level_Nutrition.png"),
    level_Nutrition_plot,
    dpi = 600,
    width = 130, height = 20 + nrow(bio_level_stat %>% filter(Category == "Nutrition")) * 1.5,
    units = "mm"
)
ggsave(
    file.path(plot_dir, "level_Pathology.png"),
    level_Pathology_plot,
    dpi = 600,
    width = 130, height = 20 + nrow(bio_level_stat %>% filter(Category == "Pathology")) * 1.5,
    units = "mm"
)
ggsave(
    file.path(plot_dir, "level_Toxin.png"),
    level_Toxin_plot,
    dpi = 600,
    width = 130, height = 25 + nrow(bio_level_stat %>% filter(Category == "Toxin")) * 1.5,
    units = "mm"
)
ggsave(
    file.path(plot_dir, "level_Hormone.png"),
    level_Hormone_plot,
    dpi = 600,
    width = 130, height = 10 + nrow(bio_level_stat %>% filter(Category == "Hormone")) * 1.5,
    units = "mm"
)

bio_risk <- fread(bio_risk_path, sep = "\t") %>%
    column_to_rownames("Indicator") %>%
    filter(!if_all(.cols = everything(), .fns = is.na)) %>%
    rownames_to_column("Indicator")

bio_risk_stat <- bio_risk %>%
    pivot_longer(cols = -Indicator, names_to = "Sample", values_to = "Risk") %>%
    filter(!is.na(Risk)) %>%
    mutate(
        Risk = recode(Risk, !!!risk_mapping)
    ) %>%
    mutate_at(
        vars(Risk),
        factor,
        levels = rev(c("Normal", "Low Risk", "Medium Risk", "High Risk"))
    ) %>%
    group_by(Indicator, Risk) %>%
    summarise(
        n = n()
    ) %>%
    group_by(Indicator) %>%
    mutate(percent = n / sum(n)) %>%
    ungroup() %>%
    left_join(bio_info, by = "Indicator") %>%
    mutate(
        Category = recode(Category, !!!category_mapping)
    ) %>%
    filter(!Category %in% c("Questionnaire"))

fwrite(bio_risk_stat, file.path(stat_dir, "bio_risk_stat.txt"), sep = "\t")

risk_plot <- function(data, indicator_type) {
    ggplot(
        data = data,
        mapping = aes(
            x = percent,
            y = Indicator,
            fill = Risk
        )
    ) +
        geom_bar(stat = "identity", position = "stack") +
        scale_x_continuous(expand = c(0, 0)) +
        scale_fill_manual(values = c(
            "Normal" = "#5DC863",
            "Low Risk" = "#FDE725",
            "Medium Risk" = "#F37041",
            "High Risk" = "#D73027"
        )) +
        labs(
            x = "Number of sample",
            y = paste0(indicator_type, " Indicator"),
            fill = "Risk Level"
        ) +
        wdy_theme(
            base_size = 10
        )
}

risk_Nutrition_plot <- risk_plot(bio_risk_stat %>% filter(Category == "Nutrition"), "Nutrition")
risk_Pathology_plot <- risk_plot(bio_risk_stat %>% filter(Category == "Pathology"), "Pathology")
risk_Toxin_plot <- risk_plot(bio_risk_stat %>% filter(Category == "Toxin"), "Toxin")
risk_Hormone_plot <- risk_plot(bio_risk_stat %>% filter(Category == "Hormone"), "Hormone")
risk_Imaging_plot <- risk_plot(bio_risk_stat %>% filter(Category == "Imaging"), "Imaging")

ggsave(
    file.path(plot_dir, "risk_Nutrition.png"),
    risk_Nutrition_plot,
    dpi = 600,
    width = 150, height = 20 + nrow(bio_risk_stat %>% filter(Category == "Nutrition")) * 1.5,
    units = "mm"
)
ggsave(
    file.path(plot_dir, "risk_Pathology.png"),
    risk_Pathology_plot,
    dpi = 600,
    width = 150, height = 20 + nrow(bio_risk_stat %>% filter(Category == "Pathology")) * 1.5,
    units = "mm"
)
ggsave(
    file.path(plot_dir, "risk_Toxin.png"),
    risk_Toxin_plot,
    dpi = 600,
    width = 150, height = 25 + nrow(bio_risk_stat %>% filter(Category == "Toxin")) * 1.5,
    units = "mm"
)
ggsave(
    file.path(plot_dir, "risk_Hormone.png"),
    risk_Hormone_plot,
    dpi = 600,
    width = 150, height = 10 + nrow(bio_risk_stat %>% filter(Category == "Hormone")) * 1.5,
    units = "mm"
)
ggsave(
    file.path(plot_dir, "risk_Imaging.png"),
    risk_Imaging_plot,
    dpi = 600,
    width = 150, height = 25 + nrow(bio_risk_stat %>% filter(Category == "Imaging")) * 1.5,
    units = "mm"
)
