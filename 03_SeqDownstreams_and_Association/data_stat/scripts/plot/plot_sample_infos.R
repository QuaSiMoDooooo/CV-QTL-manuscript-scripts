#! Rscript
# -------------
# FileName     : plot_sample_infos
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-12-23 17:08
# Last Modified: 2025-06-06 21:47
# Modified By  : EastsunW
# -------------
# Description  :
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(cnmap))
suppressPackageStartupMessages(library(sf))
suppressPackageStartupMessages(library(proj4))
suppressPackageStartupMessages(library(scatterpie))
suppressPackageStartupMessages(library(patchwork))
custom_colors <- c(
    "#FA7F6F",
    "#8ECFC9",
    "#FFBE7A",
    "#82B0D2",
    "#BEB8DC",
    "#E7DAD2",
    "#999999"
)
source("https://gitee.com/eastsunw/personal_code_notebook/raw/master/plot_tools/wdy_theme.r")
setwd("/home/wangdy/Projects/Weibin/Downstreams/data_stat")



sample_with_wgs <- c(
    "HN-13404e8FO8", "HN-134e409GB7", "HN-1356c6bey0", "HN-1356c93Wl9", "HN-137124bcI6",
    "HN-13727bfRY1", "HN-13798fdIp5", "HN-1383204uA9", "HN-138505eAx6", "HN-138c29akC3",
    "HN-1392bacGb0", "HN-1397e97xW2", "HN-1399116Jn0", "HN-139cb9aBE3", "HN-1400c41Nt1",
    "HN-1404784Ih1", "HN-1405a43dG3", "HN-1408d50AL7", "HN-140adbdRH8", "HN-1412b66Ry9",
    "HN-1419016PY7", "HN-1419569fS6", "HN-141c330Kv2", "HN-141c6e4Af3", "HN-141dc9cCu9",
    "HN-141e05bIo4", "HN-141f675Oo4", "HN-1426ee0XV9", "HN-1428cdbPi9", "HN-142b979hD0",
    "HN-142bd90Bc9", "HN-143ba80Cu5", "HN-143bf20qQ4", "HN-143cd56AU8", "HN-143dd5eqG0",
    "HN-143e947wK1", "HN-1448a15LZ4", "HN-144a245ND1", "HN-144a68aem7", "HN-144afd0mJ3",
    "HN-145396cZD3", "HN-14567a4Ri9", "HN-1457374Wz6", "HN-1458136hr2", "HN-14597d3io3",
    "HN-14690d0PW5", "HN-1473195nq9", "HN-147a63bNI7", "HN-147b740qY1", "HN-148144bBW1",
    "HN-1487842xq3", "HN-148ac21lu7", "HN-148cbf2lb2", "HN-149cb31WL8", "HN-1500504lh3",
    "HN-1504facDP8", "HN-150c0b5pW7", "HN-150f7e9qz0", "HN-1512f80bE4", "HN-151aa02yB0",
    "HN-151b5acCm4", "HN-153af75OB5", "HN-153c945UR1", "HN-15460d6AG6", "HN-154a506mm4",
    "HN-154bd19ym8", "HN-1558daanj1", "HN-15719d6cl4", "HN-15860a2Sy1", "HN-158e55aav7",
    "HN-1595cc6ld2", "HN-1599109dk7", "HN-1606851VS2", "HN-16078dbSy4", "HN-1693943EM7",
    "HN-23560d8VI1", "HN-237175auy8", "HN-2374084HY5", "HN-2376676TQ8", "HN-237a1e5lK9",
    "HN-237d72ckd7", "HN-238254dDq3", "HN-2384412OP9", "HN-2386765aT4", "HN-238bb17MV0",
    "HN-23998ffVr5", "HN-239c5aekY4", "HN-239d84dIA4", "HN-239df08dq7", "HN-239edf7di1",
    "HN-24022d9Xu2", "HN-2402801lC8", "HN-24084f4hk3", "HN-240c0e1Br7", "HN-240e8bfZA2",
    "HN-24111afoC6", "HN-2411cdczf2", "HN-241b94ddB4", "HN-242534fAz7", "HN-2428333Zi6",
    "HN-242e501nj9", "HN-2433fdanK8", "HN-243511dzK1", "HN-243b996wK3", "HN-243e704pS4",
    "HN-24405c9CO7", "HN-2440b00JR2", "HN-244177fep8", "HN-2446e5aoB6", "HN-244a280kA7",
    "HN-244afc1SQ5", "HN-244bf78Tg2", "HN-2454941HG8", "HN-24563e5wo4", "HN-2458bbeUO2",
    "HN-245ee65tl8", "HN-246b1ccTL1", "HN-246e51aaf5", "HN-247bd10ew2", "HN-247dbcajv5",
    "HN-248fdadAw1", "HN-24963a7zC5", "HN-249bacaBC3", "HN-250267bNo7", "HN-250571eNa4",
    "HN-250b061XQ9", "HN-250bf71Ih2", "HN-250cf9ail9", "HN-250db08aG1", "HN-2518833Pe7",
    "HN-2525054Lj1", "HN-252c782OD0", "HN-2534279Fo9", "HN-25366b9LJ7", "HN-2538f69ux0",
    "HN-253b903Rg3", "HN-254b06dsB6", "HN-255f6adro4", "HN-2569578pp2", "HN-2574851qj7",
    "HN-25769edPg3", "HN-2593aa8TO4", "HN-25944f8ki4", "HN-25999derI2", "HN-2607b4aUS2",
    "HN-260f6eaTc9", "HN-261350aYi6", "HN-26406eaBQ3"
)

sample_info <- fread("data/sample_info.txt") %>%
    select(ID, `性别`, `年龄`, `身份证号`) %>%
    setnames(c("id", "gender", "age", "id_num")) %>%
    filter(
        id %in% sample_with_wgs,
        gender %in% c("男", "女"),
        age > 0
    ) %>%
    mutate(
        age = as.integer(age),
        gender = ifelse(gender == "男", "male", "female"),
        region_code = ifelse(str_length(id_num) == 18, as.integer(substr(id_num, 1, 2)), NA)
    )

region_count <- table(sample_info[, c("region_code", "gender")]) %>%
    as.data.frame() %>%
    pivot_wider(
        id_cols = "region_code",
        names_from = "gender",
        values_from = "Freq"
    )

map_region <- getMap(
    code = "100000",
    subRegion = TRUE,
    returnClass = "sf"
) %>%
    filter(adcode != 10) %>%
    mutate(adcode = substr(adcode, 1, 2)) %>%
    mutate(
        lontitude = apply(., 1, function(row) as.numeric(row[["center"]][1])),
        latitude = apply(., 1, function(row) as.numeric(row[["center"]][2]))
    ) %>%
    left_join(region_count, by = c("adcode" = "region_code")) %>%
    mutate_at(vars(c("male", "female")), ~ ifelse(is.na(.), 0, .)) %>%
    mutate(
        total = male + female,
        colored = ifelse(total > 0, "yes", "no")
    )

plot_region <- map_region %>%
    as.data.frame() %>%
    select(lontitude, latitude, male, female) %>%
    mutate(
        x = mapply(
            function(lon, lat) {
                project(c(lon, lat), "+proj=laea +lat_0=40 +lon_0=104")[1]
            },
            lontitude,
            latitude
        ),
        y = mapply(
            function(lon, lat) {
                project(c(lon, lat), "+proj=laea +lat_0=40 +lon_0=104")[2]
            },
            lontitude,
            latitude
        ),
        total = male + female
    )

# 主地图
main_map <- ggplot() +
    geom_sf(
        data = map_region,
        mapping = aes(fill = colored),
        color = "gray50"
    ) +
    geom_scatterpie(
        data = plot_region,
        mapping = aes(
            x = x,
            y = y,
            r = log2(total + 1) * 2e4 + 8e4
        ),
        size = 0,
        color = NA,
        cols = c("male", "female")
    ) +
    geom_text(
        data = plot_region %>% filter(total > 0),
        mapping = aes(
            x = x,
            y = y,
            label = total
        ),
        size = 9 / .pt,
        color = "black"
    ) +
    coord_sf(
        ylim = c(-2387082, 1654989),
        crs = "+proj=laea +lat_0=40 +lon_0=104"
    ) +
    scale_fill_manual(
        breaks = c("male", "female"),
        values = c(
            "yes" = "gray80",
            "no" = "white",
            "male" = custom_colors[4],
            "female" = custom_colors[1]
        )
    ) +
    guides(
        color = "fill",
        override.aes = list(fill = c(custom_colors[4], custom_colors[1]))
    ) +
    labs(fill = NULL) +
    wdy_theme(
        10,
        legend.position = "inside",
        legend.direction = "horizontal",
        legend.position.inside = c(0.5, 1),
        legend.justification = c(0.5, 1),
        legend.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()
    )

submap <- ggplot() +
    geom_sf(
        data = map_region,
        color = "gray50"
    ) +
    coord_sf(
        xlim = c(117131.4, 2115095),
        ylim = c(-4028017, -1877844),
        crs = "+proj=laea +lat_0=40 +lon_0=104"
    ) +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(
            fill = NA,
            color = "grey10",
            linetype = 1,
            linewidth = 1
        ),
        plot.margin = unit(c(0, 0, 0, 0), "mm")
    )

pdf("results/plot/cohort_sample_region.pdf", width = 4, height = 3.5, onefile = TRUE)
main_map + inset_element(submap, left = 0.01, bottom = 0.01, right = 0.25, top = 0.25)
dev.off()


# 画年龄分布
pdf("results/plot/cohort_sample_age.pdf", width = 3.5, height = 3.5, onefile = TRUE)
print(paste0("Min age: ", min(sample_info$age), ", Max age: ", max(sample_info$age)))
age_stat <- sample_info %>%
    mutate(age_group = case_when(
        age < 40 ~ "<40",
        age <= 50 ~ "40~50",
        .default = ">50"
    )) %>%
    group_by(age_group) %>%
    summarise(n = n()) %>%
    mutate(prop = round(n / sum(n) * 100, 1)) %>%
    mutate(age_group = factor(age_group))
# fwrite(age_stat, "results/stat/clinic/age_stat.txt", sep = "\t")
age_stat %>%
    ggplot(aes(x = "", y = prop, fill = age_group)) +
    geom_col(width = 1, color = "black") +
    geom_text(
        aes(label = paste0(age_group, "\n(", n, ", ", prop, "%)")),
        position = position_stack(vjust = 0.5),  # 自动计算标签位置
        color = "#000000",
        size = 10 / .pt
    ) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = custom_colors) +
    wdy_theme(
        legend.position = "none",
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank()
    )
dev.off()
