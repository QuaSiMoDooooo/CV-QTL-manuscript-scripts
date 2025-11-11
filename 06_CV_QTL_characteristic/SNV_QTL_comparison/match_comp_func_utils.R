#! Rscript

# -------------
# FileName     : match_comp_func_utils.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Utility functions for QTL comparison analysis including allele alignment, scatter plotting, and QTL matching
# -------------

align_alleles <- function(row) {
    complementary_bases <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
    if (is.na(row$oa) || is.na(row$ea)) {
        row$match_type <- "no_allele_info"
        return(row)
    }
    a2_a1 <- paste(row$a2, row$a1, sep = "-")
    oa_ea <- paste(row$oa, row$ea, sep = "-")
    ea_oa <- paste(row$ea, row$oa, sep = "-")
    
    if (complementary_bases[row$a2] == row$oa && complementary_bases[row$a1] == row$ea) {
        row$a2 <- row$oa
        row$a1 <- row$ea
        row$match_type <- "flip"
    }
    else if (a2_a1 == ea_oa) {
        temp <- row$a1
        row$a1 <- row$a2
        row$a2 <- temp
        row$beta.x <- -row$beta.x
        row$match_type <- "switch"
    }
    else if (complementary_bases[row$a2] == row$ea && complementary_bases[row$a1] == row$oa) {
        row$a2 <- row$ea
        row$a1 <- row$oa
        row$beta.x <- -row$beta.x
        row$match_type <- "flip_switch"
    }
    else {
        row$match_type <- "no_match"
    }
    return(row)
}

plot_scatter <- function(data_deal, cor.value, dir.ratio) {
    source1 = data_deal$source.x
    source2 = data_deal$source.y
    p <- ggplot(data_deal, aes(x = beta.x.scaled, y = beta.y.scaled)) +
        geom_point(aes(color = scaled_same_direction), alpha = 0.5) +
        scale_color_manual(values = c("TRUE" = "#FB7F72", "FALSE" = "#82B1D1")) +
        geom_abline(slope = 1, intercept = 0, color = "#FFC55A", linetype = "dashed") +
        geom_hline(yintercept = 0, linetype = "dashed") +
        geom_vline(xintercept = 0, linetype = "dashed") +
        annotate("text", x = -0.8, y = 0.8, label = sprintf("Cor: %.2f\nDir: %.1f%%", cor.value, dir.ratio), hjust = 0, size = 8) +
        labs(x = paste0("Scaled effect size in ", source1), y = paste0("Scaled effect size in ", source2)) +
        theme_minimal() +
        xlim(-1, 1) + ylim(-1, 1) +
        theme(axis.title.x = element_text(size = 22),
            axis.title.y = element_text(size = 22),
            axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20),
            legend.position = "none")
    p_path_prefix = paste0("scatter_",source1[1],"_",source2[1],"_rep")
    ggsave(paste0(p_path_prefix,".png"), plot = p, width = 8, height = 8)
    ggsave(paste0(p_path_prefix,".pdf"), plot = p, width = 8, height = 8)
    save(p,file = paste0(p_path_prefix,".plot.RData"))
    print(paste0("scatter plot saved to ", p_path_prefix))
}

match_qtl <- function(dat1=NULL,dat2=NULL) {
    matched_file_path = paste0("matched_",dat1$source[1],"_",dat2$source[1],".tsv")
    
    complementary_bases <- c(A = "T", T = "A", C = "G", G = "C")
    tmp = nrow(dat1)
    dat1 <- dat1 %>%
            filter(!(a2 == "A" & a1 == "T") & !(a2 == "T" & a1 == "A") &
            !(a2 == "C" & a1 == "G") & !(a2 == "G" & a1 == "C"))
    print(paste0("remove palindrome alleles: ", tmp - nrow(dat1)))
    dat1 <- dat1 %>% filter(nchar(a1) == 1 & nchar(a2) == 1)
    dat2 <- dat2 %>% filter(nchar(oa) == 1 & nchar(ea) == 1)
    
    data_join = inner_join(dat1, dat2, by = c("var_chr_pos", "molphe_chr_pos"))
    
    data_join_wait_match = data_join[!data_join$var_chr_pos_a2_a1 == data_join$var_chr_pos_oa_ea,]
    data_deal_noneed_match_df = data_join[data_join$var_chr_pos_a2_a1 == data_join$var_chr_pos_oa_ea,]
    data_deal_noneed_match_df$match_type = "match"
    print(paste0("alleles totally matched: ", nrow(data_deal_noneed_match_df)))
    print(paste0("alleles wait for match: ", nrow(data_join_wait_match)))
    
    if ( is.null(data_join_wait_match) ) {
        data_deal_wait_match <- mclapply(1:nrow(data_join_wait_match), function(i) {
                                tryCatch({
                                    align_alleles(data_join_wait_match[i,])
                                }, error = function(e) {
                                    print(paste0("Error in row ", i, ": ", e$message))
                                })
                                }, mc.cores = mc)
        data_deal_wait_match_df = do.call(rbind, data_deal_wait_match)
        print("match result:")
        print(table(data_deal_wait_match_df$match_type))
        print("remove no_match alleles")
        data_deal_wait_match_df = data_deal_wait_match_df[data_deal_wait_match_df$match_type != "no_match",]
        data_deal = rbind(data_deal_noneed_match_df, data_deal_wait_match_df)
    } else {
        data_deal = data_deal_noneed_match_df
    }

    fwrite(data_deal, matched_file_path, sep = "\t")
    print(paste0("matched result saved to ", matched_file_path))

    data_deal$beta.x.direct = ifelse(data_deal$beta.x >= 0, "+", "-")
    data_deal$beta.y.direct = ifelse(data_deal$beta.y >= 0, "+", "-")
    data_deal$same_direction = data_deal$beta.x.direct == data_deal$beta.y.direct
    same_dir_ratio = table(data_deal$same_direction)["TRUE"] / dim(data_deal)[1]
    same_dir_ratio_round = round(same_dir_ratio*100, 2)
    print(paste0("same direction ratio: ", same_dir_ratio))
    
    normalize <- function(x) {
        max_pos <- max(x[x > 0], na.rm = TRUE)
        max_neg <- min(x[x < 0], na.rm = TRUE)
        res <- x
        res[x > 0] <- x[x > 0] / max_pos
        res[x < 0] <- -x[x < 0] / max_neg
        return(res)
    }
    data_deal$beta.x.scaled = normalize(data_deal$beta.x)
    data_deal$beta.y.scaled = normalize(data_deal$beta.y)
    data_deal$beta.x.scaled.direct = ifelse(data_deal$beta.x.scaled >= 0, "+", "-")
    data_deal$beta.y.scaled.direct = ifelse(data_deal$beta.y.scaled >= 0, "+", "-")
    data_deal$scaled_same_direction = data_deal$beta.x.scaled.direct == data_deal$beta.y.scaled.direct
    scaled_same_dir_ratio = table(data_deal$scaled_same_direction)["TRUE"] / dim(data_deal)[1]
    scaled_same_dir_ratio_round = round(scaled_same_dir_ratio*100, 2)
    print(paste0("scaled same direction ratio: ", scaled_same_dir_ratio))

    cor_res = cor.test(data_deal$beta.x, data_deal$beta.y, method = "spearman")
    cor_pvalue = cor_res$p.value
    cor_value_round = round(cor_res$estimate, 2)
    scaled_cor_res = cor.test(data_deal$beta.x.scaled, data_deal$beta.y.scaled, method = "spearman")
    scaled_cor_value_round = round(scaled_cor_res$estimate, 2)
    print(paste0("correlation p-value: ", cor_pvalue))
    print(paste0("correlation coefficient: ", cor_value_round))
    print(paste0("scaled correlation coefficient: ", scaled_cor_value_round))

    plot_scatter(data_deal, cor_value_round, same_dir_ratio_round)
}