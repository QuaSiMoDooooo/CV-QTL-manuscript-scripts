#!/usr/bin/env python3

# -------------
# FileName     : get_hotspot_area.py
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Identify trans-QTL hotspot regions by analyzing SNP density across chromosomes
# -------------

import pandas as pd
import numpy as np
import os

os.chdir("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/08_trans_QTL_hotspot/06_hotsplot_enrich/01_get_hotsplot_area")

qtl_dir = "/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/08_trans_QTL_hotspot/04_test/01_hotspot"
qtl_list = ["eQTL","sQTL","apaQTL","meQTL"]
qtl_name_list = ["expression","splicing","APA","methylation"]

outdir = "01_hotspot_area"

count_thres_list = [3,5,36,4]

for i in range(4):
    type = "trans"
    qtl = qtl_list[i]
    qtl_name = qtl_name_list[i]
    count_thres = count_thres_list[i]
    
    # Create empty DataFrame with chr, start, end columns
    hotspot_area = pd.DataFrame(columns=["chr","start","end"])
    
    for chrm in range(1,23):
        chrm = str(chrm)
        
        # Read hotspot data
        df = pd.read_csv(f"{qtl_dir}/{type}-{qtl_name}/output/{type}-{qtl}.{chrm}.hotspot.tsv", sep="\t")
        
        # Filter by count threshold
        df = df[df["count"] >= count_thres]
        
        # Group by locus and count occurrences
        df_locus = df.groupby("locus").size().reset_index(name="count")
        
        # Filter loci with at least 3 occurrences
        df_locus = df_locus[df_locus["count"] >= 3]
        
        # Select loci that meet the criteria
        df1 = df[df["locus"].isin(df_locus["locus"])]
        
        # Sort by locus, chromosome, and position
        df1 = df1.sort_values(by=["locus","snp_chr","snp_pos"])
        
        # Calculate start and end positions with 1Mb extension
        df1["start"] = (df1.groupby("locus")["snp_pos"].transform("min") - 1000000).clip(lower=1)
        df1["end"] = df1.groupby("locus")["snp_pos"].transform("max") + 1000000
        
        # Extract and rename columns
        df2 = df1[["locus","snp_chr","start","end"]].drop_duplicates()
        df2 = df2.rename(columns={"snp_chr":"chr"})
        df2 = df2[["chr","start","end"]]
        
        # Merge with main DataFrame
        hotspot_area = pd.concat([hotspot_area, df2])
    
    # Sort by chromosome and position
    hotspot_area = hotspot_area.sort_values(by=["chr","start","end"])
    
    # Save results
    hotspot_area.to_csv(f"{outdir}/{type}-{qtl}.hotspot_area.tsv", sep="\t", index=False, header=False)