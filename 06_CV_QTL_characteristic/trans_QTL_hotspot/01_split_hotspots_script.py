#!/usr/bin/env python3

# -------------
# FileName     : 01_split_hotspots_script.py
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Assign loci to trans-QTL hotspots by grouping SNPs within 1Mb regions
# -------------

import pandas as pd
import numpy as np
import os
import sys

os.chdir("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/08_trans_QTL_hotspot/04_test")

# Receive external parameters
chrm = sys.argv[1]
qtl = sys.argv[2]
qtl_name = sys.argv[3]
data_path = sys.argv[4]

type = "trans"

def assign_locus(chrm, data_path):
    # Read data
    df = pd.read_csv(data_path, sep='\t')
    
    # Sort by count in descending order
    df = df.sort_values(by='count', ascending=False).reset_index(drop=True)
    
    locus_counter = 1
    
    for index, row in df.iterrows():
        # If locus is empty
        if pd.isna(row['locus']):
            # Define current locus name
            current_locus = f'{chrm}_locus{locus_counter}'
            locus_counter += 1
            
            # Get current SNP position and chromosome
            snp_chr = row['snp_chr']
            snp_pos = row['snp_pos']
            
            # Find SNPs within 1Mb range and assign to current locus
            in_range = (df['snp_chr'] == snp_chr) & \
                       (df['snp_pos'] >= snp_pos - 500000) & \
                       (df['snp_pos'] <= snp_pos + 500000) & \
                       (pd.isna(df['locus']))
            
            df.loc[in_range, 'locus'] = current_locus
    
    return df

df_locus = assign_locus(chrm, data_path)

output_dir = f"01_hotspot/{type}-{qtl_name}/output/"
output_path = f"{output_dir}/{type}-{qtl}.{chrm}.hotspot.tsv"

df_locus.to_csv(output_path, sep='\t', index=False)
print(f"Output saved to {output_path}")