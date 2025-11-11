#!/usr/bin/env python3

# -------------
# FileName     : 01_split_hotsports_loop.py
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Process trans-QTL data to identify hotspot regions by analyzing SNP density across chromosomes
# -------------

import pandas as pd
import numpy as np
import os

os.chdir("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/08_trans_QTL_hotspot/04_test")

qtl_dir = "/home/wangdy/Projects/Weibin/Downstreams/phenotype_qtl/results"
qtl_list = ["eQTL", "sQTL", "apaQTL", "meQTL"]
qtl_name_list = ["expression", "splicing", "APA", "methylation"]

for i in range(4):
    type = "trans"
    qtl = qtl_list[i]
    qtl_name = qtl_name_list[i]

    # Load trans-QTL data
    data_path = f"{qtl_dir}/SNP-{qtl_name}/QTL_results/trans.filtered.txt.gz"
    data = pd.read_csv(data_path, sep="\t", compression="gzip")
    
    # Extract chromosome and position from Variant column
    data['chr'] = data['Variant'].str.split('_').str[0]
    data['snp_pos'] = data['Variant'].str.split('_').str[1]
    
    # Extract molecular phenotype chromosome and position
    data['molphe_chr'] = data['Phenotype_position'].str.split(':').str[0]
    data['molphe_pos'] = data['Phenotype_position'].str.split(':').str[1].str.split('-').str[0]

    # Select relevant columns
    keep_cols = ["Variant", "chr", "snp_pos", "Phenotype", "molphe_chr", "molphe_pos", "P"]
    data = data[keep_cols]
    
    # Rename columns
    data = data.rename(columns={
        "Variant": "snp", 
        "chr": "chr", 
        "snp_pos": "snp_pos", 
        "Phenotype": "molphe", 
        "molphe_chr": "molphe_chr", 
        "molphe_pos": "molphe_pos", 
        "P": "pvalue"
    })

    # Process chromosome names and positions
    data['chr'] = data['chr'].replace('chrX', 'chr23')
    data['snp_pos'] = data['snp_pos'].astype(int)
    
    # Sort data by chromosome and position
    data_sorted = data.sort_values(by=['chr', 'snp_pos'])
    data_sorted.rename(columns={"chr": "snp_chr"}, inplace=True)
    
    # Standardize chromosome names
    data_sorted['snp_chr'] = data_sorted['snp_chr'].str.replace('chr', '')
    data_sorted['snp_chr'] = data_sorted['snp_chr'].replace('X', '23')
    
    # Calculate SNP frequency for hotspot identification
    df = data_sorted.groupby(['snp', 'snp_chr', 'snp_pos']).size().reset_index(name='count')
    df = df.sort_values(by=['snp_chr', 'snp_pos'])
    df.reset_index(drop=True, inplace=True)
    df['locus'] = "NA"

    data_input = df

    # Create output directories
    input_dir = f"01_hotspot/{type}-{qtl_name}/input/"
    out_dir = f"01_hotspot/{type}-{qtl_name}/output/"
    
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(input_dir, exist_ok=True)

    # Save chromosome-specific input files
    for chrm in data_input['snp_chr'].unique():
        data_input_chr = data_input[data_input['snp_chr'] == chrm]
        output_path = f'{input_dir}/{type}-{qtl}.{chrm}.independent.tsv'
        
        if not os.path.exists(output_path):
            data_input_chr.to_csv(output_path, sep='\t', index=False)
            print(f"{chrm} saved")

    # Process each chromosome with hotspot script
    for chrm in data_input['snp_chr'].unique(): 
        log_path = f"{out_dir}{chrm}.log"
        os.system(f"nohup python 01_split_hotspots_script.py {chrm} {qtl} {qtl_name} {input_dir}/{type}-{qtl}.{chrm}.independent.tsv > {log_path} 2>&1 &")