# -*- coding: UTF-8 -*-
#
# FileName     : stat_variant_sample_increase
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-11-28 17:27
# Last Modified: 2024-11-28 17:27
# Modified By  : EastsunW
# -------------
# Description  : 统计随样本增加的变异位点
# -------------

import gzip
import pandas as pd

def process_vcf(file_path):
    if file_path.endswith(".gz"):
        VCF_IN = gzip.open(file_path, 'rt')
    else:
        VCF_IN = open(file_path)
    all_SNPs = set()
    all_Indels = set()
    for line in VCF_IN:
        if line.startswith('#'):
            continue
        line_splited = line.strip().split('\t')
        if not ";" in line_splited[2]:
            if line_splited[2] == ".":
                ID = "_".join([
                    line_splited[0],
                    line_splited[1],
                    line_splited[3],
                    line_splited[4]
                ])
                if len(line_splited[3]) == len(line_splited[4]):
                    all_SNPs.add(ID)
                else:
                    all_Indels.add(ID)
            else:
                ID = line_splited[2]
                chr, pos, ref, alt = ID.split("_")
                if len(ref) == len(alt):
                    all_SNPs.add(ID)
                else:
                    all_Indels.add(ID)
        else:
            IDs = line_splited[2].split(";")
            for ID in IDs:
                chr, pos, ref, alt = ID.split("_")
                if len(ref) == len(alt):
                    all_SNPs.add(ID)
                else:
                    all_Indels.add(ID)
    return all_SNPs, all_Indels

if __name__ == "__main__":
    all_variant = {
        "SNP": set(),
        "InDel": set()
    }
    result_dict = {
        str(i): {
            "SNP": 0,
            "InDel": 0
        }
        for i in range(1, len(snakemake.input)+1)
    }
    for i in range(len(snakemake.input)):
        n_SNPs, n_Indels = process_vcf(snakemake.input[i])
        all_variant["SNP"].update(n_SNPs)
        all_variant["InDel"].update(n_Indels)
        result_dict[str(i+1)]["SNP"] = len(all_variant["SNP"])
        result_dict[str(i+1)]["InDel"] = len(all_variant["InDel"])
df = pd.DataFrame.from_dict(result_dict, orient="index")
df.index.name = "n_sample"
df.to_csv(snakemake.output[0], sep="\t")
