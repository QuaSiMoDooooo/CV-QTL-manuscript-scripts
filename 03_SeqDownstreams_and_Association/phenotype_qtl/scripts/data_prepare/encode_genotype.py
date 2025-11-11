# -*- coding: UTF-8 -*-
#
# FileName     : encode_genotype
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-09-19 16:59
# Last Modified: 2024-09-19 17:00
# Modified By  : EastsunW
# -------------
# Description  : 将基因型文件中的样本编号转换为新编码，并删除2个重复样本
# -------------

import pandas as pd

# 删除的样本
sample_to_delete = ["wholeblood_161B", "wholeblood_38B"]

# 导入样本信息
sample_info = pd.read_csv(snakemake.input["sample_info"], sep="\t", na_values=["NA", ""])

def rename_sample(names):
    new_names = []
    for oldname in names:
        if oldname.endswith("A"):
            new_name = sample_info.loc[sample_info["blood_1"] == oldname, "ID"].values
        elif oldname.endswith("B"):
            new_name = sample_info.loc[sample_info["blood_2"] == oldname, "ID"].values
        else:
            new_name = ["unknown"]
        new_names.append(new_name[0] if len(new_name) > 0 else "unknown")
    return new_names

# 编码SNP
genotype = pd.read_csv(snakemake.input["genotype"], sep="\t", na_values=["NA", ""])

# 选择ID和以wholeblood开头的列
genotype = genotype[["ID"] + [col for col in genotype.columns if col.startswith("wholeblood")]]

# 删除指定的样本
genotype = genotype.drop(columns=sample_to_delete, errors='ignore')

# 重命名样本
genotype.columns = ["ID"] + rename_sample(genotype.columns[1:])

# 选择ID和样本信息中存在的列
genotype = genotype[["ID"] + [col for col in sample_info["ID"] if col in genotype.columns]]

# 写入文件
genotype.to_csv(snakemake.output[0], sep="\t", na_rep="NA", index=False, quoting=False)
