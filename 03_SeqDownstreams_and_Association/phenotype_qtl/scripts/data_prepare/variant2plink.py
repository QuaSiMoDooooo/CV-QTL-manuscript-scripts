# -------------
# FileName     : variant2plink
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-09-24 09:45
# Last Modified: 2024-09-24 16:41
# Modified By  : EastsunW
# -------------
# Description  :
# -------------

import pandas as pd


def parse_genotype(genotype, ref, alt, allele):
    if "," in ref or "," in alt:
        ref = ref.replace(",", "_")
        alt = alt.replace(",", "_")
    if "<" in ref or ">" in alt:
        ref = ref.replace("<", "").replace(">", "")
        alt = alt.replace("<", "").replace(">", "")
    if allele == 1:
        match genotype:
            case 0:
                return ref
            case 1:
                return ref
            case 2:
                return alt
            case _:
                return "0"
    elif allele == 2:
        match genotype:
            case 0:
                return ref
            case 1:
                return alt
            case 2:
                return alt
            case _:
                return "0"


# 导入变异数据
sample_info = pd.read_csv(str(snakemake.input["sample_info"]), sep="\t")
genotype = pd.read_csv(str(snakemake.input["genotype"]), sep="\t").\
    set_index("ID")
position = pd.read_csv(str(snakemake.input["position"]), sep="\t")

genotype_samples = [i for i in genotype.columns if i.startswith(
    snakemake.params["sample_prefix"])]
sample_info_filtered = sample_info[sample_info[snakemake.params["id_colname"]].
                                   isin(genotype_samples)].\
    set_index(snakemake.params["id_colname"]).\
    reindex(genotype_samples).\
    reset_index()
gender_map = dict([x.split("=") for x in snakemake.params["gender_map"].split(",")])
# 构建ped数据框
df_ped = pd.DataFrame({
    'FID': sample_info_filtered[snakemake.params["id_colname"]],
    'IID': sample_info_filtered[snakemake.params["id_colname"]],
    'PID': 0,
    'MID': 0,
    'SEX': sample_info_filtered[snakemake.params["gender_colname"]].apply(lambda x: gender_map[x] if x in gender_map else -1),
    'PHENO': 2
})
df_ped.set_index("FID", inplace=True)
# 从基因型数据构建ped数据框的基因型部分
columns = {}
for index, row in genotype.iterrows():
    chr, pos, ID, ref, alt = index.split(":")
    if ref == ".":
        ref = "REF"
    if alt == ".":
        alt = "ALT"
    columns[f"{ID}_A1"] = row.apply(parse_genotype, args=(ref, alt, 1))
    columns[f"{ID}_A2"] = row.apply(parse_genotype, args=(ref, alt, 2))
temp_ped = pd.concat([df_ped, pd.DataFrame(columns)], axis=1)
temp_ped.to_csv(
    snakemake.output["ped_file"],
    sep="\t",
    index=True,
    header=False
)

df_map = pd.DataFrame({
    "CHR": position["CHR"],
    "Variant": position["ID"],
    "GD": 0,
    "BPP": position["START"]
})
df_map.to_csv(
    snakemake.output["map_file"],
    sep="\t",
    index=False,
    header=False
)
