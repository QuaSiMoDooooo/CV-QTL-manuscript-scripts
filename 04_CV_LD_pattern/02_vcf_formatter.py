#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# -------------
# FileName     : 02_vcf_formatter.py
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Format VCF files for LD pattern analysis by standardizing variant IDs and genotypes
# -------------


import os

def process_small_variants(infile, outfile, var_type):
    """Process SNP/InDel/MNV variants by standardizing IDs and genotypes"""
    with open(infile, "r") as fr, open(outfile, "w") as fw:
        for line in fr:
            if line.startswith("#"):
                fw.write(line)
                continue

            cols = line.rstrip("\n").split("\t")

            if var_type == "SNP":
                cols[2] = f"{cols[0]}_{cols[1]}_SNP"

            elif var_type in ("InDel", "MNV"):
                suffix = "INDEL" if var_type == "InDel" else "MNV"
                cols[2] = f"{cols[0]}_{cols[1]}_{suffix}"
                cols[3] = "A"
                cols[4] = "T"

            fw.write("\t".join(cols) + "\n")


def process_sv(infile, outfile):
    """Process SV variants by standardizing IDs and keeping only GT information"""
    with open(infile, "r") as fr, open(outfile, "w") as fw:
        for line in fr:
            if line.startswith("#"):
                fw.write(line)
                continue

            cols = line.rstrip("\n").split("\t")
            out = [""] * len(cols)

            out[0] = cols[0]
            out[1] = cols[1]
            sv_type = cols[2].split("_")[1]
            out[2] = f"{cols[0]}_{cols[1]}_{sv_type}"
            out[3] = "A"
            out[4] = "T"
            out[5] = cols[5]
            out[6] = cols[6]
            out[7] = "."
            out[8] = "GT"

            out[9:] = [
                "./." if "." in field[:3] else field[:3]
                for field in cols[9:]
            ]

            fw.write("\t".join(out) + "\n")


def main():
    input_dir  = "01_flt_vcf"
    output_dir = "02_vcf_format_simplified"
    os.makedirs(output_dir, exist_ok=True)

    files = {
        "SNP":   "SNP.biochem_common.vcf",
        "InDel": "InDel.biochem_common.vcf",
        "MNV":   "MNV.biochem_common.vcf",
        "SV":    "SV.biochem_common.vcf",
    }

    for vtype, fname in files.items():
        src = os.path.join(input_dir, fname)
        dst = os.path.join(output_dir, fname)

        print(f"Processing {fname} → {vtype}")

        if vtype == "SV":
            process_sv(src, dst)
        else:
            process_small_variants(src, dst, vtype)

    print(f"\nCompleted! Results saved to {output_dir}")


if __name__ == "__main__":
    main()