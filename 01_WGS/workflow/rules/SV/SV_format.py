# -*- coding: UTF-8 -*-
#
# FileName     : 1_formating
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-06-12 16:27
# Last Modified: 2024-06-12 16:27
# Modified By  : EastsunW
# -------------
# Description  : 1. 将每个SV鉴定软件的原始结果的格式进行统一化，保留的信息包括：INFO部分的SVTYPE、SVLEN、CHR2、END、RE；FORMAT部分的GT、DR、DV; 2. 基因型会根据DR和DV重新定义（0~0.2;0.2~0.8;0.8~1）; 3. 此脚本针对单个样本，不可用于多样本的VCF
# -------------

import gzip
from utils import vcf_variant_info_parser, vcf_sample_info_parser


def genotype_redefine(dr, dv, threshold=0.2):
    """根据支持alt等位基因的reads的比率确定SV的基因型

    Args:
        dr (int): 支持ref等位基因的reads数
        dv (int): 支持alt等位基因的reads数
        threshold (float): SV的基因型的判断阈值，默认为0.2

    Returns:
        str: SV的基因型
    """
    VAB = dv / (dr + dv)
    if VAB <= threshold:
        return '0/0'
    elif threshold < dv/(dr+dv) <= 1-threshold:
        return '0/1'
    else:
        return '1/1'


def convert_vcf_format(input_vcf, output_vcf, software, af_threshold=0.2):
    with gzip.open(input_vcf, "rt") as VCF_IN:
        # if not os.path.exists(Path(output_vcf).parent):
        #     os.makedirs(Path(output_vcf).parent)
        with gzip.open(output_vcf, "wb") as VCF_OUT:
            for line in VCF_IN:
                if line.startswith("#"):
                    VCF_OUT.write(line.encode('utf-8'))
                    continue
                line_splited = line.strip().split("\t")
                # 获得所有的SV信息
                vcf_infos = vcf_variant_info_parser(line_splited[7])
                # 删掉不准确的SV
                if "IMPRECISE" in vcf_infos:
                    continue
                # chr,pos,svid,ref,alt,qual,filter = line_splited[:7]
                output_list = line_splited[:7]
                # 处理BND的格式
                if vcf_infos["SVTYPE"] in ("BND", "TRA"):
                    # BND没有end，而且需要进一步归类到四大SV分类中
                    if '[' in line_splited[4]:
                        chr2 = line_splited[4].split('[')[1].split(':')[0]
                        end = line_splited[4].split('[')[1].split(':')[1]
                        # 只有pbsv报告了同一染色体上的BND
                        # if chr2 == line_splited[0]:
                        #     if line_splited[4].split('[')[-1] in ['A', 'T', 'C', 'G', 'N']:
                        #         svtype = 'INV'
                        #         svlen = vcf_infos['MATEDIST']
                        #     if line_splited[4].split('[')[0] in ['A', 'T', 'C', 'G', 'N']:
                        #         svtype = 'DEL'
                        #         svlen = '-' + vcf_infos['MATEDIST']
                        # else:
                        #     svtype = 'BND'
                        #     svlen = '0'
                        svtype = 'BND'
                        svlen = '0'
                    elif ']' in line_splited[4]:
                        chr2 = line_splited[4].split(']')[1].split(':')[0]
                        end = line_splited[4].split(']')[1].split(':')[1]
                        # if chr2 == line_splited[0]:
                        #     if line_splited[4].split(']')[0] in ['A', 'T', 'C', 'G', 'N']:
                        #         svtype = 'INV'
                        #         svlen = vcf_infos['MATEDIST']
                        #     if line_splited[4].split(']')[-1] in ['A', 'T', 'C', 'G', 'N']:
                        #         svtype = 'DUP'
                        #         svlen = vcf_infos['MATEDIST']
                        # else:
                        #     svtype = 'BND'
                        #     svlen = '0'
                        svtype = 'BND'
                        svlen = '0'
                else:
                    svtype = vcf_infos["SVTYPE"]
                    svlen = vcf_infos["SVLEN"]
                    chr2 = line_splited[0]
                    end = vcf_infos["END"] if vcf_infos["SVTYPE"] != "INS" else str(int(vcf_infos["END"]) + int(vcf_infos["SVLEN"]) - 1)
                sample_infos = vcf_sample_info_parser(line_splited[8], line_splited[9])
                match software:
                    case "sniffles2":
                        supp_reads = vcf_infos["SUPPORT"]
                        info_str = f"SVTYPE={svtype};SVLEN={svlen};CHR2={chr2};END={end};RE={supp_reads}"
                        if sample_infos["DR"] == ".":
                            sample_infos["DR"] = "0"
                        if sample_infos["DV"] == ".":
                            sample_infos["DV"] = "0"
                        if int(sample_infos["DV"]) < 2:
                            continue
                        genotype = genotype_redefine(int(sample_infos["DR"]), int(sample_infos["DV"]), af_threshold)
                        sample_str = f"{genotype}:{sample_infos['DR']}:{sample_infos['DV']}"
                    case "cutesv":
                        supp_reads = vcf_infos["RE"]
                        info_str = f"SVTYPE={svtype};SVLEN={svlen};CHR2={chr2};END={end};RE={supp_reads}"
                        if sample_infos["DR"] == ".":
                            sample_infos["DR"] = "0"
                        if sample_infos["DV"] == ".":
                            sample_infos["DV"] = "0"
                        if int(sample_infos["DV"]) < 2:
                            continue
                        genotype = genotype_redefine(int(sample_infos["DR"]), int(sample_infos["DV"]), af_threshold)
                        sample_str = f"{genotype}:{sample_infos['DR']}:{sample_infos['DV']}"
                    case "pbsv":
                        supp_reads = sample_infos["AD"].split(",")[1]
                        info_str = f"SVTYPE={svtype};SVLEN={svlen};CHR2={output_list[0]};END={end};RE={supp_reads}"
                        sample_infos["DR"], sample_infos["DV"] = sample_infos['AD'].split(",")
                        if sample_infos['AD'].split(',')[0] == ".":
                            sample_infos["DR"] = "0"
                        if sample_infos['AD'].split(',')[1] == ".":
                            sample_infos["DV"] = "0"
                        if int(sample_infos["DV"]) < 2:
                            continue
                        genotype = genotype_redefine(int(sample_infos["DR"]), int(sample_infos["DV"]), af_threshold)
                        sample_str = f"{genotype}:{sample_infos['DR']}:{sample_infos['DV']}"
                    case "svision":
                        supp_reads = vcf_infos["SUPPORT"]
                        info_str = f"SVTYPE={svtype};SVLEN={svlen};CHR2={chr2};END={end};RE={supp_reads}"
                        if sample_infos["DR"] == ".":
                            sample_infos["DR"] = "0"
                        if sample_infos["DV"] == ".":
                            sample_infos["DV"] = "0"
                        if int(sample_infos["DV"]) < 2:
                            continue
                        genotype = genotype_redefine(int(sample_infos["DR"]), int(sample_infos["DV"]), af_threshold)
                        sample_str = f"{genotype}:{sample_infos['DR']}:{sample_infos['DV']}"
                    case "svisionpro":
                        if "DUP" in svtype:
                            svtype = "DUP"
                        supp_reads = vcf_infos["SUPPORT"]
                        info_str = f"SVTYPE={svtype};SVLEN={svlen};CHR2={chr2};END={end};RE={supp_reads}"
                        if not "DR" in sample_infos:
                            sample_infos["DR"] = "."
                        if not "DV" in sample_infos:
                            sample_infos["DV"] = "."
                        genotype = sample_infos["GT"]
                        sample_str = f"{genotype}:{sample_infos['DR']}:{sample_infos['DV']}"
                format_str = "GT:DR:DV"
                output_list.extend([info_str, format_str, sample_str])
                VCF_OUT.write(("\t".join(output_list) + "\n").encode('utf-8'))

convert_vcf_format(
    snakemake.input[0],
    snakemake.output[0],
    snakemake.params.software,
    snakemake.params.af_threshold
)
