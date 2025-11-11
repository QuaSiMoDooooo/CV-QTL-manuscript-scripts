# -*- coding: UTF-8 -*-
#
# FileName     : vcf_sample_info_parser
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-06-12 16:43
# Last Modified: 2024-06-12 16:49
# Modified By  : EastsunW
# -------------
# Description  : 解析vcf文件中的sample信息，返回sample信息的字典
# -------------

from typing import List
import json

def vcf_sample_info_parser(field_str:str, sample_strs:List[str]):
    """解析vcf文件中的sample信息，返回sample信息的列表，列表中的每一项都是一个样本的信息，为一个字典，字典的键为fomat_str中的字段，值为相应字段的值

    Args:
        fomat_str (str): VCF中的format字段的字符串
        sample_strs (List[str]): 样本信息的列表，每个元素为一个样本的信息字符串

    Returns:
        dict: 样本信息，字典的键为样本名，值为样本信息的字典，键为fomat_str中的字段，值为sample_strs中的字段
    """
    sample_info_results = []
    format_fields = field_str.split(":")
    for sample_str in sample_strs:
        if len(sample_str.split(";")) == 1:
            format_values = sample_str.split(":")
            assert len(format_fields) == len(format_values), f"Format string {field_str} and sample string {sample_str} do not match"
            # 解析样本信息
            info_dict = {}
            for field_idx, field in enumerate(format_fields):
                info_dict[field] = format_values[field_idx]
            sample_info_results.append(info_dict)
        # 针对SV的情况，SV有可能是多个SV合并的结果
        elif len(sample_str.split(";")) > 1:
            info_dicts = []
            for sub_sv_str in sample_str.split(";"):
                format_values = sub_sv_str.split(":")
                assert len(format_fields) == len(format_values), f"Format string {field_str} and sample string {sample_str} do not match"
                # 解析样本信息
                info_dict = {}
                for field_idx, field in enumerate(format_fields):
                    info_dict[field] = format_values[field_idx]
                info_dicts.append(info_dict)
            sample_info_results.append(info_dicts)
    return sample_info_results

# 单元测试
if __name__ == "__main__":
    field_str = "GT:DR:DV:PL:GQ"
    sample_strs = ["1/1:0:12:115,31,0:30", "0/1:14:10:34,0,72:34"]
    print(f"测试format：{field_str}\n")
    for idx, result in enumerate(vcf_sample_info_parser(field_str, sample_strs)):
        print(f"测试info：{sample_strs[idx]}")
        print(f"测试结果：{json.dumps(result, indent=2)}\n")
    print(vcf_sample_info_parser(field_str, sample_strs))
