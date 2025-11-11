#! python
# -*- coding: UTF-8 -*-
# 
# FileName     : vcf_helper
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-11-22 09:51
# Last Modified: 2024-12-02 21:01
# Modified By  : EastsunW
# -------------
# Description  : 用来解析VCF文件的辅助函数
# -------------

from collections import defaultdict


def vcf_variant_info_parser(info_str: str):
    """解析VCF中的INFO字段，将其解析为字典的格式，对于flag类型的字段，值为True

    Args:
        info_str (str): VCF中的INFO字段

    Returns:
        dict: 解析后的INFO字典，键为INFO字段的key，值为INFO字段的value
    """
    info_dict = defaultdict(str)
    for item in info_str.split(';'):
        if len(item.split('=')) == 2:
            key, value = item.split('=')
            info_dict[key] = value
        elif len(item.split('=')) == 1:
            key = item.split('=')[0]
            info_dict[key] = None
    return info_dict


def vcf_sample_info_parser(field_str, sample_str):
    """解析vcf文件中的sample信息，返回sample信息的列表，列表中的每一项都是一个样本的信息，为一个字典，字典的键为fomat_str中的字段，值为相应字段的值

    Args:
        fomat_str (str): VCF中的format字段的字符串
        sample_strs (List[str]): 样本信息的列表，每个元素为一个样本的信息字符串

    Returns:
        dict: 样本信息，字典的键为样本名，值为样本信息的字典，键为fomat_str中的字段，值为sample_strs中的字段
    """
    format_fields = field_str.split(":")
    format_values = sample_str.split(":")
    assert len(format_fields) == len(format_values), \
        f"Format string {field_str} and sample string {sample_str} do not match"
    # 解析样本信息
    info_dict = {}
    for field_idx, field in enumerate(format_fields):
        info_dict[field] = format_values[field_idx]
    return info_dict
