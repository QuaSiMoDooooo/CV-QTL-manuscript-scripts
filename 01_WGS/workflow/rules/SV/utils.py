# -*- coding: UTF-8 -*-
#
# FileName     : utils
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-06-23 22:16
# Last Modified: 2024-06-23 22:16
# Modified By  : EastsunW
# -------------
# Description  : 一些帮助函数
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


def modify_vcf_info(vcf_info: str, key, value):
    """
    修改 VCF 信息。

    Args:
        vcf_info (str): VCF 信息字符串。
        key (str or list): 要修改的键，可以是字符串或字符串列表。
        value (str, list, None): 要修改的值，可以是字符串、字符串列表或 None。

    Returns:
        str: 修改后的 VCF 信息字符串。

    Raises:
        AssertionError: 当 key 和 value 的长度不一致时抛出。
        ValueError: 当 value 不是字符串、字符串列表或 None 时抛出。
    """
    info_dict = vcf_variant_info_parser(vcf_info)
    if isinstance(value, list):
        assert isinstance(key, list) and len(key) == len(
            value), "key 和 value 的长度必须相同"
        for k, v in zip(key, value):
            if v is None:
                if k in info_dict:
                    del info_dict[k]
            else:
                info_dict[k] = v
    elif isinstance(value, str):
        assert isinstance(key, str), "当 value 是字符串时，key 必须是字符串"
        if value is None:
            if key in info_dict:
                del info_dict[key]
        else:
            info_dict[key] = value
    elif value is None:
        assert isinstance(key, str) or isinstance(
            key, list), "当 value 是 None 时，key 必须是字符串或字符串列表"
        if isinstance(key, str):
            key = [key]
        for k in key:
            if k in info_dict:
                del info_dict[k]
    else:
        raise ValueError(
            "value 必须是字符串、字符串列表或 None，而 key 必须是字符串或字符串列表且与 value 的长度相同")
    return ";".join([f"{k}={v}" if v else k for k, v in info_dict.items()])
