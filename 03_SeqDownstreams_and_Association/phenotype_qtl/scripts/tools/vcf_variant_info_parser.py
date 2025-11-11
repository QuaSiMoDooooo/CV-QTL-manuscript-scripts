# -*- coding: UTF-8 -*-
#
# FileName     : vcf_info_parser
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-06-12 16:39
# Last Modified: 2024-06-12 16:39
# Modified By  : EastsunW
# -------------
# Description  : 解析VCF文件中的INFO字段，返回一个字典，对于非key=value类型的字段，value为True
# -------------

from collections import defaultdict

# region: helper functions
def vcf_variant_info_parser(info_str:str):
    info_dict = defaultdict(str)
    for item in info_str.split(';'):
        if len(item.split('=')) == 2:
            key, value = item.split('=')
            info_dict[key] = value
        elif len(item.split('=')) == 1:
            key = item.split('=')[0]
            info_dict[key] = True
    return info_dict
