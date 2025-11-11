#! python
# -*- coding: UTF-8 -*-
# 
# FileName     : reformat_INS
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2025-01-06 15:46
# Last Modified: 2025-01-06 16:03
# Modified By  : EastsunW
# -------------
# Description  : 
# -------------

def process_indel_bed(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            fields = line.strip().split('\t')
            chrom, start, end, name = fields
            parts = name.split('_')
            ref = parts[2]
            alt = parts[3]
            if len(alt) > len(ref):
                # INS: 将坐标改为插入点的坐标
                new_start = start
                new_end = str(int(start) + 1)
            else:
                new_start = start
                new_end = end
            # 写入新文件
            outfile.write(f"{chrom}\t{new_start}\t{new_end}\t{name}\n")

def process_svins_bed(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            fields = line.strip().split('\t')
            chrom, start, end, name = fields
            prefix,svtype,ID = name.split('_')
            if svtype == 'INS':
                new_start = start
                new_end = str(int(start) + 1)
            else:
                new_start = start
                new_end = end
            # 写入新文件
            outfile.write(f"{chrom}\t{new_start}\t{new_end}\t{name}\n")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 4:
        print("Usage: python process_indel_bed.py <input_bed_file> <output_bed_file>")
        sys.exit(1)

    file_type = sys.argv[1]
    input_file = sys.argv[2]
    output_file = sys.argv[3]
    if file_type == 'indel':
        process_indel_bed(input_file, output_file)
    elif file_type == 'sv':
        process_svins_bed(input_file, output_file)
