# -*- coding: UTF-8 -*-
#
# FileName     : coverage_calculate_click.py
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-01-14 23:26
# Last Modified: 2024-01-14 23:26
# Modified By  : EastsunW
# -------------
# Description  : 使用Click库重写的脚本
# -------------

import click
import sys


@click.command()
@click.option('--input', 'input_file', required=True, type=click.Path(exists=True),
              help='测序覆盖度的统计信息，由genomeCoverageBed(bedtools genomecov)软件生成的bedgraph')
@click.option('--output', 'output_file', required=True, type=click.Path(),
              help='各染色体和总的测序覆盖度统计表格输出路径')
@click.option('--chr_len', 'chr_len_file', required=True, type=click.Path(exists=True),
              help='各染色体的长度文件，可以使用参考基因组的索引代替，第一列是染色体号，第二列是染色体长度')
def main(input_file, output_file, chr_len_file):
    if len(sys.argv) == 1:
        # 当没有任何命令行参数传递时，显示帮助信息
        ctx = click.get_current_context()
        click.echo(ctx.get_help())
        ctx.exit()
    chr_length_dict = {}
    # 获得各染色体长度
    with open(chr_len_file, 'r') as chr_length_file:
        for line in chr_length_file:
            fields = line.split()
            chr_length_dict[fields[0]] = int(fields[1])
    total_reference = sum(chr_length_dict.values())

    # 获得各染色体测序碱基数
    sequencing_base_dict = {key: 0 for key in chr_length_dict.keys()}
    with open(input_file, 'r') as input_file:
        for line in input_file:
            fields = line.split()
            if fields[0] in sequencing_base_dict.keys():
                sequencing_base_dict[fields[0]
                                     ] += int(fields[2]) - int(fields[1])
    total_sequenced = sum(sequencing_base_dict.values())

    # 计算测序覆盖度
    with open(output_file, 'w') as output_file:
        output_file.write("chr\tsequenced\treference\tcoverage\n")
        for key in chr_length_dict.keys():
            coverage = sequencing_base_dict[key] / chr_length_dict[key]
            output_file.write(
                f"{key}\t{sequencing_base_dict[key]}\t{chr_length_dict[key]}\t{coverage:.6f}\n")
        total_coverage = total_sequenced / total_reference
        output_file.write(
            f"Total\t{total_sequenced}\t{total_reference}\t{total_coverage:.6f}\n")


if __name__ == '__main__':
    # 创建一个命令对象
    cmd = click.Command('main', callback=main, params=main.params)
    # 检查是否有命令行参数
    if len(sys.argv) == 1:
        # 显示帮助信息
        click.echo(cmd.get_help(click.Context(cmd)))
    else:
        # 调用主函数
        cmd()
