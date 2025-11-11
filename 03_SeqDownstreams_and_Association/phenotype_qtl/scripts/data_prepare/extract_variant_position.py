# -*- coding: UTF-8 -*-
#
# FileName     : extract_position_from_genotype
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-09-23 10:24
# Last Modified: 2024-09-23 10:24
# Modified By  : EastsunW
# -------------
# Description  : 提取变异的位置，用于QTL，输入筛选后的基因型文件，输出变异的位置，染色体，开始，结束，ID，这里的ID是变异的原ID，由基因型文件的ID拆分而来
# -------------

import click
import gzip


def prepare_input(path):
    if path.endswith(".gz"):
        VCF_IN = gzip.open(path, 'rt')
    else:
        VCF_IN = open(path, 'rt')
    return VCF_IN


@click.group(invoke_without_command=True)
# common options here
@click.pass_context
def cmd_group(ctx):
    if not ctx.invoked_subcommand:
        click.echo(ctx.get_help())
    ctx.ensure_object(dict)
# endregion


# region: snp command
@cmd_group.command(
    name='SNP',
    help='从SNP基因型的文件中提取变异的位置，主要是拆分ID来获得信息',
)
@click.option(
    '-i', '--input', 'input',
    type=click.Path(exists=True, file_okay=True, dir_okay=False,
                    readable=True, resolve_path=True),
    default=None,
    help='输入的SNP基因型文件',
    required=True
)
@click.option(
    '-o', '--output', 'output',
    type=click.Path(exists=False, file_okay=True, dir_okay=False,
                    writable=True, resolve_path=True),
    default="genotype.txt",
    help='输出的snp位置信息',
    required=True
)
@click.pass_context
def snp(ctx, input, output):
    FILE_IN = prepare_input(input)
    FILE_OUT = open(output, 'wt')
    FILE_OUT.write("ID\tCHR\tSTART\tEND\tREF\tALT\n")
    for line in FILE_IN:
        if line.startswith("ID"):
            continue
        line_splited = line.strip().split("\t")
        chr, pos, ID, ref, alt = line_splited[0].split(":")
        if not chr.startswith("chr"):
            chr = "chr" + chr
        start, end = pos.split("-")
        FILE_OUT.write(
            "\t".join([ID, chr, start, end, ref, alt]) + "\n")
# endregion

# region: indel command


@cmd_group.command(
    name='InDel',
    help='从indel基因型的文件中提取变异的位置，主要是拆分ID来获得信息',
)
@click.option(
    '-i', '--input', 'input',
    type=click.Path(exists=True, file_okay=True, dir_okay=False,
                    readable=True, resolve_path=True),
    default=None,
    help='输入的indel基因型文件',
    required=True
)
@click.option(
    '-o', '--output', 'output',
    type=click.Path(exists=False, file_okay=True, dir_okay=False,
                    writable=True, resolve_path=True),
    default="genotype.txt",
    help='输出的indel位置信息',
    required=True
)
@click.pass_context
def indel(ctx, input, output):
    FILE_IN = prepare_input(input)
    FILE_OUT = open(output, 'wt')
    FILE_OUT.write("ID\tCHR\tSTART\tEND\tREF\tALT\n")
    for line in FILE_IN:
        if line.startswith("ID"):
            continue
        line_splited = line.strip().split("\t")
        chr, pos, ID, ref, alt = line_splited[0].split(":")
        if not chr.startswith("chr"):
            chr = "chr" + chr
        start, end = pos.split("-")
        FILE_OUT.write(
            "\t".join([ID, chr, start, end, ref, alt]) + "\n")
# endregion

# region: sv command


@cmd_group.command(
    name='SV',
    help='从SV基因型的文件中提取变异的位置，主要是拆分ID来获得信息',
)
@click.option(
    '-i', '--input', 'input',
    type=click.Path(exists=True, file_okay=True, dir_okay=False,
                    readable=True, resolve_path=True),
    default=None,
    help='输入的SV基因型文件',
    required=True
)
@click.option(
    '-o', '--output', 'output',
    type=click.Path(exists=False, file_okay=True, dir_okay=False,
                    writable=True, resolve_path=True),
    default="genotype.txt",
    help='输出的SV位置信息',
    required=True
)
@click.pass_context
def sv(ctx, input, output):
    FILE_IN = prepare_input(input)
    FILE_OUT = open(output, 'wt')
    FILE_OUT.write("ID\tCHR\tSTART\tEND\tREF\tALT\n")
    for line in FILE_IN:
        if line.startswith("ID"):
            continue
        line_splited = line.strip().split("\t")
        chr, pos, ID, ref, alt = line_splited[0].split(":")
        if not chr.startswith("chr"):
            chr = "chr" + chr
        start, end = pos.split("-")
        FILE_OUT.write("\t".join([
            ID,
            chr,
            start,
            end,
            ref,
            alt
        ]) + "\n")
# endregion


# region: mnv command
@cmd_group.command(
    name='MNV',
    help='从MNV基因型的文件中提取变异的位置，主要是拆分ID来获得信息',
)
@click.option(
    '-i', '--input', 'input',
    type=click.Path(exists=True, file_okay=True, dir_okay=False,
                    readable=True, resolve_path=True),
    default=None,
    help='输入的MNV基因型文件',
    required=True
)
@click.option(
    '-o', '--output', 'output',
    type=click.Path(exists=False, file_okay=True, dir_okay=False,
                    writable=True, resolve_path=True),
    default="genotype.txt",
    help='输出的MNV位置信息',
    required=True
)
@click.pass_context
def mnv(ctx, input, output):
    FILE_IN = prepare_input(input)
    FILE_OUT = open(output, 'wt')
    FILE_OUT.write("ID\tCHR\tSTART\tEND\tREF\tALT\n")
    for line in FILE_IN:
        if line.startswith("ID"):
            continue
        line_splited = line.strip().split("\t")
        chr, pos, ID, ref, alt = line_splited[0].split(":")
        if not chr.startswith("chr"):
            chr = "chr" + chr
        sites = pos.split(",")
        start = sites[0]
        end = sites[-1]
        FILE_OUT.write("\t".join([
            ID,
            chr,
            start,
            end,
            ref,
            alt
        ]) + "\n")


if __name__ == '__main__':
    # use cmd_group.add_command(<command function name>) to add commands
    cmd_group.add_command(snp)
    cmd_group.add_command(indel)
    cmd_group.add_command(sv)
    cmd_group.add_command(mnv)
    cmd_group()
