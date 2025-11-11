# -*- coding: UTF-8 -*-
#
# FileName     : prepare_gene_position
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-10-22 23:14
# Last Modified: 2024-10-22 23:18
# Modified By  : EastsunW
# -------------
# Description  :
# -------------

import re
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


# region: GENE command
@cmd_group.command(
    name='expression',
    help='从GTF文件中提取基因的位置信息'
)
@click.option(
    '-i', '--input', 'input',
    type=click.Path(exists=True, file_okay=True, dir_okay=False,
                    readable=True, resolve_path=True),
    default=None,
    help='原始GTF文件',
    required=True
)
@click.option(
    '-o', '--output', 'output',
    type=click.Path(exists=False, file_okay=True, dir_okay=False,
                    writable=True, resolve_path=True),
    default="gene.position.txt",
    help='输出的基因位置信息，包含ID、chr、start、end、strand、gene',
    required=True
)
@click.option(
    '--exclude_YM', 'exclude_YM',
    is_flag=True,
    help='是否排除YM染色体',
    default=True,
    required=True
)
@click.pass_context
def expression(ctx, input, output, exclude_YM):
    FILE_IN = prepare_input(input)
    FILE_OUT = open(output, 'wt')
    FILE_OUT.write("ID\tCHR\tSTART\tEND\tSTRAND\tGENE\n")
    for line in FILE_IN:
        if line.startswith("#"):
            continue
        line = line.strip().split("\t")
        gene_id_match = re.match(r'.*gene_id "([^"]+)"', line[8])
        gene_name_match = re.match(r'.*gene_name "([^"]+)"', line[8])
        if line[2] == "gene" and gene_id_match and gene_name_match:
            gene_id = gene_id_match.group(1)
            gene_name = gene_name_match.group(1)
            chr = line[0] if line[0].startswith("chr") else "chr" + line[0]
            if exclude_YM:
                if chr in ["chrY", "chrM"]:
                    continue
            FILE_OUT.write("\t".join([
                gene_id,
                chr,
                line[3],
                line[4],
                line[6],
                gene_name
            ]) + "\n")
# endregion


# region: APA command
@cmd_group.command(
    name='APA',
    help='从QAPA鉴定结果中提取apa事件的位置信息'
)
@click.option(
    '-i', '--input', 'input',
    type=click.Path(exists=True, file_okay=True, dir_okay=False,
                    readable=True, resolve_path=True),
    default=None,
    help='原始鉴定结果文件',
    required=True
)
@click.option(
    '-o', '--output', 'output',
    type=click.Path(exists=False, file_okay=True, dir_okay=False,
                    writable=True, resolve_path=True),
    default="gene.position.txt",
    help='输出的APA位置信息，包含ID、chr、start、end、strand和对应基因名，正链基因的PAS位点是UTR3.End列，负链的PAS位点是UTR3.Start列，如果是区域，则都为UTR3的Start和End',
    required=True
)
@click.option(
    '--exclude_YM', 'exclude_YM',
    is_flag=True,
    help='是否排除YM染色体',
    default=True,
    required=True
)
@click.pass_context
def apa(ctx, input, output, exclude_YM):
    FILE_IN = prepare_input(input)
    FILE_OUT = open(output, 'wt')
    FILE_OUT.write("ID\tCHR\tSTART\tEND\tSTRAND\tGENE\n")
    for line in FILE_IN:
        if line.startswith("APA_ID"):
            continue
        line = line.strip().split("\t")
        chr = line[4] if line[4].startswith("chr") else "chr" + line[4]
        if exclude_YM:
            if chr in ["chrY", "chrM"]:
                continue
        FILE_OUT.write("\t".join([
            line[0],
            chr,
            line[8],
            line[9],
            line[7],
            line[3],
        ]) + "\n")
# endregion


# region: splice command
@cmd_group.command(
    name='splicing',
    help='从splice鉴定结果中提取可变剪切的位置信息'
)
@click.option(
    '-i', '--input', 'input',
    type=click.Path(exists=True, file_okay=True, dir_okay=False,
                    readable=True, resolve_path=True),
    default=None,
    help='原始鉴定结果文件',
    required=True
)
@click.option(
    '-o', '--output', 'output',
    type=click.Path(exists=False, file_okay=True, dir_okay=False,
                    writable=True, resolve_path=True),
    default="gene.position.txt",
    help='输出的splic位置信息，包含ID、chr、start、end、strand，正链基因的PAS位点是UTR3.End列，负链的PAS位点是UTR3.Start列，如果是区域，则都为UTR3的Start和End',
    required=True
)
@click.option(
    '--exclude_YM', 'exclude_YM',
    is_flag=True,
    help='是否排除YM染色体',
    default=True,
    required=True
)
@click.pass_context
def splicing(ctx, input, output, exclude_YM):
    FILE_IN = prepare_input(input)
    FILE_OUT = open(output, 'wt')
    FILE_OUT.write("ID\tCHR\tSTART\tEND\tSTRAND\tGENE\n")
    for line in FILE_IN:
        if line.startswith("#"):
            continue
        line = line.strip().split("\t")
        chr,start,end,_,_,strand = line[3].replace("_", ":").split(":")
        if not chr.startswith("chr"):
            chr = "chr" + chr
        if exclude_YM:
            if chr in ["chrY", "chrM"]:
                continue
        FILE_OUT.write("\t".join([
            line[3],
            chr,
            start,
            end,
            strand,
        ]) + "\n")
# endregion


if __name__ == '__main__':
    # use cmd_group.add_command(<command function name>) to add commands
    cmd_group.add_command(expression)
    cmd_group.add_command(apa)
    cmd_group.add_command(splicing)
    cmd_group()
