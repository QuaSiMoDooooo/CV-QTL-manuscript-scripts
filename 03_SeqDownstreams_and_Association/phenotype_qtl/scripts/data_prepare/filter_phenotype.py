# -*- coding: UTF-8 -*-
#
# FileName     : filter_phenotype.py
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-09-20 21:40
# Last Modified: 2024-11-01 17:30
# Modified By  : EastsunW
# -------------
# Description  : 对基因表达、APA、AS进行质量控制，对于基因表达，去掉低表达基因，对于APA和AS，去掉缺失率过高的样本，表型的ID就是原ID，没有其他信息
# -------------

import click
import pandas as pd


@click.group(invoke_without_command=True)
@click.pass_context
def cmd_group(ctx):
    if not ctx.invoked_subcommand:
        click.echo(ctx.get_help())
    ctx.ensure_object(dict)


@cmd_group.command(
    name='expression',
    help='对基因表达进行质量控制，去掉低表达基因'
)
@click.option(
    '-i', '--input', 'input_file',
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False,
        readable=True,
        resolve_path=True
    ),
    required=True,
    help='输入的基因表达文件'
)
@click.option(
    '--average',
    type=float,
    default=1,
    help='基因平均表达量的阈值，默认为1',
    required=False
)
@click.option(
    '--median',
    type=float,
    default=0,
    help='基因中位表达量的阈值，默认为1',
    required=False
)
@click.option(
    '-o', '--output', 'output_file',
    type=click.Path(
        exists=False, file_okay=True, dir_okay=False,
        writable=True,
        resolve_path=True
    ),
    required=True,
    help='输出的基因表达文件'
)
def filter_expression(input_file, output_file, average, median):
    expression_raw = pd.read_csv(input_file, sep="\t", na_values=["", "NA"])
    expression_expressed = expression_raw.set_index("ID")
    expression_expressed["median"] = expression_expressed.median(axis=1)
    expression_expressed["mean"] = expression_expressed.mean(axis=1)
    expression_expressed = expression_expressed[(
        expression_expressed["median"] > median) & (expression_expressed["mean"] > average)]
    expression_expressed = expression_expressed.drop(
        columns=["median", "mean"])
    expression_expressed = expression_expressed.reset_index()
    expression_expressed.to_csv(
        output_file, sep="\t", index=False, quoting=False)


@cmd_group.command(
    name='APA',
    help='对APA进行质量控制，去掉缺失率过高的样本'
)
@click.option(
    '-i', '--input', 'input_file',
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False,
        readable=True,
        resolve_path=True
    ),
    required=True,
    help='输入的APA文件'
)
@click.option(
    '--na',
    type=float,
    default=0.05,
    help='max missing rate, default is 0.05',
    required=False
)
@click.option(
    '-o', '--output', 'output_file',
    type=click.Path(
        exists=False, file_okay=True, dir_okay=False,
        writable=True,
        resolve_path=True
    ),
    required=True,
    help='输出的APA文件'
)
def filter_APA(input_file, output_file, na):
    apa_raw = pd.read_csv(input_file, sep="\t", na_values=["", "NA"])
    apa_raw = apa_raw[(apa_raw["Gene_Name"] != "NA") &
                      (apa_raw["Gene_Name"].notna())]
    apa_raw = apa_raw.rename(columns={"APA_ID": "ID"})
    apa_raw = apa_raw.loc[:, [
        "ID"] + [col for col in apa_raw.columns if col.endswith(".PAU")]]
    apa_raw.columns = ["ID"] + [col.replace(".PAU", "")
                                for col in apa_raw.columns if col.endswith(".PAU")]
    apa_raw = apa_raw.set_index("ID")
    apa_raw["missing_rate"] = apa_raw.isna().sum(axis=1) / \
        (apa_raw.shape[1] - 1)
    apa_filtered = apa_raw[apa_raw["missing_rate"]
                           < na].drop(columns=["missing_rate"])
    apa_filtered = apa_filtered.fillna(0).reset_index()
    apa_filtered.to_csv(output_file, sep="\t", index=False, quoting=False)


@cmd_group.command(
    name='splicing',
    help='对AS进行质量控制，去掉缺失率过高的样本'
)
@click.option(
    '-i', '--input', 'input_file',
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False,
        readable=True,
        resolve_path=True
    ),
    required=True,
    help='输入的AS文件'
)
@click.option(
    '--na',
    type=float,
    default=0.05,
    help='max missing rate, default is 0.05',
    required=False
)
@click.option(
    '-o', '--output', 'output_file',
    type=click.Path(
        exists=False, file_okay=True, dir_okay=False,
        writable=True,
        resolve_path=True
    ),
    required=True,
    help='输出的AS文件'
)
def filter_splicing(input_file, output_file, na):
    as_raw = pd.read_csv(input_file, sep="\t", na_values=["", "NA"])
    as_raw = as_raw.drop(columns=["#Chr", "start", "end"])
    as_raw = as_raw.set_index("ID")
    as_raw["missing_rate"] = as_raw.isna().sum(axis=1) / (as_raw.shape[1] - 1)
    as_filtered = as_raw[as_raw["missing_rate"]
                         < na].drop(columns=["missing_rate"])
    as_filtered = as_filtered.fillna(0).reset_index()
    as_filtered.to_csv(output_file, sep="\t", index=False, quoting=False)


if __name__ == '__main__':
    cmd_group.add_command(filter_expression)
    cmd_group.add_command(filter_APA)
    cmd_group.add_command(filter_splicing)
    cmd_group()
