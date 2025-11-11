#! python
# -*- coding: UTF-8 -*-
# 
# FileName     : extract_gt_from_vcf copy
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2025-03-14 20:02
# Last Modified: 2025-03-14 20:04
# Modified By  : EastsunW
# -------------
# Description  : 
# -------------

# -*- coding: UTF-8 -*-
#
# FileName     : extract_gt_from_vcf
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-08-27 10:39
# Last Modified: 2024-08-27 10:39
# Modified By  : EastsunW
# -------------
# Description  : 从VCF中提取出基因型信息, 变异的ID的格式为chr:start-end:原ID:ref:alt:
# -------------

from pathlib import Path
import sys
import click
import gzip
from typing import List
parent_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(parent_dir))
from tools.vcf_sample_info_parser import vcf_sample_info_parser
from tools.vcf_variant_info_parser import vcf_variant_info_parser


class Interval:
    def __init__(self, chrom, start, end, name=None):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.name = name

    def __repr__(self):
        return f"{self.chrom}\t{self.start}\t{self.end}"


class Bed:
    def __init__(self, path, merge=True):
        self.__read_bed(path, merge)

    def __read_bed(self, file_path, merge=True):
        if file_path.endswith(".bed.gz"):
            file = gzip.open(file_path, 'rt')
        else:
            file = open(file_path, 'r')
        intervals = []
        for line in file:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            name = parts[3] if len(parts) > 3 else None
            intervals.append(Interval(chrom, start, end, name))
        self.intervals = intervals
        self.intervals.sort(key=lambda x: (x.chrom, x.start))
        file.close()
        if merge:
            self.__merge_continuous()

    def __merge_continuous(self):
        assert self.intervals, "No intervals to merge"
        merged = []
        for interval in self.intervals:
            if not merged or merged[-1].chrom != interval.chrom or merged[-1].end < interval.start:
                merged.append(interval)
            else:
                merged[-1].end = max(merged[-1].end, interval.end)
        merged.sort(key=lambda x: (x.chrom, x.start))
        self.intervals = merged

    def check_overlap(self, interval: Interval):
        intervals_to_check = self.intervals
        intervals = sorted(
            [i for i in intervals_to_check if i.chrom == interval.chrom],
            key=lambda x: x.start
        )
        if not intervals or intervals[-1].start > interval.end:
            return False
        low, high = 0, len(intervals) - 1
        while low <= high:
            mid = (low + high) // 2
            if intervals[mid].start <= interval.end:
                if intervals[mid].end > interval.start:
                    return True
                high = mid - 1
            else:
                low = mid + 1
        return False


def variant_filter_region(FILE, OUTPUT, exclude_regions: List[Bed], exclude_contigs: List[str]):
    for line in FILE:
        if line.startswith("#"):
            OUTPUT.write(line.encode('utf-8'))
        else:
            line_splited = line.strip().split("\t")
            chr, range, id, ref, alt = line_splited[0].split(":")
            start, end = range.split("-")
            # 筛选染色体
            if len(exclude_contigs) > 0 and chr in exclude_contigs:
                continue
            else:
                if len(exclude_regions) > 0:
                    interval_check = Interval(
                        chrom=chr,
                        start=int(start) - 1,
                        end=int(end)
                    )
                    if not any([
                        region_set.check_overlap(interval_check)
                        for region_set in exclude_regions
                    ]):
                        OUTPUT.write(line.encode('utf-8'))
                else:
                    OUTPUT.write(line.encode('utf-8'))


def genotype_judge_snp(sample_info: dict):
    genotypes = sample_info["GT"].split("/")
    num_1 = genotypes.count("1")
    num_0 = genotypes.count("0")
    num_NA = genotypes.count(".")
    if num_1 != 0:
        return str(num_1)
    elif num_0 != 0:
        return "0"
    else:
        return "NA"


def genotype_judge_sv(sample_info: list[dict] | dict):
    if not isinstance(sample_info, list):
        genotypes = sample_info["GT"].split("/")
        num_1 = genotypes.count("1")
        num_0 = genotypes.count("0")
        num_NA = genotypes.count(".")
        if num_1 != 0:
            return str(num_1)
        elif num_0 != 0:
            return "0"
        else:
            return "NA"
    # 有的SV是多个SV合并的结果
    elif len(sample_info) > 1:
        # 复合SV的基因型根据每个子SV的DR和DV重新计算基因型
        DRs = []
        DVs = []
        for sub_sv_dict in sample_info:
            DVs.append(int(sub_sv_dict["DV"]))
            DRs.append(int(sub_sv_dict["DR"]))
        if sum(DVs) == 0:
            return "0"
        elif sum(DVs)/(sum(DVs) + sum(DVs)) < 0.2:
            return "0"
        elif sum(DVs)/(sum(DVs) + sum(DVs)) > 0.8:
            return "2"
        else:
            return "1"


def genotype_judge_mnv(sample_list: list, sample_index_str: str):
    sample_index_list = [int(i) for i in sample_index_str.split(",")]
    assert max(sample_index_list) + 1 <= len(sample_list) * \
        2, f"sample_index_str is out of range: {sample_index_str}"
    # 产生一个列表，都为0，个数等于sample_list的长度
    genotype_list = [0] * len(sample_list)
    for sample_index in sample_index_list:
        genotype_list[int(sample_index)//2] += 1
    return [str(i) for i in genotype_list]


def prepare_input(path):
    if path.endswith(".gz"):
        VCF_IN = gzip.open(path, 'rt')
    else:
        VCF_IN = open(path, 'rt')
    return VCF_IN


@click.group(invoke_without_command=True)
@click.pass_context
def cmd_group(ctx):
    if not ctx.invoked_subcommand:
        click.echo(ctx.get_help())
    ctx.ensure_object(dict)


@cmd_group.command(
    name='SNP',  # 命令的名字
    help='适用于SNP的基因型提取，提供了一些质量控制选项',  # 命令的帮助信息
)
@click.option(
    '-i', '--input', 'input',
    type=click.Path(exists=True, file_okay=True, dir_okay=False,
                    readable=True, resolve_path=True),
    default=None,
    help='待提取基因型的VCF文件，可以压缩',
    required=True
)
@click.option(
    '-o', '--output', 'output',
    type=click.Path(exists=False, file_okay=True, dir_okay=False,
                    writable=True, resolve_path=True),
    default="genotype.txt",
    help='输出的基因型文件',
    required=True
)
@click.option(
    '--hwe',
    type=float,
    default=None,
    help='哈温平衡的筛选阈值，默认不筛选',
)
@click.option(
    '--maf',
    type=float,
    default=None,
    help='MAF的筛选阈值，默认不筛选',
)
@click.option(
    '--na',
    type=float,
    default=None,
    help='缺失率的筛选阈值，默认不筛选',
)
@click.option(
    '--no_multi', 'multi',
    is_flag=True,
    default=True,
    help='是否删除复等位变异',
)
@click.option(
    '--exclude_YM', 'exclude_YM',
    is_flag=True,
    default=True,
    help='是否删除Y和MT染色体的变异',
)
@click.option(
    '--exclude_region', 'exclude_region',
    type=click.Path(exists=True, file_okay=True, dir_okay=False,
                    readable=True, resolve_path=True),
    help='需要排除的基因组区域（bed格式）',
)
@click.pass_context
def snp(ctx, input, output, hwe, maf, na, multi, exclude_YM, exclude_region):
    VCF_IN = prepare_input(input)
    OUTPUT = open(output, 'wt')
    for line in VCF_IN:
        if line.startswith("##"):
            continue
        elif line.startswith("#CHROM"):
            samples = (line[1:].strip().split("\t"))[9:]
            OUTPUT.write("\t".join(["ID"] + samples) + "\n")
        else:
            line_splited = line.strip().split("\t")
            variant_infos = vcf_variant_info_parser(line_splited[7])
            sample_infos = vcf_sample_info_parser(
                line_splited[8], line_splited[9:])
            genotypes = [genotype_judge_snp(sample_info)
                         for sample_info in sample_infos]
            is_multi = ";" in line_splited[2] if multi else False
            if is_multi:
                continue
            if exclude_region:
                exclude_regions = Bed(exclude_region)
                interval_check = Interval(
                    chrom=line_splited[0],
                    start=int(line_splited[1]) - 1,
                    end=int(line_splited[1]) - 1 + len(line_splited[3])
                )
                if exclude_regions.check_overlap(interval_check):
                    continue
            hwe_pass = float(variant_infos["HWE"]
                             ) > hwe if hwe != None else True
            maf_pass = float(variant_infos["MAF"]
                             ) > maf if maf != None else True
            na_pass = genotypes.count(
                "NA")/len(genotypes) < na if na != None else True
            if exclude_YM and line_splited[0] in ["chrY", "chrM"]:
                continue
            if hwe_pass and maf_pass and na_pass and not is_multi:
                OUTPUT.write("\t".join(
                    [
                        ":".join([
                            line_splited[0],
                            "-".join([line_splited[1], line_splited[1]]),
                            line_splited[2],
                            line_splited[3],
                            line_splited[4],
                        ])
                    ] +
                    genotypes
                ) + "\n")


# using click_command snippet to add commands
# region: command_name command
@cmd_group.command(
    name='SV',  # 命令的名字
    help='适用于SV的基因型提取，提供了一些质量控制选项',  # 命令的帮助信息
)
@click.option(
    '-i', '--input',
    type=click.Path(exists=True, file_okay=True, dir_okay=False,
                    readable=True, resolve_path=True),
    default=None,
    help='待提取基因型的VCF文件，可以压缩',
    required=True
)
@click.option(
    '-o', '--output',
    type=click.Path(exists=False, file_okay=True, dir_okay=False,
                    writable=True, resolve_path=True),
    default="genotype.txt",
    help='输出的基因型文件',
    required=True
)
@click.option(
    '--maf',
    type=float,
    default=None,
    help='最小次等位基因频率，默认不筛选',
)
@click.option(
    '--na',
    type=float,
    default=None,
    help='缺失率的筛选阈值，默认不筛选',
)
@click.option(
    '--no_bnd',
    is_flag=True,
    default=True,
    help='是否删除BND类型的变异，默认删除',
)
@click.option(
    '--exclude_YM', 'exclude_YM',
    is_flag=True,
    default=True,
    help='是否删除Y和MT染色体的变异',
)
@click.option(
    '--exclude_region', 'exclude_region',
    type=click.Path(exists=True, file_okay=True, dir_okay=False,
                    readable=True, resolve_path=True),
    help='需要排除的基因组区域（bed格式）',
)
@click.pass_context
def sv(ctx, input, output, maf, na, no_bnd, exclude_YM, exclude_region):
    print(no_bnd)
    VCF_IN = prepare_input(input)
    OUTPUT = open(output, 'wt')
    for line in VCF_IN:
        if line.startswith("##"):
            continue
        elif line.startswith("#CHROM"):
            samples = (line[1:].strip().split("\t"))[9:]
            OUTPUT.write("\t".join(["ID"] + samples) + "\n")
        else:
            line_splited = line.strip().split("\t")
            variant_infos = vcf_variant_info_parser(line_splited[7])
            sample_infos = vcf_sample_info_parser(
                line_splited[8], line_splited[9:])

            if no_bnd and variant_infos["SVTYPE"] == "BND":
                continue
            if exclude_region:
                exclude_regions = Bed(exclude_region)
                interval_check = Interval(
                    chrom=line_splited[0],
                    start=int(line_splited[1]) - 1,
                    end=int(variant_infos["END"])
                )
                if exclude_regions.check_overlap(interval_check):
                    continue
            genotypes = [genotype_judge_sv(sample_info)
                         for sample_info in sample_infos]
            genotypes = ["0" if g == "NA" else g for g in genotypes]
            genotypes_digit = [int(i) for i in filter(
                lambda x: x.isdigit(), genotypes)]
            line_af = sum(genotypes_digit) / (2*len(genotypes_digit))
            line_maf = min([line_af, 1-line_af])

            maf_pass = line_maf > maf if maf != None else True
            na_pass = genotypes.count(
                "NA")/len(genotypes) < na if na != None else True
            if exclude_YM and line_splited[0] in ["chrY", "chrM"]:
                continue
            if maf_pass and na_pass:
                OUTPUT.write("\t".join(
                    [":".join([
                        line_splited[0],
                        "-".join([line_splited[1], variant_infos["END"]]),
                        line_splited[2],
                        ".",
                        ".",
                    ])] +
                    genotypes
                ) + "\n")
# endregion


# using click_command snippet to add commands
# region: command_name command
@cmd_group.command(
    name='MNV',  # 命令的名字
    help='适用于MNV的基因型提取',  # 命令的帮助信息
)
@click.option(
    '-i', '--input',
    type=click.Path(exists=True, file_okay=True, dir_okay=False,
                    readable=True, resolve_path=True),
    default=None,
    help='MNV鉴定结果，一般是txt格式',
    required=True
)
@click.option(
    '-o', '--output',
    type=click.Path(exists=False, file_okay=True, dir_okay=False,
                    writable=True, resolve_path=True),
    default="genotype.txt",
    help='输出的MNV基因型文件',
    required=True
)
@click.option(
    '-v', '--vcf',
    type=click.Path(exists=False, file_okay=True, dir_okay=False,
                    writable=True, resolve_path=True),
    default=None,
    help='鉴定MNV时使用的VCF文件，可以压缩，用来提取样本信息',
    required=True
)
@click.option(
    '--maf',
    type=float,
    default=None,
    help='最小次等位基因频率，默认不筛选',
)
@click.option(
    '--exclude_YM', 'exclude_YM',
    is_flag=True,
    default=True,
    help='是否删除Y和MT染色体的变异',
)
@click.option(
    '--exclude_region', 'exclude_region',
    type=click.Path(exists=True, file_okay=True, dir_okay=False,
                    readable=True, resolve_path=True),
    help='需要排除的基因组区域（bed格式）',
)
@click.pass_context
def mnv(ctx, input, output, vcf, maf, exclude_YM, exclude_region):
    # your codes here
    VCF_IN = prepare_input(vcf)
    sample_names = []
    for line in VCF_IN:
        if line.startswith("#CHROM"):
            sample_names = (line.strip().split("\t"))[9:]
            break
    VCF_IN.close()
    MNV_IN = prepare_input(input)
    OUTPUT = open(output, 'wt')
    OUTPUT.write("\t".join(["ID"] + sample_names) + "\n")
    for line in MNV_IN:
        if not line.startswith("#"):
            line_splited = line.strip().split("\t")
            if line_splited[0].isdigit():
                line_splited[0] = "chr" + line_splited[0]
            genotypes = genotype_judge_mnv(sample_names, line_splited[11])
            line_maf = min(
                [float(line_splited[9]), 1-float(line_splited[9])])
            maf_pass = line_maf > maf if maf != None else True
            if exclude_YM and line_splited[0] in ["chrY", "chrM"]:
                continue
            if exclude_region:
                exclude_regions = Bed(exclude_region)
                interval_check = Interval(
                    chrom=line_splited[0] if line_splited[0].startswith("chr") else "chr" + line_splited[0],
                    start=int(line_splited[1].split(",")[0]) - 1,
                    end=int(line_splited[1].split(",")[-1])
                )
                if exclude_regions.check_overlap(interval_check):
                    continue
            if maf_pass:
                OUTPUT.write("\t".join(
                    [":".join([
                        line_splited[0] if line_splited[0].startswith("chr") else "chr" + line_splited[0],
                        line_splited[1],
                        line_splited[2],
                        line_splited[3],
                        line_splited[4],
                    ])] +
                    genotypes
                ) + "\n")
# endregion


if __name__ == '__main__':
    # use cmd_group.add_command(<command function name>) to add commands
    cmd_group.add_command(snp)
    cmd_group.add_command(sv)
    cmd_group.add_command(mnv)
    cmd_group()
