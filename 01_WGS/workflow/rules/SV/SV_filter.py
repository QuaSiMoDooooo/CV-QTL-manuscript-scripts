# -*- coding: UTF-8 -*-
#
# FileName     : SV_filter
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-06-30 20:57
# Last Modified: 2024-06-30 20:57
# Modified By  : EastsunW
# -------------
# Description  : 对合并后的SV进行质量控制
# -------------

import gzip
from typing import List
from utils import vcf_variant_info_parser, vcf_sample_info_parser, modify_vcf_info


def find_nearst_subsv(sv_format: str, sv_start: int, sv_end: int, sv_genotype_str: str) -> str:
    if not isinstance(sv_start, int):
        sv_start = int(sv_start)
    if not isinstance(sv_end, int):
        sv_end = int(sv_end)
    subsv_genos = sv_genotype_str.split(";")
    distances = {}
    for subsv_geno in subsv_genos:
        subsv_genotype = vcf_sample_info_parser(sv_format, subsv_geno)
        sv_ranges = subsv_genotype["RG"].split("-")
        start2 = int(sv_ranges[0].split("_")[-1])
        end2 = int(sv_ranges[1].split("_")[-1])
        distances[subsv_geno] = abs(sv_start - start2) + abs(sv_end - end2)
    sorted_distances = sorted(distances.items(), key=lambda x: x[1])
    return sorted_distances[0][0]


def sv_merge_software(line: list, keep_order: list = ["cutesv", "svisionpro", "sniffles2"], softwares=["cutesv", "svisionpro", "sniffles2"]):
    order_idx = [softwares.index(i) for i in keep_order]
    genotype_strs = line[9:]
    infos = vcf_variant_info_parser(line[7])
    genotype_infos = [vcf_sample_info_parser(
        line[8], g.split(";")[0]) for g in genotype_strs]
    if genotype_infos[order_idx[0]]["GT"] != "./.":
        # 优先级最高的软件的genotype将被保留
        kept_genotype_str = find_nearst_subsv(
            sv_format=line[8],
            sv_start=line[1],
            sv_end=infos["END"],
            sv_genotype_str=genotype_strs[order_idx[0]]
        )
        # 删除原始info中的SUPPORT和SUPPORT_VEC
        new_info_str = modify_vcf_info(
            vcf_info=line[7],
            key=["SUPPORT", "SUPPORT_VEC", "SOURCE"],
            value=[None, None, softwares[order_idx[0]]]
        )

        # 构建新的line
        new_line = "\t".join(
            line[0:7] + [new_info_str] + [line[8]] + [kept_genotype_str]
        ) + "\n"
        return new_line
    elif genotype_infos[order_idx[1]]["GT"] != "./.":
        kept_genotype_str = find_nearst_subsv(
            sv_format=line[8],
            sv_start=line[1],
            sv_end=infos["END"],
            sv_genotype_str=genotype_strs[order_idx[1]]
        )
        new_info_str = modify_vcf_info(
            vcf_info=line[7],
            key=["SUPPORT", "SUPPORT_VEC", "SOURCE"],
            value=[None, None, softwares[order_idx[1]]]
        )

        # 构建新的line
        new_line = "\t".join(
            line[0:7] + [new_info_str] + [line[8]] + [kept_genotype_str]
        ) + "\n"
        return new_line
    else:
        kept_genotype_str = find_nearst_subsv(
            sv_format=line[8],
            sv_start=line[1],
            sv_end=infos["END"],
            sv_genotype_str=genotype_strs[order_idx[2]]
        )
        new_info_str = modify_vcf_info(
            vcf_info=line[7],
            key=["SUPPORT", "SUPPORT_VEC", "SOURCE"],
            value=[None, None, softwares[order_idx[2]]]
        )

        # 构建新的line
        new_line = "\t".join(
            line[0:7] + [new_info_str] + [line[8]] + [kept_genotype_str]
        ) + "\n"
        return new_line


def sv_filter_software(FILE, OUTPUT, sample_name: str, order: list):
    softwares = []
    for line in FILE:
        if line.startswith("##"):
            OUTPUT.write(line.encode('utf-8'))
            continue
        elif line.startswith("#CHROM"):
            softwares = line.strip().split("\t")[9:]
            OUTPUT.write(("\t".join(
                line.strip().split("\t")[0:9] +
                [sample_name]
            )+"\n").encode('utf-8'))
        else:
            line_splited = line.strip().split("\t")
            sv_infos = vcf_variant_info_parser(line_splited[7])
            if int(sv_infos["SUPPORT"]) > 1:
                OUTPUT.write(sv_merge_software(
                    line=line_splited,
                    keep_order=order,
                    softwares=softwares
                ).encode('utf-8'))
            else:
                continue


def sv_filter_length(FILE, OUTPUT, thresholds: dict, field="AVGLEN"):
    for line in FILE:
        if line.startswith("#"):
            OUTPUT.write(line.encode('utf-8'))
        else:
            line_splited = line.strip().split("\t")
            infos = vcf_variant_info_parser(line_splited[7])
            try:
                match infos["SVTYPE"]:
                    case "INS":
                        if int(infos[field]) <= thresholds["INS"]:
                            OUTPUT.write(line.encode('utf-8'))
                        else:
                            continue
                    case "DEL":
                        if int(infos[field]) <= thresholds["DEL"]:
                            OUTPUT.write(line.encode('utf-8'))
                        else:
                            continue
                    case "INV":
                        if int(infos[field]) <= thresholds["INV"]:
                            OUTPUT.write(line.encode('utf-8'))
                        else:
                            continue
                    case "DUP":
                        if int(infos[field]) <= thresholds["DUP"]:
                            OUTPUT.write(line.encode('utf-8'))
                        else:
                            continue
                    case "BND":
                        OUTPUT.write(line.encode('utf-8'))
            except Exception as e:
                print(f"{line}: {e}")


def sv_filter_depth(FILE, OUTPUT, thresholds: int):
    if not isinstance(thresholds, int):
        thresholds = int(thresholds)
    for line in FILE:
        if line.startswith("#"):
            OUTPUT.write(line.encode('utf-8'))
        else:
            line_splited = line.strip().split("\t")
            genotype_info = vcf_sample_info_parser(
                line_splited[8], line_splited[9])
            if int(genotype_info["DV"]) >= thresholds:
                OUTPUT.write(line.encode('utf-8'))


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
            if not merged:
                merged.append(interval)
            elif merged[-1].chrom != interval.chrom:
                merged.append(interval)
            elif merged[-1].end < interval.start:
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
        if not intervals:
            return False
        elif intervals[0].start > interval.end:
            return False
        elif intervals[-1].end < interval.start:
            return False
        else:
            low, high = 0, len(intervals) - 1
            while low <= high:
                mid = (low + high) // 2
                if intervals[mid].start > interval.end:
                    high = mid - 1
                elif intervals[mid].end < interval.start:
                    low = mid + 1
                else:
                    return True
            return False


def sv_filter_region(FILE, OUTPUT, exclude_regions: List[Bed], exclude_contigs: list):
    for line in FILE:
        if line.startswith("#"):
            OUTPUT.write(line.encode('utf-8'))
        else:
            line_splited = line.strip().split("\t")
            sv_infos = vcf_variant_info_parser(line_splited[7])
            if len(exclude_contigs) > 0 and line_splited[0] in exclude_contigs:
                continue
            else:
                if len(exclude_regions) > 0:
                    interval_check = Interval(
                        chrom=line_splited[0],
                        start=int(line_splited[1]),
                        end=int(sv_infos["END"])
                    )
                    if not any([
                        region_set.check_overlap(interval_check)
                        for region_set in exclude_regions
                    ]):
                        OUTPUT.write(line.encode('utf-8'))
                else:
                    OUTPUT.write(line.encode('utf-8'))


if __name__ == "__main__":
    if snakemake.params.mode:
        mode = snakemake.params.mode
    assert mode in ["software", "length", "depth", "region"], \
        f"Invalid mode: {mode}, must be one of ['software', 'length', 'depth', 'region']"

    if snakemake.input:
        input = snakemake.input[0]
    assert input.endswith(".vcf") or input.endswith(".vcf.gz"), \
        f"Invalid input file: {input}, must be a vcf file"
    VCF = gzip.open(input, "rt") if input.endswith(
        ".vcf.gz") else open(input, "rt")

    if snakemake.output:
        output = snakemake.output[0]
    assert output.endswith(".vcf") or output.endswith(".vcf.gz"), \
        f"Invalid output file: {output}, must be a vcf file"

    with gzip.open(output, "w") as OUT:
        match mode:
            case "software":
                sv_filter_software(
                    FILE=VCF, OUTPUT=OUT,
                    sample_name=snakemake.params.sample_name,
                    order=snakemake.params.software_order.split(",")
                )
            case "length":
                thresholds = {
                    "INS": int(snakemake.params.max_len_ins),
                    "DEL": int(snakemake.params.max_len_del),
                    "INV": int(snakemake.params.max_len_inv),
                    "DUP": int(snakemake.params.max_len_dup),
                }
                sv_filter_length(
                    FILE=VCF, OUTPUT=OUT,
                    thresholds=thresholds,
                    field=snakemake.params.length_field
                )
            case "depth":
                sv_filter_depth(
                    FILE=VCF, OUTPUT=OUT,
                    thresholds=int(snakemake.params.min_depth)
                )
            case "region":
                filter_regions = []
                if isinstance(snakemake.params.centromeres_region, str):
                    filter_regions.append(
                        Bed(snakemake.params.centromeres_region))
                if isinstance(snakemake.params.gap_telomere_region, str):
                    filter_regions.append(
                        Bed(snakemake.params.gap_telomere_region))
                if isinstance(snakemake.input.high_depth_region, str):
                    filter_regions.append(
                        Bed(snakemake.input.high_depth_region))
                if isinstance(snakemake.params.remove_contigs, str):
                    contigs_to_remove = snakemake.params.remove_contigs.split(
                        ",")
                sv_filter_region(
                    FILE=VCF, OUTPUT=OUT,
                    exclude_regions=filter_regions,
                    exclude_contigs=contigs_to_remove
                )
            case _:
                raise ValueError(
                    f"Invalid mode: {mode}, must be one of ['length', 'depth', 'region']")
