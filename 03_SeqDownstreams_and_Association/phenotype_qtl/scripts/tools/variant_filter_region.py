#! python
# -*- coding: UTF-8 -*-
#
# FileName     : variant_filter_region
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-12-05 11:06
# Last Modified: 2024-12-05 11:13
# Modified By  : EastsunW
# -------------
# Description  :
# -------------

import gzip
from typing import List

from vcf_variant_info_parser import vcf_variant_info_parser


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


def variant_filter_region(FILE, OUTPUT, variant_type: str, exclude_regions: List[Bed], exclude_contigs: List[str]):
    for line in FILE:
        if line.startswith("#"):
            OUTPUT.write(line.encode('utf-8'))
        else:
            line_splited = line.strip().split("\t")
            chr,range,id,ref,alt = line_splited[0].split(":")
            start,end = range.split("-")
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

if __name__ == "__main__":
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
        filter_regions = []
        if isinstance(snakemake.params.centromeres_region, str):
            filter_regions.append(
                Bed(snakemake.params.centromeres_region))
        if isinstance(snakemake.params.gap_telomere_region, str):
            filter_regions.append(
                Bed(snakemake.params.gap_telomere_region))
        if isinstance(snakemake.params.tandem_repeats, str):
            filter_regions.append(
                Bed(snakemake.params.tandem_repeats))
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
