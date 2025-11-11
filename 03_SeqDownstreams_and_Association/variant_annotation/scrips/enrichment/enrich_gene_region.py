#! python
# -*- coding: UTF-8 -*-
#
# FileName     : enrich_gene_region
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2025-01-05 21:57
# Last Modified: 2025-01-05 22:02
# Modified By  : EastsunW
# -------------
# Description  :
# # 鉴定步骤：
# 1. 检测 SV 与 CDS 的交集.
# 2. 检查 SV 断点是否位于 UTR，标记 UTR 破坏性 SV.
# 3. 检查 SV 断点是否位于启动子，标记启动子破坏性 SV.
# 4. 若两个断点都在同一基因中且不满足前几项标准，标记为内含子破坏性 SV.
# 5. 未与任何基因区域相交的 SV 标记为基因间 SV.
# -------------

import sys
import pybedtools
from pybedtools.featurefuncs import five_prime
import click

def extract_promoter(gtf, upstream=1000, downstream=0, **args):
    genes = gtf.filter(lambda x: x[2] == 'gene')
    promoters = genes.each(
        lambda x: five_prime(
            x,
            upstrea=upstream,
            downstream=downstream,
            **args
        )
    )
    return promoters.merge()

def extract_exon(gtf, file=None):
    if file:
        return gtf.filter(lambda x: x[2] == 'exon').saveas(file)
    else:
        return gtf.filter(lambda x: x[2] == 'exon').merge()

def extract_cds(gtf, file=None):
    if file:
        return gtf.filter(lambda x: x[2] == 'CDS').saveas(file)
    else:
        return gtf.filter(lambda x: x[2] == 'CDS').merge()

def extract_intron(gtf, file=None):
    genes = gtf.filter(lambda x: x[2] == 'gene')
    exons = gtf.filter(lambda x: x[2] == 'exon')
    if file:
        return genes.subtract(exons).merge().saveas(file)
    else:
        return genes.subtract(exons).merge()

def extract_5utr(gtf, file=None):
    utr5 = gtf.filter(lambda x: x[2] == 'five_prime_UTR')
    return utr5.cat(utr3, postmerge=False).merge()

@click.command()
@click.option(
    '-a', '--annotation', 'annotation',
    type = click.Path(
        exists=True,
        dir_okay=False, file_okay=True,
        resolve_path=True
    ),
    default  = "/home/data_admin/Reference/Annotation/gencode_47.basic.gtf",
    help     = '注释文件（GTF文件）',
    required = False
)
@click.option(
    '-f', '--feture', 'feature',
    type = click.Path(
        exists=True,
        dir_okay=False, file_okay=True,
        resolve_path=True
    ),
    default  = None,
    help     = '检查富集的区域（QTL变异）',
    required = False
)
@click.option(
    '-b', '--background', 'background',
    type = click.Path(
        exists=True,
        dir_okay=False, file_okay=True,
        resolve_path=True
    ),
    default  = None,
    help     = '富集的背景区域（所有变异）',
    required = False
)
def main(annotation, feature, background):
    # 读取 GTF 文件和 BED 文件
    gtf_file = 
    # 提取 CDS 区域
    cds = gtf.filter(lambda x: x[2] == 'CDS').saveas()

    # 提取 UTR 区域
    utr5 = gtf.filter(lambda x: x[2] == 'five_prime_UTR').saveas()
    utr3 = gtf.filter(lambda x: x[2] == 'three_prime_UTR').saveas()
    utr = utr5.cat(utr3, postmerge=False)

    # 提取启动子区域（假设启动子为基因转录起始位点前 1kb）
    promoters = gtf.filter(lambda x: x[2] == 'gene').each(lambda x: pybedtools.create_interval_from_list(
        [x[0], int(x[3]) - 1000, x[3], x[4], x[5], x[6], x[7], x[8], x[9]]))

    # 1. 检测 SV 与 CDS 的交集
    cds_intersections = bed.intersect(cds, wa=True, wb=True)
    print("CDS intersections:")
    for interval in cds_intersections:
        print(interval)

    # 2. 检查 SV 断点是否位于 UTR
    utr_intersections = bed.intersect(utr, wa=True, wb=True)
    utr_destructive = utr_intersections.filter(
        lambda x: x not in cds_intersections)
    print("UTR destructive:")
    for interval in utr_destructive:
        print(interval)

    # 3. 检查 SV 断点是否位于启动子
    promoter_intersections = bed.intersect(promoters, wa=True, wb=True)
    promoter_destructive = promoter_intersections.filter(
        lambda x: x not in cds_intersections and x not in utr_destructive)
    print("Promoter destructive:")
    for interval in promoter_destructive:
        print(interval)

    # 4. 检查内含子破坏性 SV
    intronic_destructive = bed.intersect(gtf, wa=True, wb=True).filter(
        lambda x: x not in cds_intersections and x not in utr_destructive and x not in promoter_destructive)
    print("Intronic destructive:")
    for interval in intronic_destructive:
        print(interval)

    # 5. 未与任何基因区域相交的 SV 标记为基因间 SV
    intergenic = bed.intersect(gtf, v=True)
    print("Intergenic:")
    for interval in intergenic:
        print(interval)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python sv_annotation.py <gtf_file> <bed_file>")
        sys.exit(1)

    gtf_file = sys.argv[1]
    bed_file = sys.argv[2]
    main(gtf_file, bed_file)
