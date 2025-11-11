# -*- coding: UTF-8 -*-
#
# FileName     : SV_merge_CAST
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-06-23 21:57
# Last Modified: 2024-06-23 21:57
# Modified By  : EastsunW
# -------------
# Description  : 使用CAST算法，将多个VCF中的SV进行合并去冗余，输出一个VCF文件
# -------------

from collections import defaultdict, Counter
import gzip
import re
from typing import List
from utils import vcf_variant_info_parser, vcf_sample_info_parser
import numpy as np
import networkx as nx
import multiprocessing as mp
import time

def chromosome_sort_key(chrom_tuple):
    chrom_order = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
                   '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
                   '21', '22', 'X', 'Y', "M"]
    try:
        return chrom_order.index(chrom_tuple[0])
    except ValueError:
        # 如果染色体编号不在常规列表中，可以放在最后
        return len(chrom_order)

class SvInfo(object):
    """SV信息类, 用于存储SV的信息"""

    def __init__(
        self,
        sample='',
        sv_id='', chr='', chr2='', start=0, end=0,
        sv_type='', svlen=0,
        genotype='', dr=0, dv=0
    ):
        self.sample = sample
        self.sv_id = sv_id
        self.chr = chr
        self.chr2 = chr2
        self.start = start
        self.end = end
        self.svlen = svlen
        self.sv_type = sv_type
        self.genotype = genotype
        self.dr = dr
        self.dv = dv

def find_stab_clique(clique_list_raw, sv_graph):
    """在一组clique中找到比较稳定的clique，如果两个clique的节点之间连接的边不少于两条，就在这两个clique之间添加一条边

    Args:
        clique_list_raw (list[list[SVInfo]]): clique的列表，每个clique是一个列表，列表元素是SV的信息
        sv_graph (nx.Graph): 所有sv之间的连接关系

    Returns:
        _type_: _description_
    """
    clique_list = []
    clique_list.extend(clique_list_raw) # [[SVInfo], [], ...]
    gc = nx.Graph()
    for clique1 in clique_list:
        for clique2 in clique_list:
            if clique1 != clique2:
                # clique1, clique2: [SVInfo, ...]
                num = 0
                for node1 in clique1:
                    for node2 in clique2:
                        # node1, node2: SVInfo
                        if sv_graph.get_edge_data(node1, node2):
                            num += 1
                if num >= 2:
                    gc.add_edge(clique1, clique2)
    clusters = list(nx.connected_components(gc))
    for cluster in clusters:
        # cluster: [clique1:list, clique2, ...]
        combine = set()
        for clique in cluster:
            # clique: [sv1:SvInfo, sv2, ...]
            clique_list.remove(clique)
            # combine: 包含clique中所有不重复的SV的集合
            combine |= set(clique)
        # clique_list: [tuple[SVInfo], ...]
        clique_list.append(tuple(combine))
    clique_list = sorted(clique_list, key=lambda x: x[0].sample + x[0].sv_id)
    return clique_list

def read_vcf(files, altnames=None):
    """用来将VCF文件中的SV信息读取的函数，并且根据支持alt等位基因的reads的比率确定SV的基因型

    Args:
        files (list): 存储多个VCF文件的列表，列表元素是vcf文件的路径
        af (float): SV的基因型的判断阈值，默认为0.2
        alt_name (list): 假定files中的每一个vcf都只包含一个样本，这个name用来给每一个文件中的样本指定名字，如果不指定，则将文件名作为样本名

    Returns:
        tuple: 返回染色体的长度的set和SV的读取结果，该结果是一个字典，字典的键是SV的类型，值是一个字典，该字典的键是染色体编号，值是包含SvInfo对象的排序后的列表：
        {
            INS: {
                chr1: [SV1, SV2, ...],
            },
            DEL: {
                chr1: [SV1, SV2, ...],
            },
            ...
        }
    """
    SV_informations = {
        'INS': defaultdict(list),
        'DEL': defaultdict(list),
        'INV': defaultdict(list),
        'DUP': defaultdict(list),
        'BND': defaultdict(list)
    }
    contig_infos = set()
    samples = []
    contigs_regex = re.compile(r'##contig=<ID=chr(?P<contig>.*),length=(?P<length>\d+)>')
    # 读取所有的SV
    for vcf_idx, vcf_path in enumerate(files):
        if vcf_path.endswith('.vcf.gz'):
            vcf = gzip.open(vcf_path, 'rt')
        elif vcf_path.endswith('.vcf'):
            vcf = open(vcf_path, 'rt')
        for vcf_line in vcf:
            if vcf_line.startswith('##contig'):
                contig_infos.add(contigs_regex.match(vcf_line).groups())
            elif vcf_line.startswith('#CHROM'):
                sample = vcf_line.strip().split('\t')[9] if altnames is None else altnames[vcf_idx]
                samples.append(sample)
            elif not vcf_line.startswith('#'):
                line_splited = vcf_line.strip().split('\t')
                sv_infos = vcf_variant_info_parser(line_splited[7])
                sample_infos = vcf_sample_info_parser(line_splited[8], line_splited[9])
                if sample_infos['DR'] == '0' and sample_infos['DV'] == '0':
                    continue
                if 'SVLEN' in sv_infos:
                    svlen_filed_name = 'SVLEN'
                elif 'AVGLEN' in sv_infos:
                    svlen_filed_name = 'AVGLEN'
                else:
                    raise ValueError('SVLEN or AVGLEN not found in INFO field')
                temp_sv = SvInfo(
                    sample   = sample,
                    sv_id    = line_splited[2],
                    chr      = line_splited[0],
                    chr2     = sv_infos['CHR2'],
                    start    = int(line_splited[1]),
                    end      = int(sv_infos['END']),
                    svlen    = int(sv_infos[svlen_filed_name]) ,
                    sv_type  = sv_infos['SVTYPE'],
                    genotype = sample_infos['GT'],
                    dr       = sample_infos['DR'],
                    dv       = sample_infos['DV']
                )
                if temp_sv.sv_type == 'BND':
                    SV_informations['BND'][(temp_sv.chr, temp_sv.chr2)].append(temp_sv)
                else:
                    SV_informations[temp_sv.sv_type][temp_sv.chr].append(temp_sv)
    # 排序每个类型每个染色体的SV
    for sv_type in SV_informations:
        if sv_type == 'BND':
            continue
        else:
            for chrom in SV_informations[sv_type]:
                start = np.array(
                    [i.start for i in SV_informations[sv_type][chrom]])
                # end = start + svlen
                end = np.array([i.start + abs(i.svlen)
                               for i in SV_informations[sv_type][chrom]])
                order = np.lexsort((end, start))
                per_chrom = np.array(SV_informations[sv_type][chrom])[order]
                SV_informations[sv_type][chrom] = list(per_chrom)
    return sorted(contig_infos, key=chromosome_sort_key), samples, SV_informations

def group_by_block(all_sv_form, interval=0):
    """根据所有SV信息，将每条染色体中的所有SV划分为不重叠的block，返回每个block中的SV信息

    Args:
        all_sv_form (dict): 由read_vcf函数返回的SV信息字典
        interval (int, optional): block之间的最小距离. 默认为0.

    Returns:
        dict: 格式与all_sv_form相同，区别在于每个SV类型的每个染色体的值是一个列表的列表，每个子列表代表一个block，其中是block中的SV
    """
    all_sv_block = {
        'INS': defaultdict(list),
        'DEL': defaultdict(list),
        'INV': defaultdict(list),
        'DUP': defaultdict(list),
        'BND': defaultdict(list)
    }

    for type_sv in all_sv_form:
        # BND类型的SV不需要分块
        if type_sv == 'BND':
            for chrom in all_sv_form[type_sv]:
                all_sv_block[type_sv][chrom].append(all_sv_form[type_sv][chrom])
        else:
            for chrom in all_sv_form[type_sv]:
                starts = np.array([i.start for i in all_sv_form[type_sv][chrom]])
                ends = np.array([i.start + abs(i.svlen) for i in all_sv_form[type_sv][chrom]])
                break_idx = []
                if len(starts) != 1:
                    for s, e in zip(
                        range(1, len(starts)),
                        range(len(ends) - 1)
                    ):
                        if starts[s] - ends[e] > interval:
                            break_idx.append(s)
                else:
                    break_idx = []
                    # print('%s chromosome have only one %s!' % (chrom, type_sv))
                block = []
                if len(break_idx) != 0:
                    block.append(all_sv_form[type_sv][chrom][0:break_idx[0]])
                    if len(break_idx) > 1:
                        for b in range(len(break_idx) - 1):
                            block.append(all_sv_form[type_sv][chrom][break_idx[b]:break_idx[b + 1]])
                    block.append(all_sv_form[type_sv][chrom][break_idx[-1]:len(all_sv_form[type_sv][chrom])])
                else:
                    block = [all_sv_form[type_sv][chrom]]
                all_sv_block[type_sv][chrom] = block
    return all_sv_block

def overlap_length(start1, end1, start2, end2):
    return max(max((end2 - start1), 0) - max((end2 - end1), 0) - max((start2 - start1), 0), 0)

def sort_function_clique_max_list(per):
    per_list = []
    for i in per[0]:
        per_list.append(i.sample + i.sv_id)
    return tuple(per_list)

def block_merge_process(all_sv_block, type_sv, chrom, overlap_perc):
    per_chrom_merge = []
    for per_block in all_sv_block[type_sv][chrom]:
        g_non_tra = nx.Graph()
        for per_sv in per_block:
            g_non_tra.add_node(per_sv)
            for other_sv in per_block:
                start1 = per_sv.start
                end1 = per_sv.start + abs(per_sv.svlen)  # end = start + svlen
                start2 = other_sv.start
                end2 = other_sv.start + abs(other_sv.svlen)  # end = start + svlen
                overlap = overlap_length(start1, end1, start2, end2)
                if overlap != 0:
                    rate1 = overlap / (end1 - start1)
                    rate2 = overlap / (end2 - start2)
                    if rate1 > overlap_perc and rate2 > overlap_perc:
                        g_non_tra.add_edge(per_sv, other_sv, weight=rate1 + rate2)
        # to find the largest clique
        max_clique_list = []
        temp_nodes = set(g_non_tra.nodes)
        while temp_nodes:
            temp_graph = nx.Graph.subgraph(g_non_tra, list(temp_nodes))
            clique_list = [(i, len(i)) for i in nx.algorithms.clique.find_cliques(temp_graph)]
            clique_list = sorted(clique_list, key=lambda x: x[1], reverse=True)
            # take largest clique
            clique_max_num = clique_list[0][1]
            clique_max_list = [i for i in clique_list if i[1] == clique_max_num]
            # sort the elements in the group
            clique_max_list_2 = []
            for per in clique_max_list:
                per = list(per)
                per[0] = sorted(per[0], key=lambda x: x.sample + x.sv_id)
                per = tuple(per)
                clique_max_list_2.append(per)
            # sort group
            clique_max_list = sorted(clique_max_list_2, key=sort_function_clique_max_list)
            # take the first largest
            clique_max = clique_max_list[0][0]
            # 这里append后，max_clique_list中的元素会变成SVInfo对象的元组，不知道是为什么
            max_clique_list.append(tuple(clique_max))
            temp_nodes -= set(clique_max)

        # merge clique
        final_merge = []
        if len(max_clique_list) > 1:
            last_clique_list, stab_clique_list = [], []
            stab_clique_list.extend(max_clique_list)
            while last_clique_list != stab_clique_list:
                last_clique_list = []
                last_clique_list.extend(stab_clique_list)
                stab_clique_list = []
                stab_clique_list.extend(find_stab_clique(last_clique_list, g_non_tra))
            # a group containing only one element
            clique_list_one = [i[0] for i in stab_clique_list if len(i) == 1]
            clique_list_one = sorted(clique_list_one, key=lambda x: x.sample + x.sv_id)
            # the elements in the group are greater than two
            clique_list_two = [i for i in stab_clique_list if len(i) >= 2]
            # sort the elements in the group
            clique_list_two_2 = []
            for per in clique_list_two:
                per = list(per)
                per = sorted(per, key=lambda x: x.sample + x.sv_id)
                per = tuple(per)
                clique_list_two_2.append(per)
            # sort group
            clique_list_two_2 = sorted(clique_list_two_2, key=lambda x: x[0].sample + x[0].sv_id)
            # add the SV to the nearest group
            for one_sv in clique_list_one:
                max_weight, max_index = 0, 'none'
                for idx_cl, per_cli in enumerate(clique_list_two_2):
                    num_weight = 0
                    for sv in per_cli:
                        if g_non_tra.get_edge_data(sv, one_sv):
                            num_weight += g_non_tra.get_edge_data(sv, one_sv)['weight']
                    if num_weight > max_weight:
                        max_index, max_weight = idx_cl, num_weight
                if max_index == 'none':
                    clique_list_two_2.append(tuple([one_sv]))
                else:
                    new_per_cli = list(clique_list_two_2[max_index])
                    new_per_cli.append(one_sv)
                    clique_list_two_2[max_index] = tuple(new_per_cli)
            final_merge.extend(clique_list_two_2)
            per_chrom_merge.append(final_merge)
        else:
            final_merge.append(max_clique_list[0])
            per_chrom_merge.append(final_merge)
    return per_chrom_merge, chrom

def block_merge(sv_blocks, extended_dis, overlap_perc, worker=1):
    all_sv_merge = {
        'INS': defaultdict(list),
        'DEL': defaultdict(list),
        'INV': defaultdict(list),
        'DUP': defaultdict(list),
        'BND': defaultdict(list)
    }
    for chrom in sv_blocks['BND']:
        g = nx.Graph()
        # 每个SV作为一个节点，检查每两个SV之间的overlap长度，如果overlap长度大于0，则添加一条边，边的权重为首尾overlap长度之和
        for sv_1 in sv_blocks['BND'][chrom][0]:
            # BND类型的SV每一种断点关系中只有一个block，例如chr1到chr2的所有BND都只作为一个block
            g.add_node(sv_1)
            for sv_2 in sv_blocks['BND'][chrom][0]:
                if sv_1 == sv_2:
                    continue
                sv1_start_min, sv1_start_max = sv_1.start - extended_dis, sv_1.start + extended_dis
                sv1_end_min, sv1_end_max = sv_1.end - extended_dis, sv_1.end + extended_dis
                sv2_start_min, sv2_start_max = sv_2.start - extended_dis, sv_2.start + extended_dis
                sv2_end_min, sv2_end_max = sv_2.end - extended_dis, sv_2.end + extended_dis
                overlap_len_start = overlap_length(sv1_start_min, sv1_start_max, sv2_start_min, sv2_start_max)
                overlap_len_end = overlap_length(sv1_end_min, sv1_end_max, sv2_end_min, sv2_end_max)
                if overlap_len_start > 0 and overlap_len_end > 0:
                    g.add_edge(
                        sv_1, sv_2,
                        weight=overlap_len_start + overlap_len_end
                    )
        # 从网络中找到最大的子图
        max_clique_list_BND = [] # [tuple(list, num), ...]
        temp_nodes_BND = set(g.nodes)
        while temp_nodes_BND:
            # 定义现有整个网络为一个子图
            temp_graph_BND = nx.Graph.subgraph(g, list(temp_nodes_BND))
            # 找出这个图中的所有子图
            clique_list_BND = [(i, len(i)) for i in nx.algorithms.clique.find_cliques(temp_graph_BND)]
            # 按照子图的节点数量降序排列
            clique_list_BND = sorted(clique_list_BND, key=lambda x: x[1], reverse=True)
            # 找出节点数量等于最大子图的所有子图
            clique_max_BND_num = clique_list_BND[0][1]
            clique_max_BND_list = [i for i in clique_list_BND if i[1] == clique_max_BND_num]
            # 将clique_max_BND_list中的每个clique中的节点按照sample和sv_id排序
            clique_max_BND_list_2 = []
            for per in clique_max_BND_list:
                # per: (clique, len(clique))
                # clique: [sv1:SvInfo, sv2, ...]
                per = list(per)
                per[0] = sorted(per[0], key=lambda x: x.sample + x.sv_id)
                per = tuple(per)
                clique_max_BND_list_2.append(per)
            clique_max_BND_list = sorted(clique_max_BND_list_2, key=sort_function_clique_max_list)
            # 将排在第一位的clique作为最大clique提取出来，继续循环进行下一轮寻找
            # clique_max_BND: [sv1:SvInfo, sv2, ...]
            clique_max_BND = clique_max_BND_list[0][0]
            # max_clique_list_BND: [clique:list[SVInfo]]
            max_clique_list_BND.append(tuple(clique_max_BND))
            temp_nodes_BND -= set(clique_max_BND)

        # 合并每一个clique中的SV
        final_merge_BND = []
        if len(max_clique_list_BND) > 1:
            last_clique_list_BND, stab_clique_list_BND = [], []
            stab_clique_list_BND.extend(max_clique_list_BND)
            # 一个迭代，直到最后的clique_list_BND不再变化
            while last_clique_list_BND != stab_clique_list_BND:
                last_clique_list_BND = []
                last_clique_list_BND.extend(stab_clique_list_BND)
                stab_clique_list_BND = []
                stab_clique_list_BND.extend(find_stab_clique(last_clique_list_BND, g))
            # a group containing only one element
            clique_list_one = [i[0] for i in stab_clique_list_BND if len(i) == 1]
            clique_list_one = sorted(clique_list_one, key=lambda x: x.sample + x.sv_id)
            # the elements in the group are greater than two
            clique_list_two = [i for i in stab_clique_list_BND if len(i) >= 2]
            # sort the elements in the group
            clique_list_two_2 = []
            for per in clique_list_two:
                per = list(per)
                per = sorted(per, key=lambda x: x.sample + x.sv_id)
                per = tuple(per)
                clique_list_two_2.append(per)
            # sort group
            clique_list_two_2 = sorted(clique_list_two_2, key=lambda x: x[0].sample + x[0].sv_id)
            # add the SV to the nearest group
            for one_sv in clique_list_one:
                max_weight, max_index = 0, 'none'
                for idx_cl, per_cli in enumerate(clique_list_two_2):
                    num_weight = 0
                    for sv in per_cli:
                        if g.get_edge_data(sv, one_sv):
                            num_weight += g.get_edge_data(sv, one_sv)['weight']
                    if num_weight > max_weight:
                        max_index, max_weight = idx_cl, num_weight
                if max_index == 'none':
                    clique_list_two_2.append(tuple([one_sv]))
                else:
                    new_per_cli = list(clique_list_two_2[max_index])
                    new_per_cli.append(one_sv)
                    clique_list_two_2[max_index] = tuple(new_per_cli)
            final_merge_BND.extend(clique_list_two_2)
            all_sv_merge['BND'][chrom].append(final_merge_BND)
        else:
            final_merge_BND.append(max_clique_list_BND[0])
            all_sv_merge['BND'][chrom].append(final_merge_BND)
    for type_sv in sv_blocks:
        if type_sv != 'BND':
            # print('Start %s: block2clique2merge' % type_sv)
            process_pool = mp.Pool(worker)
            multi_result = [
                process_pool.apply_async(
                    func=block_merge_process,
                    args=(sv_blocks, type_sv, chrom, overlap_perc)
                )
                for chrom in sv_blocks[type_sv]
            ]
            process_pool.close()
            process_pool.join()
            for result in multi_result:
                merged_result, chrom = result.get()
                all_sv_merge[type_sv][chrom] = merged_result
            # for chrom in sv_blocks[type_sv]:
            #     merged_result, chrom = block_merge_process(sv_blocks, type_sv, chrom, overlap_perc)
            #     all_sv_merge[type_sv][chrom] = merged_result
    return all_sv_merge

def vcf_header(contig_lens:List[tuple], samples:List[str]):
    help_info_base = [
        '##fileformat=VCFv4.2',
        '##source=merged',
        f'##fileDate="{time.strftime("%Y/%m/%d %H:%M", time.localtime())}"'
    ]
    help_info_contigs = [f'##contig=<ID=chr{contig},length={length}>' for contig, length in contig_lens]
    help_info_svtype = [
        '##ALT=<ID=INS,Description="Insertion">',
        '##ALT=<ID=DEL,Description="Deletion">',
        '##ALT=<ID=DUP,Description="Duplication">',
        '##ALT=<ID=INV,Description="Inversion">',
        '##ALT=<ID=BND,Description="Breakend; Translocation">'
    ]
    help_info_info = [
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV">',
        '##INFO=<ID=AVGLEN,Number=1,Type=Integer,Description="Average length of merged structural variation">',
        '##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for breakend">',
        '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Interval around POS for merged SV">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the SV">',
        '##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Interval around END for merged SV">',
        '##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of samples supporting the structural variation">',
        '##INFO=<ID=SUPPORT_VEC,Number=1,Type=String,Description="Indexes of samples supporting the structural variation">'
    ]
    help_info_format = [
        # GT:LN:DR:DV:RG
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=LN,Number=1,Type=Integer,Description="SV length">',
        '##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">',
        '##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">',
        '##FORMAT=<ID=RG,Number=1,Type=String,Description="SV range represented as CHR_POS-CHR2_END">',
    ]
    header = "\t".join(
        [
        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"
        ] + samples
    )
    return "\n".join(help_info_base + help_info_contigs + help_info_svtype + help_info_info + help_info_format + [header]) + "\n"

def output_vcf(output, header, all_sv_merge, samples):
    total_wait_dict = defaultdict(list)
    for sv_type in all_sv_merge:
        for chrom in all_sv_merge[sv_type]:
            for sv_block in all_sv_merge[sv_type][chrom]:
                for clique in sv_block:
                    clique_svs = sorted(list(clique), key=lambda x: x.sample + x.sv_id)
                    clique_start = [sv.start for sv in clique_svs]
                    clique_end = [sv.end for sv in clique_svs]
                    # 合并后的SV的start和end取聚类中出现最多的SV的start和end，并给出cipos和ciend
                    most_common_start_idx = clique_start.index(Counter(clique_start).most_common(1)[0][0])
                    most_common_sv = clique_svs[most_common_start_idx]
                    start = most_common_sv.start
                    end = most_common_sv.end
                    chr = most_common_sv.chr
                    chr2 = most_common_sv.chr2
                    cipos = [0, 0]
                    max_start = int(np.max(np.array(clique_start)))
                    min_start = int(np.min(np.array(clique_start)))
                    cipos[0] = min_start - start
                    cipos[1] = max_start - start
                    cipos = ','.join([str(i) for i in cipos])
                    ciend = [0, 0]
                    max_end = int(np.max(np.array(clique_end)))
                    min_end = int(np.min(np.array(clique_end)))
                    ciend[0] = min_end - end
                    ciend[1] = max_end - end
                    ciend = ','.join([str(i) for i in ciend])
                    if sv_type != 'BND':
                        mean_len = int(np.array([abs(i.svlen) for i in clique_svs]).mean())
                        if sv_type == 'DEL':
                            mean_len = -mean_len
                    else:
                        mean_len = 0
                    # 产生每个clique的支持样本信息
                    supp_vec = ['0'] * len(samples)
                    # 先根据clique中的所有sv，生成这个clique的样本信息
                    # {sample: [sv1_info, sv2_info, ...]}
                    clique_formatstr_allsv = defaultdict(list)
                    for sv in clique_svs:
                        if sv.sample in samples:
                            sample_idx = samples.index(sv.sample)
                            supp_vec[sample_idx] = '1'
                            genotype = sv.genotype
                            length = abs(sv.svlen)
                            dr = sv.dr
                            dv = sv.dv
                            sv_range = f'{sv.chr}_{sv.start}-{sv.chr2}_{sv.end}'
                            supp_line = f'{genotype}:{length}:{dr}:{dv}:{sv_range}'
                            clique_formatstr_allsv[sv.sample].append(supp_line)
                    supp_vec = ''.join(supp_vec)
                    sample_supported = supp_vec.count('1')
                    # 接着按照样本顺序，生成所有样本中的这个clique的信息
                    # [sample1_clique1_info, sample2_clique1_info, ...]
                    clique_formatstrs_allsample = []
                    for sample in samples:
                        if sample not in clique_formatstr_allsv:
                            clique_formatstrs_allsample.append(
                                './.:0:0:0:.'
                            )
                        else:
                            clique_formatstrs_allsample.append(';'.join(clique_formatstr_allsv[sample]))
                    clique_final_format_str = '\t'.join(clique_formatstrs_allsample)
                    new_info_str = ';'.join([
                        f'SVTYPE={sv_type}',           # SV类型
                        f'AVGLEN={mean_len}',          # clique中所有sv的平均长度
                        f'CHR2={chr2}',                # BND类型的断点
                        f'END={end}',                  # SV的end位置
                        f'CIPOS={cipos}',              # 起始位置范围（以出现最多的SV为基准）
                        f'CIEND={ciend}',              # 位置范围（以出现最多的SV为基准）
                        f'SUPPORT={sample_supported}', # 支持clique的样本数量
                        f'SUPPORT_VEC={supp_vec}',     # 具体支持clique的样本向量
                    ])
                    newline = [
                        chr, start, 'sv_id', '.', f'<{sv_type}>', '.', 'PASS',
                        new_info_str,
                        'GT:LN:DR:DV:RG',
                        clique_final_format_str
                    ]
                    total_wait_dict[chrom].append(newline)
    # sort by start
    for chrom in total_wait_dict:
        # 按照start位置排序
        sorted_chrom = sorted(total_wait_dict[chrom], key=lambda x: x[1])
        total_wait_dict[chrom] = sorted_chrom
    # give ID
    n = 0
    with gzip.open(output, 'wb') as VCF_OUT:
        VCF_OUT.write(header.encode('utf-8'))
        for chrom in total_wait_dict:
            for sv in total_wait_dict[chrom]:
                sv_id = f'HN_{sv[4][1:-1]}_{n:05d}'
                n += 1
                sv[2] = sv_id
                sv[1] = str(sv[1])
                VCF_OUT.write(('\t'.join(sv) + '\n').encode('utf-8'))

def test():
    vcf_list = [
        "/home/wangdy/Projects/HZAU-Weibin/WGS_whole_blood/backup_20240624/WGS_wholeblood_20240427/SV/formated/wholeblood_7B.cutesv.vcf.gz",
        "/home/wangdy/Projects/HZAU-Weibin/WGS_whole_blood/backup_20240624/WGS_wholeblood_20240427/SV/formated/wholeblood_7B.pbsv.vcf.gz",
        "/home/wangdy/Projects/HZAU-Weibin/WGS_whole_blood/backup_20240624/WGS_wholeblood_20240427/SV/formated/wholeblood_7B.sniffles.vcf.gz"
    ]
    sample_names = ['cuteSV', 'pbsv', 'sniffles2']
    contig_lengths, _, SV_informations = read_vcf(vcf_list, sample_names)
    header_str = vcf_header(contig_lengths, sample_names)
    all_sv_block = group_by_block(SV_informations)
    all_sv_merged = block_merge(
        sv_blocks=all_sv_block,
        extended_dis=100,
        overlap_perc=0.5,
        worker=1
    )
    output_vcf(
        output="/home/wangdy/Projects/HZAU-Weibin/WGS_whole_blood/backup_20240624/WGS_wholeblood_20240427/SV/merged/wholeblood_7B.merged.vcf.gz",
        header=header_str,
        samples=sample_names,
        all_sv_merge=all_sv_merged
    )


def main():
    file_names = list(snakemake.input.keys())
    file_paths = list(snakemake.input)
    if len(file_names) != 0:
        contig_lengths, _, SV_informations = read_vcf(file_paths, file_names)
    else:
        contig_lengths, file_names, SV_informations = read_vcf(file_paths)
    header_str = vcf_header(contig_lengths, file_names)
    all_sv_block = group_by_block(SV_informations)
    all_sv_merged = block_merge(
        sv_blocks=all_sv_block,
        extended_dis=snakemake.params["extend_len"],
        overlap_perc=snakemake.params["overlap_rate"],
        worker=snakemake.threads
    )
    output_vcf(
        output=snakemake.output[0],
        header=header_str,
        samples=file_names,
        all_sv_merge=all_sv_merged
    )

if __name__ == "__main__":
    main()
