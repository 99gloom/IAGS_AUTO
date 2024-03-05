from utils.pipline.processGraphFilter.graphAlgorithm.SelfSyntenyGraph import SelfSyntenyGraph
from utils.pipline.processGraphFilter.graphAlgorithm.MergedSyntenyGraph import MergedSyntenyGraph
from utils.pipline.processGraphFilter.findBack import FindBack
from utils.pipline.processGraphFilter.graphBlockFilter import GraphFilter
import os

def syntenyGraphMerge(species, ratio_path, order, drimmBlocks_dir, res_dir):
    root_name = 'root'
    virtual_name = 'virtual'
    # 获取物种拷贝数
    unsorted_sp_ratio = {}
    with open(ratio_path, 'r') as f:
        for line in f:
            temp = line.strip()
            temp = temp.split('\t')
            unsorted_sp_ratio[temp[0]] = int(temp[1])
    sp_copy_number = {}
    for i in species:
        sp_copy_number[i] = unsorted_sp_ratio[i]

    # 自匹配
    syn = {}
    for i in species:
        g = SelfSyntenyGraph(i, drimmBlocks_dir, sp_copy_number[i])
        syn[i] = g
    # 自底向上匹配
    for k, v in order.items():
        g = MergedSyntenyGraph(syn[v[0]], syn[v[1]])
        syn[k] = g
    # 存储结果
    syn['root'].out_res(res_dir, root_name, virtual_name)

def processGraphFilter(species, ratio_path, order, drimmBlocks_dir, merge_dir, sequences_dir,
                       synteny_path, manual_block_dir, unfiltered_block_dir, char_shape):
    # merge
    syntenyGraphMerge(species, ratio_path, order, drimmBlocks_dir, merge_dir)

    # merge结果映射回block
    finder = FindBack(merge_dir, drimmBlocks_dir, species)

    # 通过拷贝数过滤，得到finalBlock
    filterProcess = GraphFilter(ratio_path, finder.get_blocks_lists(), sequences_dir,
                                finder.get_expand_manual_synteny(os.path.join(synteny_path)), species, manual_block_dir,
                                unfiltered_block_dir, char_shape)
    filterProcess.matchLCS()
    filterProcess.processGenenumber()


    return filterProcess


