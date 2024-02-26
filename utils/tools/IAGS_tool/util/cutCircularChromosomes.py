import numpy as np
import pandas as pd

from ..inferringAncestorGenomeStructure.BlockMatchingOptimization import BlockMatchingOptimization


class StatisticsAdjacencies:

    def __init__(self, ancestor_file, sp_filelist, ancestor_name, splist):
        ancestor_adj = self.__sequence2Adj(ancestor_file)
        sp_adjs = []
        for i in sp_filelist:
            sp_adjs.append(self.__sequence2Adj(i))
        ancestor_adj_dir = {}
        for i in ancestor_adj:
            item = sorted(i)
            key = item[0] + '@' + item[1]
            if key not in ancestor_adj_dir.keys():
                ancestor_adj_dir[key] = 1
            else:
                ancestor_adj_dir[key] += 1
        sp_adjs_dir = []
        for i in sp_adjs:
            sp_adj_dir = {}
            sp_adj = i
            for j in sp_adj:
                item = sorted(j)
                key = item[0] + '@' + item[1]
                if key not in sp_adj_dir.keys():
                    sp_adj_dir[key] = 1
                else:
                    sp_adj_dir[key] += 1
            sp_adjs_dir.append(sp_adj_dir)
        self.__statistic = []
        columns = [ancestor_name] + splist
        index = list(ancestor_adj_dir.keys())
        for i in ancestor_adj_dir.keys():
            item = [ancestor_adj_dir[i]]
            for j in range(len(sp_adjs_dir)):
                if i in sp_adjs_dir[j].keys():
                    item.append(sp_adjs_dir[j][i])
                else:
                    item.append(0)
            self.__statistic.append(item)
        self.__statistic = pd.DataFrame(self.__statistic, index=index, columns=columns)

    def out_ev_file(self, outfile):
        self.__statistic.to_csv(outfile, sep='\t')

    def __sequence2Adj(self, sequence_file):
        sequence = []
        with open(sequence_file, 'r') as sf:
            while True:
                line = sf.readline()[:-2]
                if not line:
                    break
                itemset = line.split(' ')
                sequence.append(itemset)
        adj = []
        for i in sequence:
            chr_type = i[0]
            block_sequence = i[1:]
            last = ''
            start = ''
            for j in range(len(block_sequence)):
                if j == 0:
                    if chr_type == 's':
                        block = block_sequence[j]
                        if block.startswith('-'):
                            adj.append(['$', block[1:] + 'b'])
                            last = block[1:] + 'a'
                        else:
                            adj.append(['$', block + 'a'])
                            last = block + 'b'
                    else:
                        block = block_sequence[j]
                        if block.startswith('-'):
                            last = block[1:] + 'a'
                            start = block[1:] + 'b'
                        else:
                            last = block + 'b'
                            start = block + 'a'

                else:
                    block = block_sequence[j]
                    if block.startswith('-'):
                        adj.append([last, block[1:] + 'b'])
                        last = block[1:] + 'a'
                    else:
                        adj.append([last, block + 'a'])
                        last = block + 'b'
            if chr_type == 's':
                adj.append([last, '$'])
            else:
                adj.append([last, start])
        return adj


def readSequence(file):
    chr = []
    with open(file, 'r') as rf:
        while True:
            line = rf.readline()[:-2]
            if not line:
                break
            itemset = line.split(' ')
            chr.append(itemset)
    return chr


def reverse(item):
    if item[-1] == 'a':
        return item[:-1] + 'b'
    else:
        return item[:-1] + 'a'


def reassembly(adjacencies, startpoint, endpoint):
    adjs = {}
    for i in adjacencies:
        if i[0] not in adjs.keys():
            adjs[i[0]] = i[1]
        if i[1] not in adjs.keys():
            adjs[i[1]] = i[0]
    path = []
    if startpoint[-1] == 'a':
        path.append(startpoint[:-1])
    else:
        path.append('-' + startpoint[:-1])
    start = reverse(startpoint)
    while True:
        next = adjs[start]
        if next[-1] == 'a':
            path.append(next[:-1])
        else:
            path.append('-' + next[:-1])
        start = reverse(next)
        if start == endpoint:
            break
    return path


def outSequence(sequence, outfile):
    outfile = open(outfile, 'w')
    for i in sequence:
        outfile.write('s ')
        for j in i[1:]:
            item = j.split('_')
            outfile.write(item[0] + ' ')
        outfile.write('\n')


def transformToAdjacency_for_scoj(file):
    adjacency_list = []
    with open(file) as df:
        while True:
            line = df.readline()[:-2]
            if not line:
                break
            item = line.split(' ')
            chr_type = item[0]
            last = ''
            start = ''
            sequence = item[1:]
            for j in range(len(sequence)):
                if j == 0:
                    if chr_type == 's':
                        if sequence[j].startswith('-'):
                            adjacency_list.append(['$', sequence[j][1:] + 'b'])
                            last = sequence[j][1:] + 'a'
                        else:
                            adjacency_list.append(['$', sequence[j] + 'a'])
                            last = sequence[j] + 'b'
                    else:
                        if sequence[j].startswith('-'):
                            last = sequence[j][1:] + 'a'
                            start = sequence[j][1:] + 'b'
                        else:
                            last = sequence[j] + 'b'
                            start = sequence[j] + 'a'
                else:
                    if sequence[j].startswith('-'):
                        adjacency_list.append([last, sequence[j][1:] + 'b'])
                        last = sequence[j][1:] + 'a'
                    else:
                        adjacency_list.append([last, sequence[j] + 'a'])
                        last = sequence[j] + 'b'
            if chr_type == 's':
                adjacency_list.append([last, '$'])
            else:
                adjacency_list.append([last, start])
    new_adjacency_list = []
    for j in adjacency_list:
        new_adjacency_list.append(j[0] + '@' + j[1])
    return new_adjacency_list


def cutCircularChromosomes(ancestor_file, ancestor_copy_number, ancestor_name,
                           speciesAndCopyList, outdir):
    """
    Firstly, calculating ancestral adjacencies support table.
    All species used for calculating this ancestor should first match with this ancestor by BMO integer programming formulations
    and then counting the number of block adjacencies appeared in all species, respectively.
    Secondly, ancestors may appear circular chromosomes.
    Here, IAGS allows user to cut one adjacency in circular chromosomes with minimum number of support to make circular to strings.
    :param ancestor_file: block sequence file for inferred ancestor
    :param ancestor_copy_number: target copy number of ancestor
    :param ancestor_name: ancestor name
    :param speciesAndCopyList: all species block sequences file,target copy number and species name
    :param outdir: output directory
    """
    matchingFileList = []
    speciesName = []
    ancestormatchingFile = ''
    sumcopynumber = 0
    # matching with ancestor
    for i in speciesAndCopyList:
        sumcopynumber += i[1]
        mo = BlockMatchingOptimization(i[0],
                                       ancestor_file,
                                       matching_dim1=i[1],
                                       matching_dim2=ancestor_copy_number,
                                       relation1=ancestor_copy_number / ancestor_copy_number,
                                       relation2=i[1] / ancestor_copy_number)
        mo.optimization()
        mo.matching_relation()
        output_sequence_file_list = [outdir + i[2] + '.ev.relabel.block',
                                     outdir + ancestor_name + '.ev.relabel.block']
        mo.out_relabel_sequence(output_sequence_file_list)
        matchingFileList.append(outdir + i[2] + '.ev.relabel.block')
        speciesName.append(i[2])
        ancestormatchingFile = outdir + ancestor_name + '.ev.relabel.block'

    # calculating ancestral adjacencies support table
    ev = StatisticsAdjacencies(ancestormatchingFile,
                               matchingFileList,
                               ancestor_name, speciesName)
    ev.out_ev_file(outdir + 'ev_' + ancestor_name + '.xls')
    ev_table = pd.read_csv(outdir + 'ev_' + ancestor_name + '.xls', sep='\t', index_col=0)
    adjs = ev_table.index.tolist()
    ev_table = np.asarray(ev_table)
    adjs_weight = {}
    for i in range(len(adjs)):
        adjs_weight[adjs[i]] = np.sum(ev_table[i][1:])

    # cut circular chromosomes
    ancestorSequence = readSequence(ancestormatchingFile)
    stringSequence = []
    cycleSequence = []
    for i in ancestorSequence:
        if i[0] == 's':
            stringSequence.append(i)
        else:
            cycleSequence.append(i)
    for i in cycleSequence:
        adjacency_list = []
        last = ''
        start = ''
        sequence = i[1:]
        for j in range(len(sequence)):
            if j == 0:
                if sequence[j].startswith('-'):
                    last = sequence[j][1:] + 'a'
                    start = sequence[j][1:] + 'b'
                else:
                    last = sequence[j] + 'b'
                    start = sequence[j] + 'a'
            else:
                if sequence[j].startswith('-'):
                    adjacency_list.append([last, sequence[j][1:] + 'b'])
                    last = sequence[j][1:] + 'a'
                else:
                    adjacency_list.append([last, sequence[j] + 'a'])
                    last = sequence[j] + 'b'
        adjacency_list.append([last, start])
        # print(adjacency_list)
        new_adjacency_list = []
        minkey = ''
        minweight = 100000
        startpoint = ''
        endpoint = ''
        for j in adjacency_list:
            key = sorted(j)
            key = key[0] + '@' + key[1]
            if adjs_weight[key] < minweight:
                minweight = adjs_weight[key]
                minkey = key
        for j in adjacency_list:
            key = sorted(j)
            key = key[0] + '@' + key[1]
            if key != minkey:
                new_adjacency_list.append(j)
            else:
                new_adjacency_list.append(['$', j[0]])
                startpoint = j[0]
                new_adjacency_list.append(['$', j[1]])
                endpoint = j[1]
        path = reassembly(new_adjacency_list, startpoint, endpoint)
        stringSequence.append(['s'] + path)
    outSequence(stringSequence, outdir + ancestor_name + '.cutcycle.block')
