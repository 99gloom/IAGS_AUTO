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


def estimatedAccuracyModel(CRB_ratio, mode_type):
    if mode_type == 'GMP':
        acc = 0.6325 * CRB_ratio * CRB_ratio + 0.3116 * CRB_ratio
    elif mode_type == 'GGHP':
        acc = 0.5292 * CRB_ratio * CRB_ratio + 0.3669 * CRB_ratio
    elif mode_type == 'MultiCopyGMP':
        acc = 0.4310 * CRB_ratio * CRB_ratio + 0.5420 * CRB_ratio
    elif mode_type == 'MultiCopyGGHP':
        acc = 0.6217 * CRB_ratio * CRB_ratio + 0.2740 * CRB_ratio
    else:
        acc = 0
        print('input error')
    return acc


def calculatedCRBrateAndEstimationAccuracy(matching_target_file, matching_target_copy_number, matching_target_name,
                                           speciesAndCopyList, outdir, model_type):
    """
    All species should first match with this a target species by BMO integer programming formulations
    Target species is a species with small copy number in input species.
    (for example,
    species 1 with copy number 4 and species 2 with 2, the target species is species 1.
    species 1 and species 2 with both 2 copy number, then target species can be 1 or 2.)
    Then, IAGS calculates completely rearranged breakpoints ratio and
    obtains estimation accuracy by accuracy estimation function.
    :param matching_target_file: block sequence file for target species
    :param matching_target_copy_number: target copy number of target species
    :param matching_target_name: target species name
    :param speciesAndCopyList: all species block sequences file,target copy number and species name
    :param outdir: output directory
    :param model_type: model used for obtaining ancestor, including GMP, GGHP, MultiCopyGMP and MultiCopyGGHP
    """
    matchingFileList = []
    sumcopynumber = 0
    # matching with target species
    for i in speciesAndCopyList:
        sumcopynumber += i[1]
        mo = BlockMatchingOptimization(i[0],
                                       matching_target_file,
                                       matching_dim1=i[1],
                                       matching_dim2=matching_target_copy_number,
                                       relation1=matching_target_copy_number / matching_target_copy_number,
                                       relation2=i[1] / matching_target_copy_number)
        mo.optimization()
        mo.matching_relation()
        output_sequence_file_list = [outdir + i[2] + '.ev.relabel.block',
                                     outdir + matching_target_name + '.ev.relabel.block']
        mo.out_relabel_sequence(output_sequence_file_list)
        matchingFileList.append(outdir + i[2] + '.ev.relabel.block')

    endpoints = {}
    for i in matchingFileList:
        adj = transformToAdjacency_for_scoj(i)
        for j in adj:
            endpoint1 = j.split('@')[0]
            endpoint2 = j.split('@')[1]
            if endpoint1 not in endpoints.keys():
                endpoints[endpoint1] = [endpoint2]
            else:
                if endpoint2 not in endpoints[endpoint1]:
                    endpoints[endpoint1].append(endpoint2)

            if endpoint2 not in endpoints.keys():
                endpoints[endpoint2] = [endpoint1]
            else:
                if endpoint1 not in endpoints[endpoint2]:
                    endpoints[endpoint2].append(endpoint1)
    rate = 0

    for i in endpoints.keys():
        if i == '$':
            continue
        if len(endpoints[i]) == sumcopynumber / matching_target_copy_number:
            rate += 1
    rate = round(rate / (len(endpoints.keys()) - 1), 4)
    acc = estimatedAccuracyModel(rate, model_type)
    ev_file = open(outdir + 'ev.txt', 'w')
    ev_file.write('CRB ratio: ' + str(rate * 100) + '%' + '\n')
    ev_file.write('Mode: ' + model_type + '\n')
    ev_file.write('Estimated acc: ' + str(100 - acc * 100) + '%' + '\n')
    ev_file.close()
