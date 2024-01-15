import gurobipy as gp
from gurobipy import *

import numpy as np
import pandas as pd


class BlockMatchingOptimization:

    def __init__(self, ancestor_file, guided_file,
                 matching_dim1, matching_dim2,
                 relation1, relation2, self_matching=False):
        self.__matching_dim1 = matching_dim1
        self.__matching_dim2 = matching_dim2
        self.__relation1 = relation1
        self.__relation2 = relation2
        self.__relabel_block_sequences = []
        self.__self_matching = self_matching

        self.__ancestor_adjacency_list, self.__ancestor_block_order = self.__assumed_block_label(ancestor_file)
        self.__ancestor_compress_adjacency_matrix, self.__ancestor_endpoint_list = \
            self.__build_assumed_matrix(self.__ancestor_adjacency_list)

        self.__guided_adjacency_list, self.__guided_block_order = self.__assumed_block_label(guided_file)
        self.__guided_compress_adjacency_matrix, self.__guided_endpoint_list = \
            self.__build_assumed_matrix(self.__guided_adjacency_list)

        self.__match_pairs = []
        for i in self.__ancestor_compress_adjacency_matrix:
            match_pair = []
            adj1 = i[-1]
            compare_key1 = ''
            if adj1[0] != '$':
                item = adj1[0].split('@')
                compare_key1 += item[0]
            if adj1[1] != '$':
                item = adj1[1].split('@')
                compare_key1 += item[0]

            for j in self.__guided_compress_adjacency_matrix:
                adj2 = j[-1]
                compare_key2 = ''
                if adj2[0] != '$':
                    item = adj2[0].split('@')
                    compare_key2 += item[0]
                if adj2[1] != '$':
                    item = adj2[1].split('@')
                    compare_key2 += item[0]
                if compare_key1 == compare_key2:
                    match_pair.append(j)
            self.__match_pairs.append([i, match_pair])
        self.__k = int((len(self.__ancestor_endpoint_list) - 1) / (self.__matching_dim1 * 2))

    def optimization(self):
        try:
            self.__m = gp.Model()
            match_matrix = self.__m.addVars(self.__k,
                                            self.__matching_dim1,
                                            self.__matching_dim2,
                                            vtype=GRB.BINARY,
                                            name="matching_matrix")
            self.__m.update()
            self.__m.setObjective(gp.quicksum(
                (i[0][2] * match_matrix[int(i[0][0] / (self.__matching_dim1 * 2)),
                                        i[0][0] % self.__matching_dim1,
                                        j[0] % self.__matching_dim2] + (1 - i[0][2])) *
                (i[0][3] * match_matrix[int(i[0][1] / (self.__matching_dim1 * 2)),
                                        i[0][1] % self.__matching_dim1,
                                        j[1] % self.__matching_dim2] + 1 - i[0][3])
                for i in self.__match_pairs for j in i[1]
            ), GRB.MAXIMIZE)
            self.__m.addConstrs((
                gp.quicksum(match_matrix[i, j, l] for l in range(self.__matching_dim2)) == self.__relation1
                for i in range(self.__k)
                for j in range(self.__matching_dim1)), name='row_unique'
            )
            self.__m.addConstrs((
                gp.quicksum(match_matrix[i, l, j] for l in range(self.__matching_dim1)) == self.__relation2
                for i in range(self.__k)
                for j in range(self.__matching_dim2)), name='col_unique'
            )
            if self.__self_matching:
                self.__m.addConstrs((
                    match_matrix[i, j, j] == 0
                    for i in range(self.__k)
                    for j in range(self.__matching_dim1)), name='diagonal'
                )

                self.__m.addConstrs((
                    match_matrix[i, j, k] == match_matrix[i, k, j]
                    for i in range(self.__k)
                    for j in range(self.__matching_dim1)
                    for k in range(self.__matching_dim2)), name='symmetry'
                )

            self.__m.optimize()
            print('Obj: %g' % self.__m.objVal)
        except gp.GurobiError as e:
            print('Error code ' + str(e.errno) + ': ' + str(e))
        except AttributeError:
            print('Encountered an attribute error')

    def matching_relation(self):
        result = []
        for v in self.__m.getVars():
            result.append(v.x)
        result = np.asarray(result)
        result = result.reshape((self.__k, self.__matching_dim1, self.__matching_dim2))
        self.__match_relations = {}

        if self.__self_matching:
            for i in range(self.__k):
                block = self.__ancestor_endpoint_list[i * self.__matching_dim1 * 2 + 1].split('@')[0][:-1]
                block_index = 1
                for j in range(self.__matching_dim1):
                    for k in range(len(result[i][j][j:])):
                        if result[i][j][k + j] == 1:
                            ancestor_block = block + '_' + str(j + 1)
                            guided_block = block + '_' + str(k + j + 1)
                            self.__match_relations[guided_block] = block + '_' + str(block_index)
                            self.__match_relations[ancestor_block] = block + '_' + str(block_index)
                            block_index += 1
        else:
            for i in range(self.__k):
                block = self.__ancestor_endpoint_list[i * self.__matching_dim1 * 2 + 1].split('@')[0][:-1]
                for j in range(self.__matching_dim1):
                    for k in range(self.__matching_dim2):
                        if result[i][j][k] == 1:
                            ancestor_block = block + '_' + str(j + 1)
                            guided_block = block + '_' + str(k + 1)
                            self.__match_relations[ancestor_block] = guided_block

    def output_matching_relation(self, outfile):
        outfile = open(outfile, 'w')
        for i in self.__match_relations.keys():
            key1 = i.split('_')

            line = key1[0] + ' ' + key1[1]
            key2 = self.__match_relations[i].split('_')
            line += ' ' + key2[1]
            line += '\n'
            outfile.write(line)
        outfile.close()

    def __assumed_block_label(self, file):
        adjacency_list = []
        block_objects = {}
        relabel_block_order = []
        with open(file) as df:
            while True:
                line = df.readline()[:-2]
                if not line:
                    break
                itemset = line.split(' ')
                chr_type = itemset[0]
                new_block_order = []
                for i in itemset[1:]:
                    block = ''
                    if i.startswith('-'):
                        block = i[1:]
                        if block not in block_objects.keys():
                            block_objects[block] = 1
                            new_block = '-' + block + '_' + str(block_objects[block])
                            block_objects[block] += 1
                        else:
                            new_block = '-' + block + '_' + str(block_objects[block])
                            block_objects[block] += 1
                    else:
                        block = i
                        if block not in block_objects.keys():
                            block_objects[block] = 1
                            new_block = block + '_' + str(block_objects[block])
                            block_objects[block] += 1
                        else:
                            new_block = block + '_' + str(block_objects[block])
                            block_objects[block] += 1
                    new_block_order.append(new_block)
                last = ''
                start = ''
                for i in range(len(new_block_order)):
                    if i == 0:
                        if chr_type == 's':
                            if new_block_order[i].startswith('-'):
                                block = new_block_order[i][1:].split('_')[0]
                                copy_number = new_block_order[i][1:].split('_')[1]
                                adjacency_list.append(['$', block + 'b' + '@' + copy_number])
                                last = block + 'a' + '@' + copy_number
                            else:
                                block = new_block_order[i].split('_')[0]
                                copy_number = new_block_order[i].split('_')[1]
                                adjacency_list.append(['$', block + 'a' + '@' + copy_number])
                                last = block + 'b' + '@' + copy_number
                        else:
                            if new_block_order[i].startswith('-'):
                                block = new_block_order[i][1:].split('_')[0]
                                copy_number = new_block_order[i][1:].split('_')[1]
                                last = block + 'a' + '@' + copy_number
                                start = block + 'b' + '@' + copy_number
                            else:
                                block = new_block_order[i].split('_')[0]
                                copy_number = new_block_order[i].split('_')[1]
                                last = block + 'b' + '@' + copy_number
                                start = block + 'a' + '@' + copy_number
                    else:
                        if new_block_order[i].startswith('-'):
                            block = new_block_order[i][1:].split('_')[0]
                            copy_number = new_block_order[i][1:].split('_')[1]
                            adjacency_list.append([last, block + 'b' + '@' + copy_number])
                            last = block + 'a' + '@' + copy_number
                        else:
                            block = new_block_order[i].split('_')[0]
                            copy_number = new_block_order[i].split('_')[1]
                            adjacency_list.append([last, block + 'a' + '@' + copy_number])
                            last = block + 'b' + '@' + copy_number
                if chr_type == 's':
                    adjacency_list.append([last, '$'])
                else:
                    adjacency_list.append([last, start])
                relabel_block_order.append([chr_type] + new_block_order)
        return adjacency_list, relabel_block_order

    def __build_assumed_matrix(self, adjacency_list):
        endpoint_list = []
        for i in adjacency_list:
            for j in i:
                if j not in endpoint_list:
                    endpoint_list.append(j)
        endpoint_list = sorted(endpoint_list)
        adjacency_matrix = {}
        for i in endpoint_list:
            adjacency_matrix[i] = {}
            for j in endpoint_list:
                adjacency_matrix[i][j] = 0
        for i in adjacency_list:
            adjacency_matrix[i[0]][i[1]] += 1
            adjacency_matrix[i[1]][i[0]] += 1
        adjacency_matrix = pd.DataFrame(adjacency_matrix)
        adjacency_matrix = np.asarray(adjacency_matrix)
        compress_adjacency_matrix = []
        for i in range(len(adjacency_matrix)):
            for j in range(len(adjacency_matrix[i])):
                if adjacency_matrix[i][j] == 1:
                    adjacency = [endpoint_list[i], endpoint_list[j]]

                    adjacency = sorted(adjacency)
                    if adjacency[0] == endpoint_list[i] and adjacency[1] == endpoint_list[j]:

                        if i == 0 and j == 0:
                            compress_adjacency_matrix.append([i, j, 0, 0, adjacency])
                        if i == 0 and j != 0:
                            compress_adjacency_matrix.append([i, j - 1, 0, 1, adjacency])
                        if i != 0 and j == 0:
                            compress_adjacency_matrix.append([i - 1, j, 1, 0, adjacency])
                        if i != 0 and j != 0:
                            compress_adjacency_matrix.append([i - 1, j - 1, 1, 1, adjacency])
                    if adjacency[0] == endpoint_list[j] and adjacency[1] == endpoint_list[i]:
                        if i == 0 and j == 0:
                            compress_adjacency_matrix.append([j, i, 0, 0, adjacency])
                        if i == 0 and j != 0:
                            compress_adjacency_matrix.append([j - 1, i, 1, 0, adjacency])
                        if i != 0 and j == 0:
                            compress_adjacency_matrix.append([j, i - 1, 0, 1, adjacency])
                        if i != 0 and j != 0:
                            compress_adjacency_matrix.append([j - 1, i - 1, 1, 1, adjacency])
        return compress_adjacency_matrix, endpoint_list

    def out_relabel_sequence(self, output_sequence_file):
        if self.__self_matching:
            outfile = output_sequence_file[0]
            outfile = open(outfile, 'w')
            for i in self.__ancestor_block_order:
                chr_type = i[0]
                outfile.write(chr_type + ' ')
                block_sequence = i[1:]
                relabel_block_sequence = []
                for j in block_sequence:
                    if j.startswith('-'):
                        block = j[1:]
                        relabel_block_sequence.append('-' + self.__match_relations[block])
                    else:
                        block = j
                        relabel_block_sequence.append(self.__match_relations[block])
                for j in relabel_block_sequence:
                    outfile.write(j + ' ')
                outfile.write('\n')
            outfile.close()
        else:
            relabel_block_sequences = []
            for i in self.__ancestor_block_order:
                relabel_block_sequence = []
                chr_type = i[0]
                for j in i[1:]:
                    if j.startswith('-'):
                        block = j[1:]
                        relabel_block_sequence.append('-' + self.__match_relations[block])
                    else:
                        block = j
                        relabel_block_sequence.append(self.__match_relations[block])
                relabel_block_sequences.append([chr_type] + relabel_block_sequence)
            candidate_file = open(output_sequence_file[0], 'w')
            guided_file = open(output_sequence_file[1], 'w')
            for i in relabel_block_sequences:
                line = ''
                for j in i:
                    line += j + ' '
                line += '\n'
                candidate_file.write(line)
            candidate_file.close()
            for i in self.__guided_block_order:
                line = ''
                for j in i:
                    line += j + ' '
                line += '\n'
                guided_file.write(line)
            guided_file.close()
