import gurobipy as gp
from gurobipy import *

import numpy as np
import pandas as pd

from ..dataSturcture.adjacencyMatrix import AdjacencyMatrix


class EndpointMatchingOptimization:

    def __init__(self, ancestor_file, guided_file, matching_dim1, matching_dim2, relation1, relation2, relabel=False):
        self.__matching_dim1 = matching_dim1
        self.__matching_dim2 = matching_dim2
        self.__relation1 = relation1
        self.__relation2 = relation2
        self.__relabel_block_sequences = []
        self.__relabel = relabel
        ancestor_matrix = pd.read_csv(ancestor_file, sep='\t', index_col=0)

        self.__ancestor_adjacency_list = self.__build_adjacency_list(ancestor_matrix)
        self.__candidate_compress_adjacency_matrix, self.__candidate_endpoint_list = \
            self.__build_assumed_matrix(self.__ancestor_adjacency_list)

        if self.__relabel:
            self.__guided_adjacency_list = self.__assumed_block_label(guided_file)
        else:
            self.__guided_adjacency_list = self.__build_adjacency_list_from_sequence(guided_file)

        self.__guided_compress_adjacency_matrix, self.__guided_endpoint_list = \
            self.__build_assumed_matrix(self.__guided_adjacency_list)

        self.__match_pairs = []
        for i in self.__candidate_compress_adjacency_matrix:
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
        self.__k = int((len(self.__candidate_endpoint_list) - 1) / (self.__matching_dim1))

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
                (i[0][2] * match_matrix[int(i[0][0] / (self.__matching_dim1)),
                                        i[0][0] % self.__matching_dim1,
                                        j[0] % self.__matching_dim2] + (1 - i[0][2])) *
                (i[0][3] * match_matrix[int(i[0][1] / (self.__matching_dim1)),
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
        for i in range(len(result)):
            column = []
            index = []
            for j in range(self.__matching_dim2):
                item = self.__candidate_endpoint_list[i * self.__matching_dim2 + 1].split('@')
                column.append(item[0] + '@' + str(j + 1))
            for j in range(self.__matching_dim1):
                item = self.__candidate_endpoint_list[i * self.__matching_dim2 + 1].split('@')
                index.append(item[0] + '@' + str(j + 1))
            match = pd.DataFrame(result[i], columns=column, index=index)
            match_relation = match.to_dict()
            for j in match_relation.keys():
                for k in match_relation[j].keys():
                    item = k.split('@')
                    if match_relation[j][k] == 1:
                        self.__match_relations[k] = j

    def output_matching_relation(self, outfile):
        outfile = open(outfile, 'w')
        for i in self.__match_relations.keys():
            line = i + ' ' + self.__match_relations[i]
            line += '\n'
            outfile.write(line)
        outfile.close()

    def build_adjacency_matrix(self):
        new_adjacency_list = []
        for i in self.__ancestor_adjacency_list:
            new_adjacency = []
            for j in i:
                if j == '$':
                    new_adjacency.append('$')
                else:
                    new_adjacency.append(self.__match_relations[j])
            new_adjacency_list.append(new_adjacency)

        new_result_Matrix = {}
        for i in self.__candidate_endpoint_list:
            if i == '$':
                key1 = '$'
            else:
                item1 = i.split('@')
                key1 = item1[0][:-1] + '_' + item1[1] + item1[0][-1]
            new_result_Matrix[key1] = {}
            for j in self.__candidate_endpoint_list:
                if j == '$':
                    key2 = '$'
                else:
                    item2 = j.split('@')
                    key2 = item2[0][:-1] + '_' + item2[1] + item2[0][-1]
                new_result_Matrix[key1][key2] = 0

        for i in new_adjacency_list:
            if i[0] == '$':
                key1 = '$'
            else:
                item1 = i[0].split('@')
                key1 = item1[0][:-1] + '_' + item1[1] + item1[0][-1]
            if i[1] == '$':
                key2 = '$'
            else:
                item2 = i[1].split('@')
                key2 = item2[0][:-1] + '_' + item2[1] + item2[0][-1]
            new_result_Matrix[key1][key2] = 1
            new_result_Matrix[key2][key1] = 1
        new_result_Matrix = pd.DataFrame(new_result_Matrix)
        self.adjacency_matrix = AdjacencyMatrix()
        self.adjacency_matrix.readFromMatrix(new_result_Matrix)
        return self.adjacency_matrix

    def __build_adjacency_list(self, ancestor_matrix):
        blocks_index = ancestor_matrix.columns.tolist()
        blocks_matrix = np.asarray(ancestor_matrix)
        endpoint_count = {}
        for i in blocks_index:
            endpoint_count[i] = 1
        process_adjacency = []
        adjacency_list = []
        for i in range(len(blocks_matrix)):
            for j in range(len(blocks_matrix[i])):
                pair1 = blocks_index[i]
                pair2 = blocks_index[j]
                key = [pair1, pair2]
                key = sorted(key)
                key = key[0] + '@' + key[1]
                if key not in process_adjacency:
                    process_adjacency.append(key)
                else:
                    continue
                if round(blocks_matrix[i][j]) != 0:
                    for k in range(int(round(blocks_matrix[i][j]))):
                        if pair1 == '$' and pair2 == '$':
                            adjacency_list.append([pair1, pair2])
                        elif pair1 == '$' and pair2 != '$':
                            adjacency_list.append([pair1, pair2 + '@' + str(endpoint_count[pair2])])
                            endpoint_count[pair2] += 1
                        elif pair1 != '$' and pair2 == '$':
                            adjacency_list.append([pair1 + '@' + str(endpoint_count[pair1]), pair2])
                            endpoint_count[pair1] += 1
                        else:
                            adjacency_list.append(
                                [pair1 + '@' + str(endpoint_count[pair1]), pair2 + '@' + str(endpoint_count[pair2])])
                            endpoint_count[pair1] += 1
                            endpoint_count[pair2] += 1

        return adjacency_list

    def __build_adjacency_list_from_sequence(self, file):
        adjacency_list = []
        with open(file, 'r') as df:
            while True:
                line = df.readline()[:-2]
                if not line:
                    break
                itemset = line.split(' ')
                chr_type = itemset[0]
                new_block_order = itemset[1:]
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
        return adjacency_list

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
        return adjacency_list

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
                if round(adjacency_matrix[i][j]) != 0:
                    for k in range(adjacency_matrix[i][j]):
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
