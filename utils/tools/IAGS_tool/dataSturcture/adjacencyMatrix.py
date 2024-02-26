import numpy as np
import pandas as pd


class AdjacencyMatrix:
    def __init__(self):
        self.adjacency_matrix = None

    def readFromMatrix(self, adjacency_matrix):
        self.adjacency_matrix = adjacency_matrix

    def readFromFile(self, file):
        self.adjacency_matrix = pd.read_csv(file, sep='\t', index_col=0, dtype=int)

    def output(self, outfile):
        self.adjacency_matrix.to_csv(outfile, sep='\t')

    def __reverse(self, item):
        if item[-1] == 'a':
            return item[:-1] + 'b'
        else:
            return item[:-1] + 'a'

    def assemble(self):
        index = self.adjacency_matrix.index.tolist()
        columns = self.adjacency_matrix.columns.tolist()
        np_adjacency_matrix = np.asarray(self.adjacency_matrix)
        adjacencies = {}
        for i in range(len(index)):
            for j in range(len(index)):
                if round(np_adjacency_matrix[i][j]) == 1:
                    if '$' == index[i] or '$' == index[j]:
                        continue
                    pair = sorted([index[i], index[j]])
                    key = pair[0] + '@' + pair[1]
                    if key not in adjacencies.keys():
                        adjacencies[key] = 1
                    else:
                        adjacencies[key] += 1
        adjs = {}
        for i in adjacencies.keys():
            itemset = i.split('@')
            if itemset[0] not in adjs.keys():
                adjs[itemset[0]] = itemset[1]
            if itemset[1] not in adjs.keys():
                adjs[itemset[1]] = itemset[0]

        startpoint = []

        for j in range(len(columns)):
            if np_adjacency_matrix[0][j] == 1:
                startpoint.append(columns[j])
        markerstartpoint = []
        self.chr = []
        for i in startpoint:
            if i not in markerstartpoint:
                path = []
                if i[-1] == 'a':
                    path.append(i[:-1])
                else:
                    path.append('-' + i[:-1])
                start = self.__reverse(i)
                if start in startpoint:
                    markerstartpoint.append(start)
                    self.chr.append(path)
                else:
                    while True:
                        next = adjs[start]
                        if next[-1] == 'a':
                            path.append(next[:-1])
                        else:
                            path.append('-' + next[:-1])
                        start = self.__reverse(next)
                        if start in startpoint:
                            markerstartpoint.append(start)
                            break
                    self.chr.append(path)
        vector = []
        for i in self.chr:
            for j in i:
                if j.startswith('-'):
                    vector.append(j[1:])
                else:
                    vector.append(j)
        cyclepoint = []
        for i in adjs.keys():
            if i[:-1] not in vector:
                cyclepoint.append(i)

        self.cyclechr = []
        markercycle = []
        for i in cyclepoint:
            if i not in markercycle:
                startpoint = i
                cycle = []
                markercycle.append(i)
                start = i
                while True:
                    next = adjs[start]
                    if next[-1] == 'a':
                        cycle.append(next[:-1])
                    else:
                        cycle.append('-' + next[:-1])
                    markercycle.append(start)
                    markercycle.append(next)
                    start = self.__reverse(next)
                    if start == startpoint:
                        break
                self.cyclechr.append(cycle)
        return self.chr, self.cyclechr

    def out_assembly(self, outfile, remove_bar):
        outfile = open(outfile, 'w')
        print("string:" + str(len(self.chr)))
        for i in self.chr:
            outfile.write('s ')
            for j in i:
                if remove_bar:
                    block = j.split('_')[0]
                    outfile.write(block + ' ')
                else:
                    outfile.write(j + ' ')
            outfile.write('\n')
        # print("cycle:" + str(len(self.cyclechr)))
        for k in self.cyclechr:
            outfile.write('c ')
            min_index = -1
            min_value = 1000000
            for l in range(len(k)):
                if k[l].startswith('-'):
                    item = k[l][1:].split('_')
                    block = int(item[0])
                else:
                    item = k[l].split('_')
                    block = int(item[0])
                if block < min_value:
                    min_index = l
                    min_value = block
            if k[min_index].startswith('-'):
                half1 = k[min_index + 1:]
                half2 = k[:min_index + 1]
                new_string = half1 + half2
            else:
                half1 = k[min_index:]
                half2 = k[:min_index]
                new_string = half1 + half2
            for l in new_string:
                if remove_bar:
                    block = l.split('_')[0]
                    outfile.write(block + ' ')
                else:
                    outfile.write(l + ' ')
            outfile.write('\n')
        outfile.close()
