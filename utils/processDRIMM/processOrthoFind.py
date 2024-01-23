import os.path
from pathlib import Path

import numpy as np
import pandas as pd


class processOrthoFind:
    def __init__(self, OrthogroupsPath, gff_path_list, outputPath_tempDir, ratioPath):
        self.OrthogroupsPath = OrthogroupsPath
        self.gff_path_list = gff_path_list
        self.outputPath_tempDir = outputPath_tempDir
        self.ratioPath = ratioPath
        if (not Path(self.outputPath_tempDir).exists()):
            os.makedirs(self.outputPath_tempDir)
        self.__gffFileSort()
        self.default_dustLength = self.__get_dustLength(ratioPath)
        self.speciesAndChrLen, \
        self.sp, \
        self.rate_dir,\
        self.species_ratio= self.__execute()

    def __get_dustLength(self, ratioPath):
        sp_ratio_dir = {}
        with open(ratioPath, 'r') as f:
            for line in f:
                temp = line.rstrip('\n').rstrip(' ')
                temp = temp.split('\t')
                sp_ratio_dir[temp[0]] = temp[1]
        dustLength = 1
        for key, val in sp_ratio_dir.items():
            dustLength += int(val)
        return dustLength

    def __outSequence(self, sequence, outfile):
        outfile = open(outfile, 'w')
        for i in sequence:
            for j in i:
                outfile.write(str(j) + ' ')
            outfile.write('\n')
        outfile.close()

    def __getAllSequence(self, bed, group_dir):
        bed = pd.read_csv(bed, sep='\t', header=None)[[0, 1]]

        chrlist = bed[0].unique().tolist()
        # print(chrlist)
        new_sequence = []
        new_name_sequence = []
        for i in chrlist:
            split = bed.loc[bed[0] == i][1].tolist()

            chr = []
            chr_name = []
            for j in split:
                if j in group_dir.keys():
                    chr.append(group_dir[j])
                    chr_name.append(j)

            new_sequence.append(chr)
            new_name_sequence.append(chr_name)
        return new_sequence, new_name_sequence

    def __getFilterSequence(self, bed, group_filter_dir):
        bed = pd.read_csv(bed, sep='\t', header=None)[[0, 1]]

        chrlist = bed[0].unique().tolist()

        new_sequence = []
        for i in chrlist:
            split = bed.loc[bed[0] == i][1].tolist()

            chr = []
            for j in split:
                if j in group_filter_dir.keys():
                    chr.append(group_filter_dir[j])
            new_sequence.append(chr)
        return new_sequence

    def __gffFileSort(self):
        # 对gff文件内按chr顺序排序
        for i in self.gff_path_list:
            dataFrame = pd.read_csv(i, header=None, sep='\t')
            name = dataFrame.at[0, 0]
            number_index = 0
            for j in range(len(name) - 1, -1, -1):
                if name[j].isdigit():
                    continue
                else:
                    number_index = j + 1
                    break
            dataFrame[0] = dataFrame[0].map(lambda x: int(x[number_index:]))

            dataFrame = dataFrame.sort_values(by=[0, 2], ascending=[True, True])
            dataFrame[0] = dataFrame[0].map(lambda x: name[:number_index] + str(x))
            dataFrame[1] = dataFrame[1].map(lambda x: x.split(';')[0])
            dataFrame.to_csv(i, header=None, index=None, sep='\t')

    def __execute(self):
        sp = []
        gff_list = []
        unsorted_sp_ratio = {}
        sp_ratio = []

        for i in self.gff_path_list:
            gff_list.append(i.split('/')[-1])
            sp.append(i.split('/')[-1].split('.')[0])

        with open(self.ratioPath, 'r') as f:
            for line in f:
                temp = line.rstrip('\n').rstrip(' ')
                temp = temp.split('\t')
                unsorted_sp_ratio[temp[0]] = temp[1]
        for i in sp:
            sp_ratio.append(int(unsorted_sp_ratio[i]))
        sp_ratio_str = ':'
        sp_ratio_str = sp_ratio_str.join(list(map(str,sp_ratio)))


        ortho = pd.read_csv(self.OrthogroupsPath, sep='\t')
        ortho_columns = ortho.columns.tolist()[1:]
        ortho = ortho.fillna('')
        ortho = np.asarray(ortho)

        group_dir = {}
        for i in ortho:
            group = i[0]
            group_dir[group] = {}
            species = i[1:]
            for j in range(len(species)):
                genes = species[j].split(', ')
                if genes[0] == '':
                    group_dir[group][ortho_columns[j]] = []
                else:
                    group_dir[group][ortho_columns[j]] = genes

        rate_dir = {}
        all_rate_dir = {}
        finalGroup = {}

        for i in group_dir.keys():
            rate_list = []
            for j in sp:
                rate_list.append(len(group_dir[i][j]))
            ok = 1
            # 这里的筛选条件不一样
            for j in range(len(rate_list)):
                if rate_list[j] == 0:
                    ok = 0
            if ok == 0:
                continue
            else:
                rate = ''
                for j in rate_list:
                    rate += str(j) + ':'
                rate = rate[:-1]
                finalGroup[i] = group_dir[i]
                if rate not in all_rate_dir.keys():
                    all_rate_dir[rate] = 1
                else:
                    all_rate_dir[rate] += 1


        for i in group_dir.keys():
            rate_list = []
            for j in sp:
                rate_list.append(len(group_dir[i][j]))
            ok = 1
            # 这里的筛选条件不一样
            for j in range(len(rate_list)):
                if rate_list[j] > sp_ratio[j] or rate_list[j] == 0:
                    ok = 0
            if ok == 0:
                continue
            else:
                rate = ''
                for j in rate_list:
                    rate += str(j) + ':'
                rate = rate[:-1]
                finalGroup[i] = group_dir[i]
                if rate not in rate_dir.keys():
                    rate_dir[rate] = 1
                else:
                    rate_dir[rate] += 1


        outfile = self.outputPath_tempDir + '/group.xls'
        count = 1
        outfile = open(outfile, 'w')
        outfile.write('gene\tgroup\n')

        outfile_filter = self.outputPath_tempDir + '/filter_group.xls'
        outfile_filter = open(outfile_filter, 'w')
        outfile_filter.write('gene\tgroup\n')

        for i in group_dir.keys():
            for j in group_dir[i].keys():
                for k in group_dir[i][j]:
                    outfile.write(k + '\t' + str(count) + '\n')

            if i in finalGroup.keys():
                for j in finalGroup[i].keys():
                    for k in finalGroup[i][j]:
                        outfile_filter.write(k + '\t' + str(count) + '\n')
            count += 1

        outfile.close()
        outfile_filter.close()

        group = pd.read_csv(self.outputPath_tempDir + '/group.xls', sep='\t')
        group = np.asarray(group)
        group_dir = {}
        for i in group:
            group_dir[i[0]] = i[1]

        group_filter = pd.read_csv(self.outputPath_tempDir + '/filter_group.xls', sep='\t')
        group_filter = np.asarray(group_filter)
        group_filter_dir = {}
        for i in group_filter:
            group_filter_dir[i[0]] = i[1]

        sample_sequence_files = self.outputPath_tempDir + '/sample.sequence'
        sample_sequence_files = open(sample_sequence_files, 'w')
        speciesAndChrLen = {}

        for i in range(len(sp)):
            sequence, sequence_name = self.__getAllSequence(self.gff_path_list[i], group_dir)
            filter_sequence = self.__getFilterSequence(self.gff_path_list[i], group_filter_dir)
            speciesAndChrLen[sp[i]] = len(filter_sequence)
            item = gff_list[i].split('.')
            outfile = self.outputPath_tempDir + '/' + item[0] + '.sequence'
            self.__outSequence(filter_sequence, outfile)
            outallfile = self.outputPath_tempDir + '/' + item[0] + '.all.sequence'
            self.__outSequence(sequence, outallfile)
            outallfilename = self.outputPath_tempDir + '/' + item[0] + '.all.sequence.genename'
            self.__outSequence(sequence_name, outallfilename)

            for j in filter_sequence:
                for k in j:
                    sample_sequence_files.write(str(k) + ' ')
                sample_sequence_files.write('\n')

        sample_sequence_files.close()
        return speciesAndChrLen, sp, all_rate_dir, sp_ratio_str
