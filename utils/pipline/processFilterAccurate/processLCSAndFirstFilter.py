import os
# from LCS import LCS   # 调试用
from pathlib import Path

from utils.pipline.processFilterAccurate.LCS import LCS


class processLCSAndFirstFilter:
    def __init__(self, outdir, outputDir_rawBlock, ratioPath, sequenceDir, SyntenyPath, sp, chr_shape='s'):
        '''
        对DRIMM生成的block文件与原序列匹配，初次过滤不需要的block

        :param outdir: DRIMM拆分后的block文件夹
        :param outputDir_rawBlock: 初次过滤输出的block文件夹
        :param ratioPath: ratio比例文件路径
        :param sequenceDir: processOrthoFind得到的sequence所在的文件夹
        :param SyntenyPath: DRIMM产生的所有Synteny文件路径
        :param sp: 物种名列表
        :param chr_shape: 染色体类型 c/s
        '''
        self.outdir_synteny = outdir
        self.outdir_rawBlock = outputDir_rawBlock
        self.sp = sp
        self.ratio = ''
        self.sequenceDir = sequenceDir
        self.SyntenyPath = SyntenyPath
        self.chr_shape = chr_shape

        # 获取species列表
        # for i in gff_path_list:
        #     self.sp.append(i.split('/')[-1].split('.')[0])

        # 获取sp之间的拷贝数比例
        ratioDir = {}
        with open(ratioPath, 'r') as rf:
            for line in rf:
                content = line.rstrip('\n')
                content = content.split('\t')
                ratioDir[content[0]] = content[1]
        for i in self.sp:
            self.ratio += ratioDir[i] + ':'
        self.ratio = self.ratio[:-1]

        if not Path(self.outdir_rawBlock).exists():
            os.makedirs(self.outdir_rawBlock)

        self.__excute()

    def __readBlockSequence(self, filePath):
        # 读入各物种的原始block信息
        chr = []
        chr_list = []
        chr_count = 1
        fr = open(filePath, 'r')
        for line in fr:
            line = line[:-2]
            itemset = line.split(' ')[1:]
            chr.append(itemset)
            chr_list.append(str(chr_count))
            chr_count += 1
        chr_Dict = {}
        for i in range(len(chr_list)):
            chr_Dict[chr_list[i]] = chr[i]
        fr.close()
        return chr_Dict

    def __readOriginSequence(self, filePath):
        # 读入原始的基因ID序列
        chr = []
        chr_list = []
        chr_count = 1
        fr = open(filePath, 'r')
        for line in fr:
            line = line[:-2]
            itemset = line.split(' ')
            chr.append(itemset)
            chr_list.append(str(chr_count))
            chr_count += 1
        chr_Dict = {}
        for i in range(len(chr_list)):
            chr_Dict[chr_list[i]] = chr[i]
        fr.close()
        return chr_Dict

    def __syntenyDict(self, filePath):
        syntenyDict = {}
        with open(filePath, 'r') as sf:
            for line in sf:
                temp = line.rstrip('\n').rstrip()
                itemset = temp.split(' ')
                header = itemset[0].split(':')
                syntenyDict[header[0]] = itemset[1:]

        return syntenyDict

    def __syntenyCount(self, filePath):
        syntenyDict = {}
        with open(filePath) as sf:
            for line in sf:
                temp = line.rstrip('\n').rstrip()
                itemset = temp.split(' ')
                header = itemset[0].split(':')
                if header[0] not in syntenyDict.keys():
                    syntenyDict[header[0]] = len(itemset[1:])
                else:
                    syntenyDict[header[0]] += len(itemset[1:])

        return syntenyDict

    def __assambleDrimmSequence(self, blockSequence, synteny):
        # 将Block序列通过synteny文件还原为基因ID形式
        sequences = {}
        sequences_ID = {}
        blockCount = {}
        for i in blockSequence.keys():
            sequence = []
            sequence_ID = []
            for j in blockSequence[i]:
                if j.startswith('-'):
                    block = j[1:]
                    synteny_sequence = synteny[block][::-1]
                    if block not in blockCount.keys():
                        blockCount[block] = 1
                    else:
                        blockCount[block] += 1
                    for k in range(len(synteny_sequence)):
                        sequence.append(synteny_sequence[k])
                        sequence_ID.append('-' + block + '|' + str(blockCount[block]) + '|' + str(k))
                else:
                    block = j
                    synteny_sequence = synteny[block]
                    if block not in blockCount.keys():
                        blockCount[block] = 1
                    else:
                        blockCount[block] += 1
                    for k in range(len(synteny_sequence)):
                        sequence.append(synteny_sequence[k])
                        sequence_ID.append('+' + block + '|' + str(blockCount[block]) + '|' + str(k))
            sequences[i] = sequence
            sequences_ID[i] = sequence_ID
        return sequences, sequences_ID

    def __matchingSequence(self, species_all_sequences, species_reassamble_sequences, species_all_sequences_name,
                           species_reassamble_sequences_ID):
        # 将两个序列做匹配
        block_range = {}
        for i in species_all_sequences.keys():
            species_reassamble_sequence = species_reassamble_sequences[i]
            # print(i)
            block_range[i] = {}
            for j in species_all_sequences[i].keys():
                if j not in species_reassamble_sequence.keys():
                    continue
                else:
                    # print(j)
                    p = LCS()
                    p.input(species_all_sequences[i][j]
                            , species_reassamble_sequence[j])
                    direction_list, lcslength_list = p.Compute_LCS()
                    lcs = p.printOneLCS()
                    for k in lcs:
                        genename = species_all_sequences_name[i][j][k[0]]
                        ID = species_reassamble_sequences_ID[i][j][k[1]]
                        ID_split = ID.split('|')
                        block = ID_split[0][1:]
                        block_count = ID_split[1]
                        block_stand = ID_split[0][0]
                        if block not in block_range[i].keys():
                            block_range[i][block] = {}
                            block_range[i][block][block_count + '@' + block_stand] = [[j, k[0], genename]]
                        else:
                            if block_count + '@' + block_stand not in block_range[i][block].keys():
                                block_range[i][block][block_count + '@' + block_stand] = [[j, k[0], genename]]
                            else:
                                block_range[i][block][block_count + '@' + block_stand].append([j, k[0], genename])
        return block_range

    def __outSynteny(self, block_range, species_all_sequences_name, species_all_sequences):
        # 输出synteny文件
        for i in block_range.keys():
            outfile = self.outdir_synteny + '/' + i + '.synteny'
            outfile_name = self.outdir_synteny + '/' + i + '.synteny.genename'
            outfile = open(outfile, 'w')
            outfile_name = open(outfile_name, 'w')
            for j in block_range[i].keys():
                block = j
                for k in block_range[i][j].keys():
                    block_count = k.split('@')[0]
                    block_stand = k.split('@')[1]
                    matching_pairs = block_range[i][j][k]
                    # print(matching_pairs)
                    matching_pairs = sorted(matching_pairs, key=lambda x: x[1])
                    chr = matching_pairs[0][0]
                    start = matching_pairs[0][1]
                    end = matching_pairs[-1][1]

                    if block_stand == '-':
                        genename = species_all_sequences_name[i][chr][start:end + 1][::-1]
                        genesequence = species_all_sequences[i][chr][start:end + 1][::-1]
                    else:
                        genename = species_all_sequences_name[i][chr][start:end + 1]
                        genesequence = species_all_sequences[i][chr][start:end + 1]
                    outfile.write(block + ':' + str(block_count) + ':chr_' + chr + ':' + block_stand + ' ')
                    outfile_name.write(block + ':' + str(block_count) + ':chr_' + chr + ':' + block_stand + ' ')
                    for l in genename:
                        outfile_name.write(l + ' ')
                    outfile_name.write('\n')
                    for l in genesequence:
                        outfile.write(l + ' ')
                    outfile.write('\n')
            outfile.close()
            outfile_name.close()

    def __excute(self):

        block_dir = self.outdir_synteny
        species = self.sp
        speciesRatio = self.ratio
        syntenyFile = self.SyntenyPath
        sequences_dir = self.sequenceDir

        blockSequences = {}
        for i in species:
            blockSequences[i] = self.__readBlockSequence(block_dir + '/' + i + '.block')
        synteny = self.__syntenyDict(syntenyFile)
        species_reassamble_sequences = {}
        species_reassamble_sequences_ID = {}
        for i in blockSequences.keys():
            sequences, sequences_ID = self.__assambleDrimmSequence(blockSequences[i], synteny)
            species_reassamble_sequences[i] = sequences
            species_reassamble_sequences_ID[i] = sequences_ID

        # 读取allsequence和name
        species_all_sequences = {}
        species_all_sequences_name = {}
        for i in species:
            # print(i)
            species_all_sequences[i] = self.__readOriginSequence(sequences_dir + '/' + i + '.all.sequence')
            species_all_sequences_name[i] = self.__readOriginSequence(
                sequences_dir + '/' + i + '.all.sequence.genename')

        # 两个序列做匹配
        block_range = self.__matchingSequence(species_all_sequences, species_reassamble_sequences,
                                              species_all_sequences_name, species_reassamble_sequences_ID)
        self.__outSynteny(block_range, species_all_sequences_name, species_all_sequences)

        # 通过synteny文件找满足条件的block，过滤不满足比例的
        blocks_ratio = {}
        block_list = []
        species_block_dict = {}
        for i in species:
            species_block_dict[i] = {}
            block_synteny_file = self.outdir_synteny + '/' + i + '.synteny'
            with open(block_synteny_file, 'r') as bsf:
                for line in bsf:
                    temp = line.rstrip('\n').rstrip()
                    itemset = temp.split(' ')
                    header = itemset[0].split(':')
                    block = header[0]
                    if block not in species_block_dict[i].keys():
                        species_block_dict[i][block] = 1
                    else:
                        species_block_dict[i][block] += 1
                    if block not in block_list:
                        block_list.append(block)

        for i in block_list:
            ratio = ''
            for j in species:
                if i not in species_block_dict[j].keys():
                    ratio += '0:'
                else:
                    ratio += str(species_block_dict[j][i]) + ':'
            ratio = ratio[:-1]
            blocks_ratio[i] = ratio

        filter_block = []
        for i in blocks_ratio.keys():
            if blocks_ratio[i] == speciesRatio:
                filter_block.append(i)

        for i in species:
            out_block_sequences = self.outdir_rawBlock + '/' + i + '.block'
            out_block_sequences = open(out_block_sequences, 'w')
            block_sequences = blockSequences[i]
            for j in block_sequences.keys():
                if self.chr_shape == 's' or self.chr_shape == 'S':
                    out_block_sequences.write('s ')
                else:
                    out_block_sequences.write('c ')

                for k in block_sequences[j]:
                    if k.startswith('-'):
                        block = k[1:]
                    else:
                        block = k
                    if block in filter_block:
                        out_block_sequences.write(k + ' ')
                out_block_sequences.write('\n')
            out_block_sequences.close()

        print("The following data is for reference only: ")
        for i in species:
            block_synteny_file = self.outdir_synteny + '/' + i + '.synteny'
            block_synteny_count = self.__syntenyCount(block_synteny_file)
            filter_synteny_count = {}
            syntenygenes_number = 0
            for j in block_synteny_count.keys():
                if j not in filter_block:
                    filter_synteny_count[j] = block_synteny_count[j]
                    syntenygenes_number += filter_synteny_count[j]
            allgenes_number = 0
            for j in species_all_sequences[i].keys():
                allgenes_number += len(species_all_sequences[i][j])
            print(i)
            print(syntenygenes_number / allgenes_number)
            print(allgenes_number)
