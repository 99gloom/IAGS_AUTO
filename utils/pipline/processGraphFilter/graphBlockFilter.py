from typing import List

from utils.tools.LCS import LCS
from utils.tools.makeDirFunc import make_dir
import os

class GraphFilter():
    def __init__(self, ratioPath: str, block_lists: dict, sequences_dir:str, synteny_dict: dict,
                 species:list, manual_block_dir:str, unfilter_block_dir:str, chr_shape: str):
        self.manual_num = 100000
        self.manual_mapping = {}
        self.block_lists = block_lists
        self.synteny_dict = synteny_dict
        self.species = species
        self.chr_shape = chr_shape
        self.sequences_dir = sequences_dir
        self.manual_block_dir = manual_block_dir
        self. unfilter_block_dir = unfilter_block_dir
        self.ratio = ''
        ratioDir = {}
        with open(ratioPath, 'r') as rf:
            for line in rf:
                content = line.rstrip('\n')
                content = content.split('\t')
                ratioDir[content[0]] = content[1]
        for i in self.species:
            self.ratio += ratioDir[i] + ':'
        self.ratio = self.ratio[:-1]

        # make_dir
        make_dir(self.manual_block_dir)

    def _readMemoryBlock(self, block_list: list):
        chr = []
        chr_list = []
        chr_count = 1
        for line in block_list:
            chr.append(line)
            chr_list.append(str(chr_count))
            chr_count += 1
        chr_Dict = {}
        for i in range(len(chr_list)):
            chr_Dict[chr_list[i]] = chr[i]
        return chr_Dict

    def _readMemorySynteny(self, synteny_dict: dict):
        return synteny_dict

    def _readFileSequence(self, file_path):
        chr = []
        chr_list = []
        chr_count = 1
        fr = open(file_path, 'r')
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

    def _syntenyCount(self, synteny_dict: dict):
        syntenyDictCount = {}
        for header, itemset in synteny_dict.items():
            if header not in syntenyDictCount.keys():
                syntenyDictCount[header] = len(itemset)
            else:
                syntenyDictCount[header] += len(itemset)
        return syntenyDictCount

    def _assambleDrimmSequence(self, blockSequence,synteny):
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


    def _matchingSequence(self, species_all_sequences, species_reassemble_sequences, species_all_sequences_name,
                           species_reassemble_sequences_ID):
        # 两个序列做匹配
        block_range = {}
        for i in species_all_sequences.keys():
            species_reassemble_sequence = species_reassemble_sequences[i]
            block_range[i] = {}
            for j in species_all_sequences[i].keys():
                if j not in species_reassemble_sequence.keys():
                    continue
                else:
                    p = LCS()
                    p.input(species_all_sequences[i][j]
                            , species_reassemble_sequence[j])
                    direction_list, lcslength_list = p.Compute_LCS()
                    lcs = p.printOneLCS()
                    for k in lcs:
                        genename = species_all_sequences_name[i][j][k[0]]
                        ID = species_reassemble_sequences_ID[i][j][k[1]]
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

    def _outSynteny(self, block_range, species_all_sequences_name, species_all_sequences):
        for i in block_range.keys():
            outfile = os.path.join(self.manual_block_dir, i + '.manual.synteny')
            outfile_name = os.path.join(self.manual_block_dir, i + '.manual.synteny.genename')
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

    def processGenenumber(self):
        block_len = {}
        for i in self.species:
            with open(os.path.join(self.unfilter_block_dir, i + '.final.synteny'), 'r') as sf:
                for line in sf:
                    temp = line
                    temp = temp.rstrip('\n').rstrip()
                    block = temp.split(' ')[0].split(':')[0]
                    length = len(temp.split(' ')[1:])

                    if block not in block_len.keys():
                        block_len[block] = length
                    else:
                        if block_len[block] < length:
                            block_len[block] = length

        with open(os.path.join(self.unfilter_block_dir, 'blockindex.genenumber'), 'w') as f:
            f.write('blockID\tblockLength\n')
            for i, j in block_len.items():
                f.write(i + '\t' + str(j) + '\n')
        return

    def _get_manual_mapping_start(self, l:List[str]):
        max_val = 0
        for i in l:
            if i.isdigit():
                if max_val < int(i):
                    max_val = int(i)
        count = 0
        while (max_val != 0):
            max_val = max_val // 10
            count += 1
        return 10**(count + 1)

    def matchLCS(self):
        species_block = {}
        for i in self.species:
            species_block[i] = self._readMemoryBlock(self.block_lists[i])
        synteny = self._readMemorySynteny(self.synteny_dict)
        species_reassemble_sequences = {}
        species_reassemble_sequences_ID = {}
        for i in species_block.keys():
            sequences,sequences_ID = self._assambleDrimmSequence(species_block[i], synteny)
            species_reassemble_sequences[i] = sequences
            species_reassemble_sequences_ID[i] = sequences_ID

        # 读取allsequence和name
        species_all_sequences = {}
        species_all_sequences_name = {}
        for i in self.species:
            species_all_sequences[i] = self._readFileSequence(self.sequences_dir + i + '.all.sequence')
            species_all_sequences_name[i] = self._readFileSequence(self.sequences_dir + i + '.all.sequence.genename')

        # 两个序列做匹配
        block_range = self._matchingSequence(species_all_sequences, species_reassemble_sequences,
                                             species_all_sequences_name, species_reassemble_sequences_ID)
        self._outSynteny(block_range, species_all_sequences_name, species_all_sequences)

        #
        blocks_ratio = {}
        block_copy_number = {}
        block_list = []
        for i in self.species:
            block_copy_number[i] = {}
            block_sequences = species_block[i]
            manual_file = open((os.path.join(self.manual_block_dir, i + '.manual.block')), 'w')
            for chr_num, item in block_sequences.items():
                if self.chr_shape == 's':
                    manual_file.write('s ')
                else:
                    manual_file.write('c ')
                for j in item:
                    manual_file.write(j+' ')
                    if j.startswith('-'):
                        block = j[1:]
                    else:
                        block = j

                    if block not in block_copy_number[i].keys():
                        block_copy_number[i][block] = 1
                    else:
                        block_copy_number[i][block] += 1
                    if block not in block_list:
                        block_list.append(block)
                manual_file.write('\n')
            manual_file.close()
        for i in block_list:
            ratio = []
            for sp in self.species:
                if i not in block_copy_number[sp].keys():
                    ratio.append(0)
                else:
                    ratio.append(block_copy_number[sp][i])
            ratio_str = ':'.join(list(map(str,ratio)))
            blocks_ratio[i] = ratio_str

        required_block = []
        for i in blocks_ratio.keys():
            if blocks_ratio[i] == self.ratio:
                required_block.append(i)

        self.manual_num = self._get_manual_mapping_start(required_block)
        manual_block = sorted([i for i in required_block if not i.isdigit()])
        for i in manual_block:
            self.manual_mapping[i] = str(self.manual_num)
            self.manual_num += 1

        for i in self.species:
            out_block_sequences = os.path.join(self.unfilter_block_dir, i + '.unevaluated.block')
            out_block_sequences = open(out_block_sequences,'w')
            block_sequences = species_block[i]
            for j in block_sequences.keys():
                if self.chr_shape == 's':
                    out_block_sequences.write('s ')
                else:
                    out_block_sequences.write('c ')

                for k in block_sequences[j]:
                    if k.startswith('-'):
                        block = k[1:]
                    else:
                        block = k
                    if block in required_block:
                        out_block_sequences.write(k+' ')
                out_block_sequences.write('\n')
            out_block_sequences.close()

        for i in self.species:
            block_synteny_file_path = os.path.join(self.manual_block_dir, i + '.manual.synteny')
            block_final_synteny_file_path = os.path.join(self.unfilter_block_dir, i + '.final.synteny')
            with open(block_synteny_file_path, 'r') as bsf:
                with open(block_final_synteny_file_path, 'w') as bfsf:
                    for line in bsf:
                        line = line.strip()
                        itemset = line.split(' ')
                        block = itemset[0].split(':')[0]
                        if block in required_block:
                            for item in itemset:
                                bfsf.write(item + ' ')
                            bfsf.write('\n')

        for i in self.species:
            block_synteny_file_path = os.path.join(self.manual_block_dir, i + '.manual.synteny.genename')
            block_final_synteny_file_path = os.path.join(self.unfilter_block_dir, i + '.final.synteny.genename')
            with open(block_synteny_file_path, 'r') as bsf:
                with open(block_final_synteny_file_path, 'w') as bfsf:
                    for line in bsf:
                        line = line.strip()
                        itemset = line.split(' ')
                        block = itemset[0].split(':')[0]
                        if block in required_block:
                            for item in itemset:
                                bfsf.write(item + ' ')
                            bfsf.write('\n')


        for i in self.species:
            block_synteny_count = self._syntenyCount(self.synteny_dict)
            filter_synteny_count = {}
            syntenygenes_number = 0
            for j in block_synteny_count.keys():
                if j not in required_block:
                    filter_synteny_count[j] = block_synteny_count[j]
                    syntenygenes_number += filter_synteny_count[j]
            allgenes_number = 0
            for j in species_all_sequences[i].keys():
                allgenes_number += len(species_all_sequences[i][j])
            print(i)
            print(syntenygenes_number / allgenes_number)
            print(allgenes_number)
