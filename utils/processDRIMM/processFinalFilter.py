from copy import deepcopy


class processFinalFilter:
    def __init__(self, sp, blockDir, syntenyDir, resultDir, chr_shape='s'):
        '''
        对LCS后的结果精细过滤

        :param sp: 物种列表
        :param blockDir: 初次过滤后的block所在文件夹
        :param syntenyDir: LCS后各个物种synteny文件的文件夹
        :param resultDir: 结果输出路径
        :param chr_shape: 染色体类型 c/s
        '''

        self.sp = sp
        self.blockDir = blockDir
        self.syntenyDir = syntenyDir
        self.resultDir = resultDir
        self.chr_shape = chr_shape

        self.__excute()

    def __outputFilterBlock(self, blockDir, syntenyDir, resultDir, sp):

        # 获取synteny所有匹配成功的block以及次数
        synteny_block_position = {}
        with open(syntenyDir + '/' + sp + '.synteny', 'r') as sf:
            for i in sf:
                block_info = i
                block_info = block_info.split(' ')[0].split(':')
                if block_info[0] not in synteny_block_position.keys():
                    synteny_block_position[block_info[0]] = []
                synteny_block_position[block_info[0]].append(int(block_info[1]))

        # 获取所有block序列
        block_sequence = []
        with open(blockDir + '/' + sp + '.block', 'r') as inf:
            for i in inf:
                block_info = i
                block_info = block_info.rstrip('\n').rstrip()
                block_info = block_info.split(' ')[1:]
                block_sequence.append(block_info)

        # 初始化block出现序列为0
        raw_block_position = {}
        for i in block_sequence:
            for j in i:
                if j.startswith('-'):
                    j = j[1:]
                if j not in raw_block_position.keys():
                    raw_block_position[j] = 0
        block_list = deepcopy(raw_block_position)

        with open(resultDir + '/' + sp + '.unevaluated.block', 'w') as outf:
            for i in block_sequence:
                if self.chr_shape == 's' or self.chr_shape == 'S':
                    outf.write('s ')
                else:
                    outf.write('c ')

                for j in i:
                    if j.startswith('-'):
                        block = j[1:]
                    else:
                        block = j

                    raw_block_position[block] += 1
                    if raw_block_position[block] in synteny_block_position[block]:
                        outf.write(j + ' ')
                outf.write('\n')
        return block_list

    def __outputFilterSynteny(self, block_list, syntenyDir, resultDir, sp):
        block_list_copy = deepcopy(block_list)
        with open(resultDir + '/' + sp + '.final.synteny', 'w') as synWriteFile:
            with open(syntenyDir + '/' + sp + '.synteny', 'r') as synReadFile:
                for i in synReadFile:
                    temp = i.rstrip('\n').rstrip()
                    temp = temp.split(' ')
                    block = temp[0].split(':')[0]
                    if block in block_list.keys():
                        block_list[block] += 1
                        synWriteFile.write(block + ':' + str(block_list[block]) + ':' + temp[0].split(':')[2] + ':' +
                                           temp[0].split(':')[3] + ' ')
                        synWriteFile.write(' '.join(temp[1:]) + ' \n')

        with open(resultDir + '/' + sp + '.final.synteny.genename', 'w') as synWriteFile:
            with open(syntenyDir + '/' + sp + '.synteny.genename', 'r') as synReadFile:
                for i in synReadFile:
                    temp = i.rstrip('\n').rstrip()
                    temp = temp.split(' ')
                    block = temp[0].split(':')[0]
                    if block in block_list_copy.keys():
                        block_list_copy[block] += 1
                        synWriteFile.write(
                            block + ':' + str(block_list_copy[block]) + ':' + temp[0].split(':')[2] + ':' +
                            temp[0].split(':')[3] + ' ')
                        synWriteFile.write(' '.join(temp[1:]) + ' \n')

        return

    def processGenenumber(self, sp, resultDir):
        block_len = {}
        for i in sp:
            with open(resultDir + '/' + i + '.final.synteny', 'r') as sf:
                for line in sf:
                    temp = line
                    temp = temp.rstrip('\n').rstrip()
                    block = temp.split(' ')[0].split(':')[0]
                    length = len(temp.split(' ')[1:])
                    # print(temp.split(' ')[1:])

                    if block not in block_len.keys():
                        block_len[block] = length
                    else:
                        if block_len[block] < length:
                            # print(block)
                            # print(length)
                            # print(block_len[block])
                            block_len[block] = length
                            # print(block_len[block])
                            # print('-----------------')

        with open(resultDir + '/blockindex.genenumber', 'w') as f:
            f.write('blockID\tblockLength\n')
            for i, j in block_len.items():
                f.write(i + '\t' + str(j) + '\n')
        return

    def __excute(self):

        for i in self.sp:
            temp = self.__outputFilterBlock(self.blockDir, self.syntenyDir, self.resultDir, i)
            self.__outputFilterSynteny(temp, self.syntenyDir, self.resultDir, i)

        self.processGenenumber(self.sp, self.resultDir)
