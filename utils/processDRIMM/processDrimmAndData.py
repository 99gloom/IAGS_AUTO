import os


class processDrimm:
    def __init__(self, inputSequencePath, output_temp, cycleLength, dustLength, speciesAndChrLen, sp, drimmBlockDir,
                 chr_shape='s'):
        '''
        以命令行形式运行DRIMM，并将结果blocks.txt划分回各个物种

        :param inputSequencePath: processOrthoFind处理后得到的sample.sequence
        :param output_temp: DRIMM输出路径
        :param cycleLength: cycleLength参数
        :param dustLength: dustLength参数
        :param speciesAndChrLen: 物种和它染色体数，dict类型
        :param sp: 物种列表
        :param drimmBlockDir: DRIMM处理后，以物种为单位划分后的block存储路径
        :param chr_shape: 染色体类型 c/s
        '''

        self.inputSequencePath = inputSequencePath
        self.output_temp = output_temp
        self.cycleLength = str(cycleLength)
        self.dustLength = str(dustLength)
        self.speciesAndChrLen = speciesAndChrLen
        self.sp = sp
        self.drimmBlockDir = drimmBlockDir
        self.chr_shape = chr_shape
        self.__execute()

    def __execute(self):
        current_path = os.path.dirname(__file__)
        drimm_path = os.path.join(current_path, 'Program.exe')
        AbsSequencePath = os.path.abspath(self.inputSequencePath)
        AbsSequencePath = AbsSequencePath.replace("\\", "/")
        AbsOutPath = os.path.abspath(self.output_temp)
        AbsOutPath = AbsOutPath.replace("\\", "/")
        cmd = "mono " + drimm_path + ' ' \
              + AbsSequencePath + " " + AbsOutPath + " " + self.cycleLength + " " + self.dustLength
        drimmPrint = os.system(cmd)
        return self.__splitBySpecies()

    def __splitBySpecies(self):
        # 将DRIMM处理得到的block通过染色体数目划分为各个文件，并添上c或者s头
        readfile = open(os.path.join(self.output_temp, "blocks.txt"), 'r')
        blocks = readfile.read().split('\n')
        readfile.close()
        speciesBlock = {}
        index = 0
        output_block_files = {}
        for speciesKey in self.sp:
            speciesBlock[speciesKey] = blocks[index:index + self.speciesAndChrLen[speciesKey]]
            index += self.speciesAndChrLen[speciesKey]
        for speciesBlockItem in speciesBlock:
            output = open(os.path.join(self.drimmBlockDir, speciesBlockItem) + ".block", 'w')
            output_block_files[speciesBlockItem] = os.path.join(self.drimmBlockDir,  speciesBlockItem + ".block")

            for s in speciesBlock[speciesBlockItem]:
                if self.chr_shape == 's' or self.chr_shape == 'S':
                    output.write('s ')
                else:
                    output.write('c ')
                output.write(s)
                output.write("\n")
            output.close()
        return output_block_files
