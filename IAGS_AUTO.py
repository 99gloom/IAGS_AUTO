#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

import os
import shutil
import sys
from utils.tools.makeDirFunc import make_dir
from utils.tools.readInput import process_input
from utils.processTREE import EvolutionaryTree as CT
from utils.processDRIMM import processDrimmAndData as pD
from utils.processDRIMM import processEmptyChrAndBlockLength as pEmpty
from utils.processDRIMM import processLCSAndFirstFilter as pLCS
from utils.processDRIMM import processOrthoFind as pO
from utils.processDRIMM import processFinalFilter as pFinal
from utils.processIAGS import processIAGS as pIAGS
from utils.makeProfile.makeProfile import check_copy_number_validity
from utils.makeProfile.plot import process_plot


def main(argv):
    # 获取命令基本信息
    fileDir, \
    orthogroups_path, \
    tree_file_path, \
    gff_path_list, \
    cycleLength, \
    dustLength, \
    chr_shape, \
    manual_option = process_input(argv)

    pre_manual = True
    after_manual = True
    if manual_option == 'manual':
        after_manual = False
    elif manual_option == 'continue':
        pre_manual = False

    # 生成所有需要的文件目录

    '''
    1. 总目录Result，包含三个文件：
    Tree_File   Process_Drimm   IAGS
    
    2. 树文件目录Tree_File，包含三个文件：
    Evolutionary_tree.txt   species.ratio   model_and_outgroup.txt
    
    3. 生成block文件的目录processDrimm，包含四个文件夹：
    Process_OrthoFind   Drimm_Synteny_Output    Drimm_Blocks    Unfiltered_Empty_Chr_Blocks    Final_Blocks
    
    '''
    Result_OutputDir = make_dir(os.path.join(fileDir + 'Result'))

    out_tree_Dir = make_dir(os.path.join(Result_OutputDir + 'Tree_File'))
    out_process_Drimm = make_dir(os.path.join(Result_OutputDir + 'Process_Drimm'))
    out_IAGS = make_dir(os.path.join(Result_OutputDir + 'IAGS'))

    outputPath_process_orthofind = make_dir(os.path.join(out_process_Drimm + 'Process_OrthoFind'))
    outputPath_DRIMM_Synteny_Files = make_dir(os.path.join(out_process_Drimm + "Drimm_Synteny_Output"))
    outputPath_drimmBlocks = make_dir(os.path.join(out_process_Drimm + 'Drimm_Blocks'))
    outputPath_unfiltered_empty_chr_Blocks = make_dir(os.path.join(out_process_Drimm + 'Unfiltered_Empty_Chr_Blocks'))
    outputPath_finalBlocks = make_dir(os.path.join(out_process_Drimm + 'Final_Blocks'))
    outputPath_temp_blocks = make_dir(os.path.join(out_process_Drimm + 'Filter_temp_file'))


    # manual_option 手动添加block或者更换外族判断

    # 初始化树
    evolutionary_tree = CT.Evolutionary_tree(tree_file_path)

    if pre_manual:
        # process tree
        evolutionary_tree.output_all_file(out_tree_Dir)

        # process orthofinder
        processOrthoFind = pO.processOrthoFind(orthogroups_path,
                                               gff_path_list,
                                               outputPath_process_orthofind,
                                               out_tree_Dir + 'species.ratio')

        # check copy number
        check_copy_number_validity(processOrthoFind.species_ratio, processOrthoFind.rate_dir, Result_OutputDir)


        # if user doesn't input dustLength, using default number
        if not dustLength:
            dustLength = str(processOrthoFind.defult_dustLength)

        # make a profile for user to evaluate
        if manual_option == 'profile':
            out_prfile = make_dir(os.path.join(Result_OutputDir, 'Profile'))
            process_plot(processOrthoFind.sp, outputPath_process_orthofind, fileDir, out_prfile)
            exit()


        # exit()
        # process DRIMM
        processDRIMM = pD.processDrimm(outputPath_process_orthofind + "sample.sequence",
                                       outputPath_DRIMM_Synteny_Files,
                                       cycleLength,
                                       dustLength,
                                       processOrthoFind.speciesAndChrLen,
                                       processOrthoFind.sp,
                                       outputPath_drimmBlocks,
                                       chr_shape)

        # process LCS and fist filter
        processLCS = pLCS.processLCSAndFirstFilter(outputPath_drimmBlocks, outputPath_temp_blocks,
                                                   out_tree_Dir + 'species.ratio', outputPath_process_orthofind,
                                                   outputPath_DRIMM_Synteny_Files + "Synteny.txt",
                                                   processOrthoFind.sp, chr_shape)

        # process final filter and get block length
        processFinalFilter = pFinal.processFinalFilter(processOrthoFind.sp,
                                                       outputPath_temp_blocks,
                                                       outputPath_drimmBlocks,
                                                       outputPath_unfiltered_empty_chr_Blocks,
                                                       chr_shape)
        shutil.rmtree(outputPath_temp_blocks)

        # remove empty chromosome
        pEmpty.evaluateBlocks(processOrthoFind.sp, outputPath_unfiltered_empty_chr_Blocks, outputPath_finalBlocks)

    if after_manual:
        # process IAGS
        processIAGS = pIAGS.ProcessIAGS(out_tree_Dir, outputPath_finalBlocks, out_IAGS, evolutionary_tree)
        processIAGS.process_ancestor_block_and_evaluate()
        processIAGS.process_painting()
        processIAGS.process_Calculating_Fissions_and_Fusions()
        shutil.rmtree(processIAGS.multiplied_file_dir)




if __name__ == '__main__':
    main(sys.argv[1:])
