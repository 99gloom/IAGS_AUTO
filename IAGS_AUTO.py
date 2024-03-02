#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

import os
import shutil
import argparse
from utils.tools.makeDirFunc import make_dir
from utils.tools.readInput import process_input
from utils.pipline.processTREE import EvolutionaryTree as CT
from utils.pipline.processNormalFilter import processEmptyChrAndBlockLength as pEmpty
from utils.pipline.processNormalFilter import processLCSAndFirstFilter as pLCS, processFinalFilter as pFinal
from utils.pipline.processGraphFilter.processExpandBlockByGraph import processGraphFilter
from utils.pipline.processDRIMM import processOrthoFind as pO, processDrimmAndData as pD
from utils.pipline.processIAGS import processIAGS as pIAGS
from utils.tools.makeDotplot.makeProfile import check_copy_number_validity
from utils.tools.makeDotplot.plot import process_plot

parsers = argparse.ArgumentParser(description='The project is based on IAGS and automates process')

def main():
    parsers.add_argument('-f','--filepath', help='filePath : Please enter the filePath.', required=True, type=str)
    parsers.add_argument('-c', '--cycleLength', required=False, type=int, default=20,
                         help=('cycleLength : Minimum cycle length in DRIMM-Synteny algorithm. The default value is 20.')
                         )
    parsers.add_argument('-d', '--dustLength', required=False, type=int, default=None,
                         help=('dustLength : Maximum genetic diversity in DRIMM-Synteny algorithm. '
                               'The default value is all species copy number plus one.'))
    parsers.add_argument('-s', '--shape', required=False, type=str, default='s',choices=['s','c'],
                         help=('chromosomes shape: \'s\' represents a string chromosome, '
                               'and \'c\' represents a circular chromosome. '
                               'The default value is \'s\'.'))
    parsers.add_argument('-m', '--model', required=False, type=str, default=None,choices=['manual','continue'],
                         help=('model: \'manual\' is for manual modification of outgroup or block '
                               'and \'continue\' for continue running the program after manually modifying the file. '
                               'Default is to run automatically without interruption. '
                               ))
    parsers.add_argument('--dotplot', required=False, action='store_true',
                         help='\'dotplot\' is generate a two-by-two dotplot for each species')

    parsers.add_argument('-e','--expand', required=False, action='store_true',
                         help='\'expand\' is increasing synteny block coverage by graph algorithms.')


    fileDir,\
    orthogroups_path, \
    tree_file_path, \
    gff_path_list, \
    cycleLength, \
    dustLength, \
    chr_shape, \
    manual_option, \
    dotplot, \
    expand= process_input(parsers)

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

    Result_OutputDir = make_dir(os.path.join(fileDir, 'Result'))

    out_tree_Dir = make_dir(os.path.join(Result_OutputDir, 'Tree_File'))
    out_process_Drimm = make_dir(os.path.join(Result_OutputDir, 'Process_Drimm'))
    out_IAGS = make_dir(os.path.join(Result_OutputDir, 'IAGS'))

    outputPath_process_orthofind = make_dir(os.path.join(out_process_Drimm, 'Process_OrthoFind'))
    outputPath_DRIMM_Synteny_Files = make_dir(os.path.join(out_process_Drimm, "Drimm_Synteny_Output"))

    # not expand
    if not expand:
        out_filter = make_dir(os.path.join(out_process_Drimm, "Filter"))
        outputPath_drimmBlocks = make_dir(os.path.join(out_filter, 'Drimm_Blocks'))
        outputPath_unfilter_empty_chr_Blocks = make_dir(os.path.join(out_filter, 'Unfilter_Empty_Chr_Blocks'))
        outputPath_finalBlocks = make_dir(os.path.join(out_filter, 'Final_Blocks'))
        outputPath_temp_blocks = make_dir(os.path.join(out_filter, 'Filter_temp_file'))

    # expand
    else:
        out_filter = make_dir(os.path.join(out_process_Drimm, "Graph_Filter"))
        outputPath_drimmBlocks = make_dir(os.path.join(out_filter, 'Drimm_Blocks'))
        merge_info_dir  = make_dir(os.path.join(out_filter, 'merge_info_dir'))
        outputPath_unfilter_empty_chr_Blocks = make_dir(os.path.join(out_filter, 'Unfilter_Empty_Chr_Blocks'))
        outputPath_finalBlocks = make_dir(os.path.join(out_filter, 'Final_blocks'))


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

        # check part
        # copy number
        check_copy_number_validity(processOrthoFind.species_ratio, processOrthoFind.rate_dir, Result_OutputDir)
        # dustLength
        if not dustLength:
            dustLength = processOrthoFind.default_dustLength
        # dotplot
        if dotplot:
            out_dotplot = make_dir(os.path.join(Result_OutputDir, 'Dotplot'))
            process_plot(processOrthoFind.sp, outputPath_process_orthofind, fileDir, out_dotplot)
            exit()

        # DRIMM
        # processDRIMM = pD.processDrimm(outputPath_process_orthofind + "sample.sequence",
        #                                outputPath_DRIMM_Synteny_Files,
        #                                cycleLength,
        #                                dustLength,
        #                                processOrthoFind.speciesAndChrLen,
        #                                processOrthoFind.sp,
        #                                outputPath_drimmBlocks,
        #                                chr_shape)

        '''
        DRIMM后处理——更改位置
        '''
        if not expand:
            # LCS and fist filter
            pLCS.processLCSAndFirstFilter(outputPath_drimmBlocks, outputPath_temp_blocks,
                                                       os.path.join(out_tree_Dir, 'species.ratio'), outputPath_process_orthofind,
                                                       os.path.join(outputPath_DRIMM_Synteny_Files, "synteny.txt"),
                                                       processOrthoFind.sp, chr_shape)

            # final filter and get block length
            pFinal.processFinalFilter(processOrthoFind.sp,
                                      outputPath_temp_blocks,
                                      outputPath_drimmBlocks,
                                      outputPath_unfilter_empty_chr_Blocks,
                                      chr_shape)
            shutil.rmtree(outputPath_temp_blocks)

            pEmpty.evaluateBlocks(processOrthoFind.sp, outputPath_unfilter_empty_chr_Blocks, outputPath_finalBlocks)
            '''
            DRIMM后处理——更改位置
            '''
        else:
            processGraphFilter(processOrthoFind.sp,
                               os.path.join(out_tree_Dir, 'species.ratio'),
                               evolutionary_tree.parser_merger_order(),
                               outputPath_drimmBlocks, merge_info_dir,
                               outputPath_process_orthofind,
                               os.path.join(outputPath_DRIMM_Synteny_Files, "synteny.txt"),
                               outputPath_unfilter_empty_chr_Blocks,
                               chr_shape)


            # remove empty chromosome
            # pEmpty.evaluateBlocks(processOrthoFind.sp, outputPath_unfilter_empty_chr_Blocks, outputPath_finalBlocks)

    # if after_manual:
    #     # process IAGS
    #     processIAGS = pIAGS.ProcessIAGS(out_tree_Dir, outputPath_finalBlocks, out_IAGS, evolutionary_tree)
    #     processIAGS.process_ancestor_block_and_evaluate()
    #     processIAGS.process_painting()
    #     processIAGS.process_Calculating_Fissions_and_Fusions()
    #     shutil.rmtree(processIAGS.multiplied_file_dir)




if __name__ == '__main__':
    main()
