import getopt
import sys
import os
import argparse

def process_input(parsers: argparse.ArgumentParser):
    args = parsers.parse_args()
    fileDir = args.filepath
    cycleLength = args.cycleLength
    dustLength = args.dustLength
    shape = args.shape
    model = args.model
    dotplot = args.dotplot
    expand = args.expand
    check = args.check
    OrthogroupsPath = ''
    treePath = ''
    gff_path_list = []

    # 将路径格式修改正确
    while fileDir.find("//") != -1:
        fileDir = fileDir.replace('//', '/')
    if fileDir[-1] != '/':
        fileDir += '/'

    for root, dirs, files in os.walk(fileDir):
        for file in files:
            if (file.split('.')[-1] == "tsv"):
                OrthogroupsPath = fileDir + file
            if (file.split('.')[-1] == 'gff'):
                gff_path_list.append(file)
            if (file.split('.')[-1] == 'tree'):
                treePath = fileDir + file
        break
    if (OrthogroupsPath == ''):
        print("[ERROR]You need to include the TSV file in the -f file path\n")
        sys.exit()
    if (treePath == ''):
        print("[ERROR]You need to include the tree file in the -f file path\n")
        sys.exit()
    if (len(gff_path_list) == 0):
        print("[ERROR]You need to include the gff file in the -f file path\n")
        sys.exit()
    gff_path_list.sort()
    for i in range(len(gff_path_list)):
        temp = fileDir + gff_path_list[i]
        gff_path_list[i] = temp

    if check == None or check == 'yes':
        check = True
    else:
        check = False

    return fileDir, OrthogroupsPath, treePath, gff_path_list, cycleLength, dustLength, shape, model, dotplot, expand, check
