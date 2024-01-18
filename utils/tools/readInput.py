import getopt
import sys
import os
def printHelp():
    print('')
    print('input example:')
    print('python3 program.py -f ./dir/ -c 20')
    print('')
    print('[-f]  filePath : Please enter the filePath')
    print('[-c]  cycleLength : Minimum cycle length in DRIMM-Synteny algorithm. The default value is 20')
    print('[-d]  dustLength : Maximum genetic diversity in DRIMM-Synteny algorithm. The default value is all species copy number plus one')
    print(
        "[-s]  shape: chromosomes shape, 's' represents a string chromosome"
        " and 'c' represents a circular chromosome. The default value is 's'")
    print(
        "[-m]  manual: 'manual' is for manual modification of outgroup or block"
        " and 'continue' for continue running the program after manually modifying the file."
        " Default is to run automatically without interruption.")




def process_input(argv):
    # 填充默认值
    fileDir = './'
    OrthogroupsPath = './'
    treePath = './'
    cycleLength = "20"
    dustLength = None
    gff_path_list = []
    chr_shape = 's'
    manual_option = False

    try:
        opts, args = getopt.getopt(argv, "hf:c:d:s:m:",
                                   ['help', 'filepath=', 'cycleLength=', 'dustLength=', 'shape=', 'model='])
    except getopt.GetoptError:
        print("error")
        printHelp()
        sys.exit(2)

    if (len(opts) < 1):
        printHelp()
        sys.exit()
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            printHelp()
            sys.exit()
        elif opt in ("-f", '--filepath'):
            fileDir = arg
        elif opt in ("-c", '--cycleLength'):
            cycleLength = arg
            if not str.isdigit(cycleLength):
                printHelp()
                sys.exit()
        elif opt in ("-d", '--dustLength'):
            dustLength = arg
            if not str.isdigit(dustLength):
                printHelp()
                sys.exit()
        elif opt in ("-s", '--shape'):
            chr_shape = arg
            if chr_shape not in 'sScC':
                printHelp()
                sys.exit()
        elif opt in ('-m', '--model'):
            manual_option = arg
            if manual_option in ['manual', 'Manual', 'm', 'M']:
                manual_option = 'manual'
            elif manual_option in ['continue', 'Continue','c','C']:
                manual_option = 'continue'
            elif manual_option in ['dotplot', 'dotPlot', 'DotPlot','p', 'P']:
                manual_option = 'dotplot'
            else:
                printHelp()
                sys.exit()

    # 将路径格式修改正确
    while fileDir.find("//") != -1:
        fileDir = fileDir.replace('//', '/')
    if fileDir[-1] != '/':
        fileDir += '/'

    if (fileDir == ""):
        print("you need to input file path")
        printHelp()
        sys.exit()
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
        printHelp()
        sys.exit()
    if (treePath == ''):
        print("[ERROR]You need to include the tree file in the -f file path\n")
        printHelp()
        sys.exit()
    if (len(gff_path_list) == 0):
        print("[ERROR]You need to include the gff file in the -f file path\n")
        printHelp()
        sys.exit()
    gff_path_list.sort()
    for i in range(len(gff_path_list)):
        temp = fileDir + gff_path_list[i]
        gff_path_list[i] = temp
    return fileDir, OrthogroupsPath, treePath, gff_path_list, cycleLength, dustLength, chr_shape, manual_option
