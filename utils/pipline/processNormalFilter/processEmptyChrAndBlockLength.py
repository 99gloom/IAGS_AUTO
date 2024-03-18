import os
import shutil


def evaluateBlocks(sp, unevaluated_dir, outdir, check):
    flag = False
    for i in sp:
        remain_row = []
        empty_row_number = []
        empty_value_Threshold = 0.3
        row_number = 1
        with open(os.path.join(unevaluated_dir, i + '.unevaluated.block'), 'r') as blockFile:
            for row in blockFile:
                row_info = row.rstrip('\n').rstrip(' ').split(' ')
                if len(row_info) == 1:
                    empty_row_number.append(row_number)
                else:
                    remain_row.append(row)
                row_number += 1

        if len(empty_row_number) / (row_number - 1) > empty_value_Threshold:
            flag = True

        with open(os.path.join(outdir, i + '.final.block'), 'w') as finalBlockFile:
            for row in remain_row:
                finalBlockFile.write(row)

        if empty_row_number:
            with open(os.path.join(outdir, i + '_removed_chr.log'), 'w') as log:
                log.write('Remove empty chromosome after filter:\n')
                for num in empty_row_number:
                    log.write(str(num) + ' ')

    file_list = os.listdir(unevaluated_dir)
    for i in file_list:
        if 'synteny' in i or 'genenumber' in i:
            shutil.copy(os.path.join(unevaluated_dir, i), os.path.join(outdir, i))

    if flag and check:
        print('Program Interruption: too many empty chromosome, please adjust the cycleThreshold parameter!\n'
              'Or you can select the parameter "--check no" to ignore the check')
        exit()
