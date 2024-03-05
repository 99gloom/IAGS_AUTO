import os


def mapping(mapping_dict: dict, input_dir: str, output_dir: str):
    # mapping_dict 为manual与模拟block的映射
    file_list = os.listdir(input_dir)
    final_file = [i for i in file_list if i.endswith('.final.block')]
    genenumber_file = ''
    for i in file_list:
        if i.endswith('genenumber'):
            genenumber_file = i

    for i in final_file:
        output_info = []
        with open(os.path.join(input_dir,i), 'r') as f:
            for line in f:
                line_split_info = line.strip().split(' ')
                output_row = [line_split_info[0]]
                items = line_split_info[1:]
                for item in items:
                    if item.startswith('-'):
                        block = item[1:]
                        order = '-'
                    else:
                        block = item
                        order = '+'

                    if block in mapping_dict.keys():
                        if order == '-':
                            output_row.append('-' + mapping_dict[block])
                        else:
                            output_row.append(mapping_dict[block])
                    else:
                        output_row.append(item)
                output_info.append(output_row)

        with open(os.path.join(output_dir, i), 'w') as f:
            for row in output_info:
                for i in row:
                    f.write(i+' ')
                f.write('\n')

    with open(os.path.join(input_dir + genenumber_file), 'r') as in_dir_f:
        with open(os.path.join(output_dir + genenumber_file), 'w') as out_dir_f:
            for line in in_dir_f:
                line_split_info = line.strip().split('\t')
                block = line_split_info[0]
                length = line_split_info[1]
                if block not in mapping_dict.keys():
                    out_dir_f.write(line)
                else:
                    out_dir_f.write(mapping_dict[block] + '\t' +  length + '\n')

    with open(os.path.join(output_dir, 'mapping_dict.txt'), 'w') as f:
        for key,val in mapping_dict.items():
            f.write(key + '\t' + val + '\n')



def mapping_back(mapping_dir:str, IAGS_res_dir:str):
    mapping_back_dir = {}
    with open(os.path.join(mapping_dir, 'mapping_dict.txt'), 'r') as f:
        for line in f:
            line_info = line.strip().split('\t')
            manual_block = line_info[0]
            num = line_info[1]
            mapping_back_dir[num] = manual_block

    file_list = os.listdir(IAGS_res_dir)
    anc_file = [i for i in file_list if (i.startswith('ANC') or i.startswith('PRE'))]
    for i in anc_file:
        output_info = []
        with open(os.path.join(IAGS_res_dir,i,i+'.block'), 'r') as f:
            for line in f:
                line_split_info = line.strip().split(' ')
                output_row = [line_split_info[0]]
                items = line_split_info[1:]
                for item in items:
                    if item.startswith('-'):
                        block = item[1:]
                        order = '-'
                    else:
                        block = item
                        order = '+'

                    if block in mapping_back_dir.keys():
                        if order == '-':
                            output_row.append('-' + mapping_back_dir[block])
                        else:
                            output_row.append(mapping_back_dir[block])
                    else:
                        output_row.append(item)
                output_info.append(output_row)

        with open(os.path.join(IAGS_res_dir,i,i+'.block'), 'w') as f:
            for row in output_info:
                for i in row:
                    f.write(i+' ')
                f.write('\n')

