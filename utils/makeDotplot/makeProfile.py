import os.path


def check_copy_number_validity(sp_ratio: str, ratio_dict: dict, path: str) -> None:
    all_ratio = []
    all_ratio_count = []
    for key, val in ratio_dict.items():
        all_ratio.append(key)
        all_ratio_count.append(int(val))
    for i in range(len(all_ratio)):
        change_flag = 1
        for j in range(len(all_ratio)-1, i, -1):
            if all_ratio_count[j]>all_ratio_count[j - 1]:
                change_flag = 0
                temp_number = all_ratio_count[j - 1]
                all_ratio_count[j - 1] = all_ratio_count[j]
                all_ratio_count[j] = temp_number
                temp_ratio = all_ratio[j - 1]
                all_ratio[j - 1] = all_ratio[j]
                all_ratio[j] = temp_ratio

        if change_flag:
            break

    rank = {}
    num = 0
    last = None
    for i in range(len(all_ratio)):
        if all_ratio_count[i] != last:
            num += 1
        rank[all_ratio[i]] = str(num)
        last = all_ratio_count[i]


    threshold_tolerance_flag = 0
    for i in range(1):
        if all_ratio[i] == sp_ratio:
            threshold_tolerance_flag = 1
    print(threshold_tolerance_flag)
    if sp_ratio not in ratio_dict:
        print(f'Can\'t find the copy number ratio you input!')
        exit()

    if threshold_tolerance_flag:
        write_log_info(sp_ratio, ratio_dict, all_ratio, all_ratio_count, rank, path, 'Log_CopyNumber')
    else:
        warning_message = (f'The number of genes matching the copy number is too low, '
                             f'we recommend checking the evolution tree you input!\n')
        write_log_info(sp_ratio, ratio_dict, all_ratio, all_ratio_count, rank, path, 'WARNING_CopyNumber', warning_message)


def write_log_info(sp_ratio, ratio_dict, all_ratio, all_ratio_count, rank, path, name, add_str = ''):
    with open(os.path.join(path,name)+'.txt', 'w') as f:
        f.write('Top three copy number ratio:\n')
        for i in range(3):
            f.write(f'{all_ratio[i]} -> {str(all_ratio_count[i])}\n')
        f.write(f'The copy number ratio you input: {sp_ratio} -> {str(ratio_dict[sp_ratio])},\n'
                f'ranked No. {rank[sp_ratio]}\n')
        if add_str:
            f.write(add_str)


