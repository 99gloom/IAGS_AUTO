class FindBackProcess():
    def __init__(self, graph_dir_path, sp_drimm_dir_path, sp_name):
        B2M_file_name = 'root-B2M.txt'
        M2B_file_name = 'root-M2B.txt'
        virtual_file_name = 'virtual.txt'
        B2M_path = graph_dir_path + B2M_file_name
        M2B_path = graph_dir_path + M2B_file_name
        virtual_path = graph_dir_path + virtual_file_name

        virtual, B2M, M2B = self._read_manual_file(virtual_path, B2M_path, M2B_path)
        sp_block_info, sp_block_signal_info = self._read_block_file(sp_drimm_dir_path, sp_name)

        # 1. 首先对实际的block进行评级,通过block的权重判别manual的方向,并获得manual中block的相对顺序，再处理manual内容
        # block_priority_no_virtual: manual: ['b1', 'b2', ...]
        block_priority_no_virtual = self._count_inner_block(B2M, M2B, sp_block_info)
        # 处理物种原始block
        changed_blocks, manual_info = self._find_block_no_virtual(sp_block_info, sp_block_signal_info, B2M,
                                                                  block_priority_no_virtual)


        # 2. virtual进行额外处理，统计M的内部block出现的次数来确定方向, 插入virtual, 还原 M
        self._find_block_virtual(changed_blocks, manual_info, virtual_path, B2M)

        print(changed_blocks)

    def _strip_block(self, s):
        return s if not s.startswith('-') else s[1:]

    def reverse_block_order(self, block: str):
        if block.startswith('-'):
            return block[1:]
        else:
            return '-' + block

    def _read_block_file(self, sp_dir_path, sp_file_name):
        # 读取物种原始block,并生成对应M矩阵
        sp_block_info = {}
        sp_block_signal_info = {}
        for sp in sp_file_name:
            sp_block_info[sp] = []
            sp_block_signal_info[sp] = []
            with open(sp_dir_path + sp + '.block', 'r') as f:
                for row in f:
                    row_info = row.strip().split(' ')[1:]
                    sp_block_info[sp].append(row_info)
                    sp_block_signal_info[sp].append([''] * len(row_info))
        return sp_block_info, sp_block_signal_info

    def _read_manual_file(self, virtual_path, B2M_path, M2B_path):
        # 读取上一步生成的所有M信息

        # virtual:
        # virtual_block: [(adj_1, adj_2), sp_1, sp_2, ... ]
        virtual = {}
        with open(virtual_path, 'r') as v_f:
            for row in v_f:
                row_info = row.strip().split(':')
                virtual_block = row_info[0]
                adj = tuple(sorted(row_info[1].split(' ')))
                sp = row_info[2].split(' ')
                virtual[virtual_block] = [adj, sp]

        # B2M:
        # block: manual
        B2M = {}
        with open(B2M_path, 'r') as B2M_f:
            for row in B2M_f:
                row_info = row.strip().split(':')
                block = row_info[0]
                manual = row_info[1]
                B2M[block] = manual

        # M2B:
        # manual: block
        M2B = {}
        with open(M2B_path, 'r') as M2B_f:
            for row in M2B_f:
                row_info = row.strip().split(':')
                manual = row_info[0]
                block = row_info[1]
                M2B[manual] = set(block.split(' '))

        return virtual_block, B2M, M2B

    def _count_inner_block(self, B2M: dict, M2B: dict, sp_block_info: dict):
        # 计算M中,每个block 出现的次数,最终以降序排列

        # 初始化统计字典  manual: block: N
        manual_block_count = {}
        for manual, blocks in M2B.items():
            manual_block_count[manual] = {}
            for block in blocks:
                manual_block_count[manual][block] = 0

        for sp, all_block in sp_block_info.items():
            for row in all_block:
                for each_block in row:
                    block = each_block if not each_block.startswith('-') else each_block[1:]
                    if block in B2M.keys():
                        manual = B2M[block]
                        manual_block_count[manual][block] += 1

        manual_block_priority = {}
        # 优先级字典  manual: [b1,b2, ... ]
        for k, v in manual_block_count.items():
            manual_block_priority[k] = [i[0] for i in sorted(v.items(), key=lambda x: x[1], reverse=True)]
        return manual_block_priority

    def _find_block_no_virtual(self, sp_block_info: dict[list], sp_block_signal_info: dict[list],
                                        B2M: dict, manual_block_priority_without_virtual):
        '''
        处理非virtual block的manual

        1.获得manual和原始block的位置以及方向关系,
        2.获取manual中block的order
        3.记录邻接,为virtual做准备
        '''
        def init_manual_info(manual, manual_info_dict: dict):
            if manual not in manual_info_dict.keys():
                manual_info_dict[manual] = {}
                manual_info_dict[manual]['adj'] = {}
                manual_info_dict[manual]['block'] = []


        def get_block_strand(original_block_order, manual_block_priority_without_virtual, last):
            manual_strand = '+'
            block_strand = {}
            for i in original_block_order:
                if i.startswith('-'):
                    block_strand[self._strip_block(i)] = '-'
                else:
                    block_strand[i] = '+'
            # 获取manual的strand，选取优先级最高的block作为manual的方向
            for i in manual_block_priority_without_virtual[last]:
                if i in block_strand:
                    manual_strand = block_strand[i]
                    break
            return manual_strand
        def get_adj_strand(start, end, original_block_row) :
            adj_flag = True
            if start <= 0:
                print("left end")
                adj_flag = False
            if end >= len(original_block_row):
                print("right end")
                adj_flag = False

            if adj_flag:
                return [original_block_row[start-1], original_block_row[end+1]]
            else:
                return None
        def merge_ordered_lists(list1, list2):
            result = []
            i, j = 0, 0

            while i < len(list1) and j < len(list2):
                if list1[i] < list2[j]:
                    result.append(list1[i])
                    i += 1
                elif list1[i] > list2[j]:
                    result.append(list2[j])
                    j += 1
                else:  # 当元素相等时，只添加一个，并移动两个指针
                    result.append(list1[i])
                    i += 1
                    j += 1

            # 添加剩余的元素
            while i < len(list1):
                result.append(list1[i])
                i += 1

            while j < len(list2):
                result.append(list2[j])
                j += 1

            return result

        # 为sp_block_signal_info赋值
        sp_list = list(sp_block_info.keys())
        for sp in sp_list:
            for row_index in range(len(sp_block_info[sp])):
                for index in range(len(sp_block_info[sp][row_index])):
                    block = self._strip_block(sp_block_info[sp][row_index][index])
                    if block in B2M.keys():
                        sp_block_signal_info[sp][row_index][index] = B2M[block]

        manual_info = {}
        changed_blocks = {}
        for sp in sp_list:
            changed_blocks[sp] = []
            for row_index in range(len(sp_block_signal_info[sp])):
                changed_blocks[sp].append([])
                # manual标志序列
                signal_block_row = sp_block_signal_info[sp][row_index]
                # 原始block序列
                original_block_row = sp_block_info[sp][row_index]
                # 更改后的block序列
                changed_block_row = changed_blocks[sp][row_index]

                # 优先通过邻接获得M方向，记录邻接信息，如果已经存储过，则以以前的邻接作为方向
                collect = []
                for index in range(1, len(signal_block_row)):
                    last = signal_block_row[index - 1]
                    if last:
                        collect.append(index - 1)
                        if last != signal_block_row[index]:
                            init_manual_info(last, manual_info)
                            start = collect[0]
                            end = collect[-1]
                            original_block_order = original_block_row[start:end+1]

                            # 根据block内部优先级评判manual方向
                            manual_block_strand_info = get_block_strand(original_block_order,
                                                                        manual_block_priority_without_virtual, last)
                            # 根据adj评判manual方向
                            manual_adj_strand_info = get_adj_strand(start, end, original_block_row)

                            if manual_adj_strand_info:
                                reversed_manual_adj_strand_info = [self.reverse_block_order(i) for i in manual_adj_strand_info[::-1]]
                                if tuple(manual_adj_strand_info) in manual_info[last]['adj']:
                                    # 寻找是否有已存在邻接作为方向
                                    manual_strand = manual_info[last]['adj'][tuple(manual_adj_strand_info)]
                                else:
                                    # 如果找不到的话则选择block信息作为方向
                                    manual_strand = manual_block_strand_info
                                    if manual_block_strand_info == '+':
                                        manual_info[last]['adj'][tuple(manual_adj_strand_info)] = '+'
                                        manual_info[last]['adj'][tuple(reversed_manual_adj_strand_info)] = '-'
                                    else:
                                        manual_info[last]['adj'][tuple(manual_adj_strand_info)] = '-'
                                        manual_info[last]['adj'][tuple(reversed_manual_adj_strand_info)] = '+'
                            else:
                                manual_strand = manual_block_strand_info

                            # 记录block顺序，默认为+
                            if manual_strand == '+':
                                manual_info[last]['block'].append(original_block_order)
                                changed_block_row.append(last)
                            else:
                                manual_info[last]['block'].append([self.reverse_block_order(i) for i in original_block_order[::-1]])
                                changed_block_row.append('-' + last)
                            collect = []
                    else:
                        changed_block_row.append(
                            original_block_row[index - 1]
                        )
                        collect = []
                last = signal_block_row[-1]
                if last:
                    collect.append(len(signal_block_row) - 1)
                    init_manual_info(last, manual_info)
                    original_block_order = original_block_row[collect[0]: collect[-1] + 1]
                    manual_strand = get_block_strand(original_block_order,
                                                                manual_block_priority_without_virtual, last)
                    if manual_strand == '+':
                        manual_info[last]['block'].append(original_block_order)
                        changed_block_row.append(last)
                    else:
                        manual_info[last]['block'].append([self.reverse_block_order(i) for i in original_block_order[::-1]])
                        changed_block_row.append('-' + last)
                else:
                    changed_block_row.append(
                        original_block_row[- 1]
                    )

        # 合并manual_info中block顺序
        for m in manual_info.keys():
            res = manual_info[m]['block'][0]
            for block_list in manual_info[m]['block'][1:]:
                res = merge_ordered_lists(res, block_list)
            manual_info[m]['block'] = res

        return changed_blocks, manual_info

    def _find_block_virtual(self, changed_blocks, manual_info, virtual_path, B2M):
        with open(virtual_path, 'r') as f:
            for line in f:
                line_info = line.strip().split(':')
                virtual_name = line_info[0]
                adj = line_info[1].split(' ')
                sp_list = line_info[2].split(' ')

                for sp in sp_list:
                    res = []
                    for row in range(len(changed_blocks[sp])):
                        res.append([])
                        for col in range(len(changed_blocks[sp][row])-1):
                            left_block = changed_blocks[sp][row][col]
                            right_block = changed_blocks[sp][row][col + 1]
                            unsigned_left_block = self._strip_block(left_block)
                            unsigned_right_block = self._strip_block(right_block)
                            res[row].append(left_block)

                            if [unsigned_left_block, unsigned_right_block] == adj \
                                    or [unsigned_right_block, unsigned_left_block] == adj:
                                manual_block = B2M[virtual_name]
                                # 判断方向，已存在的邻接关系直接取值，否则默认为+并记录
                                if (left_block, right_block) in manual_info[manual_block]['adj']:
                                    manual_signal = manual_info[manual_block]['adj'][(left_block, right_block)]
                                else:
                                    manual_signal = '+'
                                    manual_info[manual_block]['adj'][(left_block, right_block)] = '+'
                                    manual_info[manual_block]['adj'][
                                        (self.reverse_block_order(right_block),
                                         self.reverse_block_order(left_block))
                                    ] = '-'

                                # virtual位置填充manual
                                if manual_signal == '+':
                                    res[row].append(manual_block)
                                else:
                                    res[row].append('-' + manual_block)
                        res[row].append(changed_blocks[sp][row][-1])
                    changed_blocks[sp] = res




if __name__ == '__main__':
    graph_dir_path = "./graph_process_res/"
    sp_name = [
        'Brachy',
        'Sorghum',
        'Maize',
        'Rice',
        'Telongatum'
    ]
    sp_drimm_dir_path = './example/drimmBlocks/'

    finder = FindBackProcess(graph_dir_path, sp_drimm_dir_path, sp_name)
