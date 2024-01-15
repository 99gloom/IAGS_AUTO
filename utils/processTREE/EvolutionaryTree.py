from treelib import Tree, Node


class Evolutionary_tree():
    def __init__(self, path: str):
        self.tree, \
        self.computing_list, \
        self.leaves_nodes_id, \
        self.non_leaves_node_id, \
        self.painting_anc_node \
            = self.__calculate_tree(path)
        self.all_nodes_id_without_root = self.leaves_nodes_id + self.non_leaves_node_id
        self.all_nodes_id_without_root.remove('root')
        self.all_nodes_id_without_root.sort()

    def __calculate_tree(self, tree_path):
        tree_format = self.__get_tree_file_info(tree_path)
        tree = self.__create_tree(tree_format)
        tree.show()
        tree.show(data_property='real')
        # 得到节点信息并输出所有的构建信息
        leaves_nodes_id, non_leaves_node_id, computing_nodes_id, painting_node_id = self.__get_nodes_info(tree)
        computing_list = self.__calculate_ancestor_nodes(tree, leaves_nodes_id, non_leaves_node_id,
                                                         computing_nodes_id)

        return tree, computing_list, leaves_nodes_id, non_leaves_node_id, painting_node_id

    def output_all_file(self, dir):
        self.__save_tree_info(dir)
        self.__save_species_ratio(dir)
        self.__save_computing_order(dir)
        self.__save_all_ratio(dir)
        self.__save_root_info(dir)

    def __save_root_info(self, dir):
        if dir[-1] != '/':
            dir = dir + '/'
        with open(dir + 'root_file_for_graph_iter.txt', 'w') as gf:
            gf.write(self.tree.root + ':' + ':'.join([i.identifier for i in self.tree.children('root')]))

    def __save_tree_info(self, dir):
        if dir[-1] != '/':
            dir = dir + '/'
        with open(dir + 'Evolutionary_tree.txt', 'w'):
            pass
        self.tree.save2file(dir + 'Evolutionary_tree.txt')
        return

    def __save_species_ratio(self, dir):
        if dir[-1] != '/':
            dir = dir + '/'
        with open(dir + 'species.ratio', 'w') as srf:
            for i in self.leaves_nodes_id:
                srf.write(i + '\t' + str(self.tree.get_node(i).data) + '\n')

    def __save_all_ratio(self, dir):
        if dir[-1] != '/':
            dir = dir + '/'
        with open(dir + 'all.ratio', 'w') as srf:
            for i in self.all_nodes_id_without_root:
                srf.write(i + '\t' + str(self.tree.get_node(i).data) + '\n')

    def __save_computing_order(self, dir):
        if dir[-1] != '/':
            dir = dir + '/'
        with open(dir + 'model_and_outgroup.txt', 'w') as cof:
            for i in self.computing_list:
                cof.write(i + '\n')

    def __get_tree_file_info(self, tree_file_path):
        with open(tree_file_path) as tf:
            for line in tf:
                tree_info = line.rstrip('\n').rstrip(' ')

        # 将树文件分开，只保留括号及其物种
        tree_format = []
        str_temp = ''
        for i in tree_info:
            if i == ',' or i == '[' or i == ']':
                if str_temp:
                    tree_format.append(str_temp)
                    str_temp = ''
                continue
            elif i == '(' or i == ')':
                if str_temp:
                    tree_format.append(str_temp)
                    str_temp = ''
                tree_format.append(i)
            else:
                str_temp += i
        # print(tree_format)
        return tree_format

    def __create_tree(self, tree_format):
        tree = Tree()
        current_ancestor_node_point = Node()
        last_ancestor_node_point = Node()
        ANC_number = 1
        ANC_prefix = 'ANC'
        WGD_prefix = 'PRE-'
        ANC = ANC_prefix + str(ANC_number)
        count = 0
        for i in tree_format:
            if count == 0:
                # 创建根节点
                current_ancestor_node_point = tree.create_node(tag='root', identifier='root', parent=None, data=0)
            else:
                if i == '(':
                    # 遇到'('向下增加需要计算的添加祖先节点
                    last_ancestor_node_point = current_ancestor_node_point
                    current_ancestor_node_point = tree.create_node(tag=ANC, identifier=ANC,
                                                                   parent=current_ancestor_node_point, data=1)
                    ANC_number += 1
                    ANC_prefix = 'ANC'
                    ANC = ANC_prefix + str(ANC_number)
                elif i == ')':
                    # 遇到')'向上回退上祖先节点
                    last_ancestor_node_point = current_ancestor_node_point
                    current_ancestor_node_point = tree.parent(current_ancestor_node_point.identifier)
                elif i == 'WGD':
                    if last_ancestor_node_point.identifier == 'root':
                        temp_tree = tree.subtree(last_ancestor_node_point.identifier)
                        for node in temp_tree.all_nodes():
                            node.data = node.data * 2
                        continue
                    # 遇到WGD创建新节点并且更改子树拷贝数
                    WGD_node_id = WGD_prefix + last_ancestor_node_point.identifier
                    WGD_node = tree.create_node(tag=WGD_node_id, identifier=WGD_node_id,
                                                parent=current_ancestor_node_point,
                                                data=1)
                    tree.move_node(last_ancestor_node_point.identifier, WGD_node.identifier)
                    temp_tree = tree.subtree(last_ancestor_node_point.identifier)
                    for node in temp_tree.all_nodes():
                        node.data = node.data * 2
                    last_ancestor_node_point = WGD_node
                else:
                    # 如果是物种则生成一个叶子节点并更新last指针
                    last_ancestor_node_point = tree.create_node(tag=i, identifier=i, parent=current_ancestor_node_point,
                                                                data=1)

            count += 1
        return tree

    def __get_nodes_info(self, tree):
        leaves_nodes_id = []
        non_leaves_node_id = []
        computing_nodes_id = []
        # 所有叶子节点
        for node in tree.leaves():
            leaves_nodes_id.append(node.identifier)
        # 所有非叶子节点
        for node in tree.all_nodes():
            if node.identifier not in leaves_nodes_id:
                non_leaves_node_id.append(node.identifier)

        # 计算节点，从叶子节点反推一层，孩子节点都是叶子结点的为计算节点
        for node_id in leaves_nodes_id:
            parent_id = tree.parent(node_id).identifier
            children_id = []
            flag = 1

            # 获取所有的子节点
            for i in tree.children(parent_id):
                children_id.append(i.identifier)

            # 判断子节点是否全都是叶子节点
            for i in children_id:
                if i not in leaves_nodes_id:
                    flag = 0

            # 得到所有满足条件的节点，并将以深度排序
            if parent_id not in computing_nodes_id and 1 == flag:
                computing_nodes_id.append(parent_id)
            if 'root' in computing_nodes_id:
                computing_nodes_id.remove('root')
            computing_nodes_id = self.__sort_nodes_by_level(tree, computing_nodes_id)

        painting_node_id = computing_nodes_id[0]
        while (tree.parent(painting_node_id).identifier != 'root'):
            painting_node_id = tree.parent(painting_node_id).identifier

        leaves_nodes_id.sort()
        non_leaves_node_id.sort()
        return leaves_nodes_id, non_leaves_node_id, computing_nodes_id, painting_node_id

    def __calculate_ancestor_nodes(self, tree, leaves_nodes_id, non_leaves_node_id, computing_nodes_id):
        computed_nodes_id = []
        max_step = 3 * (tree.level(computing_nodes_id[0]) + 1)
        outgroup_step_length = 2
        modes_list = []
        # 自底向上、由近到远计算所有的节点
        while (computing_nodes_id and outgroup_step_length < max_step):
            computing_number = 0
            max_computing_number = len(computing_nodes_id) - 1

            # 从computing_nodes_id列表中获取重建进化树的先后顺序
            while (computing_number <= max_computing_number):
                each_mode = ''
                parent_computing_node = tree.get_node(computing_nodes_id[computing_number])
                children_nodes = tree.children(parent_computing_node.identifier)
                leaves_child_node = []
                non_leaves_child_node = []
                for i in children_nodes:
                    if i.identifier in leaves_nodes_id:
                        leaves_child_node.append(i)
                    elif i.identifier in non_leaves_node_id:
                        non_leaves_child_node.append(i)
                    else:
                        print('error children')
                children_nodes = leaves_child_node + non_leaves_child_node

                # GMP和MultiCopyGMP为类似的构建模型，GGHP和MultiCopyGGHP为类似的构建模型
                # GMP和MultiCopyGMP
                if parent_computing_node.data == children_nodes[0].data:
                    if parent_computing_node.data == 1:
                        each_mode += 'GMP'
                    else:
                        each_mode += 'MultiCopyGMP'
                    outgroup_id = self.__find_outgroup(tree, parent_computing_node.identifier, outgroup_step_length,
                                                       leaves_nodes_id, computed_nodes_id)
                    if outgroup_id:
                        outgroup_times = ''
                        if tree.get_node(outgroup_id).data != parent_computing_node.data:
                            outgroup_times += '*' + str(int(parent_computing_node.data/tree.get_node(outgroup_id).data))
                        each_mode = parent_computing_node.identifier + ':' \
                                    + each_mode + ':' \
                                    + children_nodes[0].identifier + ':' \
                                    + children_nodes[1].identifier + ':' \
                                    + outgroup_id + outgroup_times
                        modes_list.append(each_mode)
                        computed_nodes_id.append(parent_computing_node.identifier)
                        computing_nodes_id[computing_number] = tree.parent(parent_computing_node.identifier).identifier
                        if computing_nodes_id[computing_number] == 'root':
                            max_computing_number -= 1
                            computing_nodes_id.pop(computing_number)
                    else:
                        computing_number += 1
                # GGHP和MultiCopyGGHP
                elif children_nodes[0].data / parent_computing_node.data == 2:
                    if parent_computing_node.data == 1:
                        each_mode += 'GGHP'
                    else:
                        each_mode += 'MultiCopyGGHP'
                    outgroup_id = self.__find_outgroup(tree, parent_computing_node.identifier, outgroup_step_length,
                                                       leaves_nodes_id, computed_nodes_id)
                    if outgroup_id:
                        outgroup_times = ''
                        if tree.get_node(outgroup_id).data != parent_computing_node.data:
                            outgroup_times += '*' + str(int(parent_computing_node.data / tree.get_node(outgroup_id).data))

                        each_mode = parent_computing_node.identifier + ':' \
                                    + each_mode + ':' \
                                    + children_nodes[0].identifier + ':' \
                                    + outgroup_id + outgroup_times
                        modes_list.append(each_mode)
                        computed_nodes_id.append(parent_computing_node.identifier)
                        computing_nodes_id[computing_number] = tree.parent(parent_computing_node.identifier).identifier
                        if computing_nodes_id[computing_number] == 'root':
                            max_computing_number -= 1
                            computing_nodes_id.pop(computing_number)
                    else:
                        computing_number += 1
                else:
                    print('Wrong Copy Number')
                    exit()
            computing_nodes_id = self.__sort_nodes_by_level(tree, computing_nodes_id)
            outgroup_step_length += 1

        # 检查是否计算完全
        computed_nodes_id.append('root')
        computed_nodes_id.sort()
        non_leaves_node_id.sort()
        if computed_nodes_id != non_leaves_node_id:
            print(computed_nodes_id)
            print(non_leaves_node_id)
            print('error')
            # exit()

        return modes_list

    def __find_outgroup(self, tree, node_id, outgroup_step_length, leaves_nodes_id, computed_nodes_id):
        up_step = 1
        remain_down_step = 0
        ancestor_node_id = tree.get_node(node_id).identifier

        while (up_step < outgroup_step_length):
            # 向上爬一层计算节点，节点爬升次数增加一步
            proto_child_node_id = ancestor_node_id
            ancestor_node_id = tree.parent(ancestor_node_id).identifier
            remain_down_step = outgroup_step_length - up_step
            acquired_leaves_nodes_id = []
            acquired_computed_nodes_id = []

            # 初始化孩子节点队列，将原本的节点，可寻找的孩子节点层数减少一步
            children_nodes_id_queue = []
            for node in tree.children(ancestor_node_id):
                children_nodes_id_queue.append(node.identifier)
            children_nodes_id_queue.remove(proto_child_node_id)
            children_nodes_id_queue.sort()

            # 层序遍历的思想，用于统计向下寻找的次数
            while (children_nodes_id_queue and remain_down_step > 0):
                temp_queue = []
                for child_node_id in children_nodes_id_queue:
                    if child_node_id in leaves_nodes_id and tree.get_node(node_id).data >= tree.get_node(child_node_id).data:
                        acquired_leaves_nodes_id.append(child_node_id)
                    elif child_node_id in computed_nodes_id and tree.get_node(node_id).data >= tree.get_node(child_node_id).data:
                        acquired_computed_nodes_id.append(child_node_id)
                    else:
                        for i in tree.children(child_node_id):
                            temp_queue.append(i.identifier)

                acquired_leaves_nodes_id.sort()
                acquired_computed_nodes_id.sort()
                if acquired_leaves_nodes_id:
                    return acquired_leaves_nodes_id[0]
                elif acquired_computed_nodes_id:
                    return acquired_computed_nodes_id[0]
                else:
                    children_nodes_id_queue = temp_queue
                remain_down_step -= 1

            up_step += 1
            if ancestor_node_id == 'root':
                break

        return None

    def __sort_nodes_by_level(self, tree, node_list_id):
        for i in range(len(node_list_id)):
            for j in range(len(node_list_id) - 1, 0, -1):
                if tree.level(node_list_id[j]) > tree.level(node_list_id[j - 1]):
                    temp = node_list_id[j - 1]
                    node_list_id[j - 1] = node_list_id[j]
                    node_list_id[j] = temp
        return node_list_id
