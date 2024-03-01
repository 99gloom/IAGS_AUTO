from copy import deepcopy
from typing import List

import networkx as nx

# from matplotlib import pyplot as plt
from .BaseSyntenyClass import BaseSyntenyClass


class SelfSyntenyGraph(BaseSyntenyClass):
    # 首先需要自匹配
    def __init__(self, sp: str, dir: str, copy_number: int):
        '''

        :param sp: 物种名列表
        :param dir: drimmBlocks文件路径
        :param copy_number: 物种拷贝数
        :return:
        '''
        super().__init__()
        self.sp = {sp}
        self.cpn[sp] = copy_number
        sp_row = []
        path = f'{dir + sp}.block'
        with open(path, 'r') as f:
            for row in f:
                row_split = row.strip().split(' ')[1:]
                for j in range(len(row_split)):
                    # 去除负号
                    if row_split[j].startswith('-'):
                        row_split[j] = row_split[j][1:]
                sp_row.append(row_split)

        # 添加节点
        for row in sp_row:
            for j in row:
                if not self.G.has_node(j):
                    self.G.add_node(j, times={sp: 1})
                else:
                    self.G.nodes[j]['times'][sp] += 1
        # 添加边
        for row in sp_row:
            for j in range(len(row) - 1):
                if not self.G.has_edge(row[j], row[j + 1]):
                    self.G.add_edge(row[j], row[j + 1], weight=1, sp_set={sp})
                else:
                    self.G[row[j]][row[j + 1]]['weight'] += 1

        copy_number_dict = {
            sp: copy_number
        }
        anchor_Points = [n for n in self.G.nodes if self.G.nodes[n]['times'] == copy_number_dict]
        self._process_node(self.G, anchor_Points)
        BaseSyntenyClass.manual_to_origin.update(self.manual_to_origin)
        BaseSyntenyClass.origin_to_manual.update(self.origin_to_manual)

    def _process_node(self, G: nx.Graph, anchor_Points: List[str]):
        visited_block = []
        for cur_Anchor_block in anchor_Points:
            node_list = []
            for i in G.neighbors(cur_Anchor_block):
                node_list.append([cur_Anchor_block, i])
            last = []
            while last != node_list:
                last = deepcopy(node_list)
                for j in range(len(node_list)):
                    if node_list[j][-1] in visited_block or node_list[j][-1] in anchor_Points:
                        continue
                    visited_block.append(node_list[j][-1])
                    neighbors = list(G.neighbors(node_list[j][-1]))
                    if len(neighbors) == 2:
                        neighbors.remove(node_list[j][-2])
                        node_list[j].append(neighbors[0])

            # 获得气泡（末端相同）
            last_blocks = {}
            for each_list in node_list:
                last_block = each_list[-1]
                if last_block not in last_blocks:
                    last_blocks[last_block] = [each_list]
                else:
                    last_blocks[last_block].append(each_list)

            for last_block, lists in last_blocks.items():
                # 长度大于4的气泡不参与统计
                if len(lists) <= 1 or max([len(i) - 2 for i in lists]) >= 5:
                    continue
                new_block_name = BaseSyntenyClass.get_manual_block()
                edge_equal_flag = True
                all_inner_block = set()
                # 遍历气泡
                sorted_list = sorted(lists, key=lambda x: x[0])
                start = lists[0][0]
                end = lists[0][-1]
                all_end_edge_weight = 0
                for i in sorted_list:
                    edge_equal_flag = self.is_equal_edge_weight(i, G)
                    all_end_edge_weight += G[i[-1]][i[-2]]['weight']
                    if len(i) != 2:
                        all_inner_block.update(i[1:-1])
                    # 为缺失型block添加虚拟block
                    else:
                        if tuple(sorted([start, end])) in BaseSyntenyClass.adj_tuple_to_virtual_block:
                            # 如果Virtual节点存在,更新虚拟节点存在于的物种信息
                            all_inner_block.add(
                                BaseSyntenyClass.adj_tuple_to_virtual_block[tuple(sorted([start, end]))])
                            virtual_node = BaseSyntenyClass.adj_tuple_to_virtual_block[tuple(sorted([start, end]))]
                            BaseSyntenyClass.virtual_block_adj_info[virtual_node][1].update(G[start][end]['sp_set'])
                        else:
                            # 如果存在缺失节点,添加Virtual节点,记录虚拟节点的物种以及的两端邻接信息
                            virtual_node = BaseSyntenyClass.get_virtual_block()
                            all_inner_block.add(virtual_node)
                            BaseSyntenyClass.virtual_block_adj_info[virtual_node] = [{start, end},
                                                                                     deepcopy(G[start][end]['sp_set'])]
                            BaseSyntenyClass.adj_tuple_to_virtual_block[tuple(sorted([start, end]))] = virtual_node
                all_times = sum(i for i in G.nodes[last_block]['times'].values())

                if all_end_edge_weight == all_times and edge_equal_flag:
                    for i in all_inner_block:
                        if i not in self.origin_to_manual.keys():
                            self.origin_to_manual[i] = [new_block_name]
                        else:
                            self.origin_to_manual[i].append(new_block_name)
                        if not i.startswith(BaseSyntenyClass.virtual_block_signal):
                            G.remove_node(i)

                    self.manual_to_origin[new_block_name] = all_inner_block
                    if G.has_edge(start, end):
                        G.remove_edge(start, end)
                    self.add_new_node(G, new_block_name, start, end, all_end_edge_weight, last_block, self.sp)

        # plt.figure(figsize=(40, 40), dpi=80)
        # pos = nx.kamada_kawai_layout(self.G)
        # nx.draw(self.G, pos, with_labels=True, width=5, node_size = 500)
        # # nodes_lable = dict((n, n+':'+str(d['times'])) for n,d in self.G.nodes(data=True))
        # # nx.draw_networkx_labels(self.G, pos,labels=nodes_lable, font_size=10)
        # edge_weights = nx.get_edge_attributes(self.G, 'weight')
        # nx.draw_networkx_edge_labels(self.G, pos, edge_labels=edge_weights, font_size=10)
        # plt.show()
        # exit()
