from copy import deepcopy
from typing import List
import networkx as nx

from .BaseSyntenyClass import BaseSyntenyClass
from .SelfSyntenyGraph import SelfSyntenyGraph


class MergedSyntenyGraph(BaseSyntenyClass):
    # 第二步需要物种之间互相匹配
    def __init__(self, sp_synteny_graph_1: SelfSyntenyGraph, sp_synteny_graph_2: SelfSyntenyGraph):
        '''

        :param sp_synteny_graph_1: 节点一实例
        :param sp_synteny_graph_2: 节点二实例
        :return:
        '''
        super().__init__()
        self.sp = set()
        self.sp.update(sp_synteny_graph_1.sp)
        self.sp.update(sp_synteny_graph_2.sp)

        self.cpn = self._cat_dict(sp_synteny_graph_1.cpn, sp_synteny_graph_2.cpn)
        # 获得所有重合的block
        connect_block = self._get_filtered_manual_block(sp_synteny_graph_1.manual_to_origin,
                                                        sp_synteny_graph_2.manual_to_origin)

        # 合并同类manual block
        new_to_old = {}
        for each_list in connect_block:
            new_block_name = BaseSyntenyClass.get_manual_block()
            block_content = set()
            new_to_old[new_block_name] = deepcopy(each_list)
            for i in each_list:
                if i in sp_synteny_graph_1.manual_to_origin:
                    block_content.update(sp_synteny_graph_1.manual_to_origin[i])
                if i in sp_synteny_graph_2.manual_to_origin:
                    block_content.update(sp_synteny_graph_2.manual_to_origin[i])
            self.manual_to_origin[new_block_name] = block_content
            for i in block_content:
                if i not in self.origin_to_manual:
                    self.origin_to_manual[i] = [new_block_name]
                else:
                    self.origin_to_manual[i].append(new_block_name)

        BaseSyntenyClass.manual_to_origin.update(self.manual_to_origin)
        BaseSyntenyClass.origin_to_manual.update(self.origin_to_manual)

        # 同类节点更名
        self.G: nx.Graph = self._relabel_and_update_attributes(sp_synteny_graph_1.G, new_to_old)
        G_temp: nx.Graph = self._relabel_and_update_attributes(sp_synteny_graph_2.G, new_to_old)

        # 将第二张图的内容，添加到self图
        for node in G_temp.nodes():
            if not self.G.has_node(node):
                self.G.add_node(node, times=G_temp.nodes[node]['times'])
            else:
                self.G.nodes[node]['times'] = self._cat_dict(self.G.nodes[node]['times'], G_temp.nodes[node]['times'])

        for edge in G_temp.edges():
            begin, end = edge
            if not self.G.has_edge(begin, end):
                self.G.add_edge(begin, end,
                                weight=G_temp[begin][end]['weight'],
                                sp_set=deepcopy(G_temp[begin][end]['sp_set']))
            else:
                self.G[begin][end]['weight'] += G_temp[begin][end]['weight']
                self.G[begin][end]['sp_set'].update(G_temp[begin][end]['sp_set'])

        anchor_Points = [n for n in self.G.nodes if self.G.nodes[n]['times'] == self.cpn]
        self._process_node(self.G, anchor_Points)
        # 更新总表
        BaseSyntenyClass.manual_to_origin.update(self.manual_to_origin)
        BaseSyntenyClass.origin_to_manual.update(self.origin_to_manual)

        # done: 添加 V 节点，记录左右锚点（左右有M的话，还原的时候一样能定位）
        # 输出邻居节点

    def _dfs(self, graph, node, visited, component):
        visited[node] = True
        component.append(node)

        for neighbor, connected in enumerate(graph[node]):
            if connected and not visited[neighbor]:
                self._dfs(graph, neighbor, visited, component)

    def _connected_components(self, graph):
        n = len(graph)
        visited = [False] * n
        components = []

        for node in range(n):
            if not visited[node]:
                component = []
                self._dfs(graph, node, visited, component)
                components.append(component)
        return components

    def _get_filtered_manual_block(self, d1: dict, d2: dict):
        all_dict = {}
        all_dict.update(d1)
        all_dict.update(d2)
        # 两两比较，获取两两之间所有有交集的M
        all_key = list(all_dict.keys())
        block_matrix = [[0] * len(all_key) for i in range(len(all_key))]
        for i in range(len(all_key)):
            for j in range(i + 1, len(all_key)):
                if set(all_dict[all_key[i]]).intersection(all_dict[all_key[j]]):
                    block_matrix[i][j] = 1
                    block_matrix[j][i] = 1
        # 使用DFS找到所有结合的气泡
        block_list = self._connected_components(block_matrix)
        connect_block = [[all_key[j] for j in i] for i in block_list if len(i) >= 2]
        return connect_block

    def _cat_dict(self, d1: dict, d2: dict):
        d = {}
        d.update(d1)
        for k, v in d2.items():
            if k in d:
                d[k] += v
            else:
                d[k] = v
        return d

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
                            # 如果Virtual节点存在,更新节点物种信息
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
                    if self.have_add_block(all_inner_block):
                        # 气泡中有Manual的情况
                        extend_inner_block = set()
                        for i in all_inner_block:
                            if i.startswith(BaseSyntenyClass.new_block_signal):
                                extend_inner_block.update(BaseSyntenyClass.manual_to_origin[i])
                            else:
                                extend_inner_block.add(i)
                        for i in extend_inner_block:
                            if i not in self.origin_to_manual.keys():
                                self.origin_to_manual[i] = [new_block_name]
                            else:
                                self.origin_to_manual[i].append(new_block_name)
                        for i in all_inner_block:
                            if not i.startswith(BaseSyntenyClass.virtual_block_signal):
                                G.remove_node(i)
                        self.manual_to_origin[new_block_name] = extend_inner_block

                    else:
                        # 普通气泡的情况
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

    def _relabel_and_update_attributes(self, G: nx.Graph, new_to_old: dict) -> nx.Graph:
        # nx.relabel_nodes在多个节点map同个的时候，方法只保留一个节点的属性，同时不累加节点和边的权重，所以需要我们自己手动额外处理
        visited_new_node = set()
        relabel_G = G.copy()
        for new_node, old_manual_list in new_to_old.items():
            # 遍历Manual以及普通block
            all_possible_node = old_manual_list + list(self.manual_to_origin[new_node])
            for old_node in all_possible_node:
                if not G.has_node(old_node):
                    continue

                if new_node not in visited_new_node:
                    relabel_G = nx.relabel_nodes(G, {old_node: new_node}, copy=False)
                    visited_new_node.add(new_node)
                else:
                    neighbors = set(G.neighbors(old_node))
                    if new_node in neighbors:
                        neighbors.remove(new_node)
                    for neighbor_node in neighbors:
                        if G.has_edge(neighbor_node, new_node):
                            G[neighbor_node][new_node]['weight'] += G[neighbor_node][old_node]['weight']
                        else:
                            G.add_edge(neighbor_node, new_node,
                                       weight=G[neighbor_node][old_node]['weight'],
                                       sp_set=deepcopy(G[old_node][neighbor_node]['sp_set']))
                    G.remove_node(old_node)

        return relabel_G
