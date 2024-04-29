from copy import deepcopy
from typing import List
import networkx as nx

class BaseSyntenyClass:
    new_block_signal = 'Manual'
    new_block_number = 0
    manual_to_origin = {} # set
    origin_to_manual = {} # list

    virtual_block_signal = 'Virtual'
    virtual_block_number = 0
    virtual_block_adj_info = {}
    adj_tuple_to_virtual_block = {}
    def __init__(self):
        self.manual_to_origin = {}
        self.origin_to_manual = {}
        self.G = nx.Graph()
        self.cpn = {}

    @classmethod
    def is_equal_edge_weight(cls, l: List[str], G: nx.Graph):
        inner_flag = True
        for i in range(1, len(l) - 1, 1):
            if G[l[i]][l[i + 1]]['weight'] != G[l[1]][l[2]]['weight']:
                inner_flag = False
        return inner_flag


    @classmethod
    def get_manual_block(cls):
        res =  BaseSyntenyClass.new_block_signal + '_' + str(BaseSyntenyClass.new_block_number)
        BaseSyntenyClass.new_block_number += 1
        return res

    @classmethod
    def get_virtual_block(cls):
        res = BaseSyntenyClass.virtual_block_signal  + '_' + str(BaseSyntenyClass.virtual_block_number)
        BaseSyntenyClass.virtual_block_number += 1
        return res

    def have_add_block(cls, l: set):
        for i in l:
            if i.startswith(BaseSyntenyClass.new_block_signal):
                return i
        return None

    def add_new_node(self, G, new_block_name, start, end, edge_weight, last_block, sp_set):
        G.add_node(new_block_name, times=G.nodes[last_block]['times'])
        G.add_edge(new_block_name, start, weight=edge_weight, sp_set=deepcopy(sp_set))
        G.add_edge(new_block_name, end, weight=edge_weight, sp_set=deepcopy(sp_set))

    def out_res(self, dir , root_name, virtual_name):

        def delete_illegal_content(m2b,b2m):
            res = set()
            for k,v in b2m.items():
                count = 0
                for i in v:
                    i:str
                    if i.startswith('Manual'):
                        count += 1
                if count > 1:
                    res.add(k)
                    res.update(v)
            res_list = list(res)
            index = 0
            while index < len(res_list):
                if res_list[index].startswith('Manual'):
                    res_list.extend(m2b[res_list[index]])
                index += 1

            manual_to_origin = {key: value for key, value in self.manual_to_origin.items() if key not in res_list}
            origin_to_manual = {key: value for key, value in self.origin_to_manual.items() if key not in res_list}
            return manual_to_origin, origin_to_manual

        self.manual_to_origin,self.origin_to_manual = delete_illegal_content(self.manual_to_origin,self.origin_to_manual)

        virtual_dict = {}
        with open(f'./{dir}/{root_name}-M2B.txt', 'w') as root_f:
            for k, v in self.manual_to_origin.items():
                root_f.write(k + ':')
                for i in v:
                    root_f.write(i + ' ')
                    if i.startswith(BaseSyntenyClass.virtual_block_signal) and i not in virtual_dict:
                        virtual_dict[i] = BaseSyntenyClass.virtual_block_adj_info[i]
                root_f.write('\n')
        with open(f'./{dir}/{root_name}-B2M.txt', 'w') as root_f:
            for k, v in self.origin_to_manual.items():
                root_f.write(k + ':')
                for i in v:
                    root_f.write(i + ' ')
                root_f.write('\n')

        with open(f'./{dir}/{virtual_name}.txt', 'w') as virtual_f:
            for k, v in virtual_dict.items():
                # num
                write_list = []
                write_list.append(k)
                # adj
                adj = v[0]
                adj_str = ''
                for j in adj:
                    adj_str += j + ' '
                write_list.append(adj_str[:-1])
                # sp
                sp = v[1]
                sp_str = ''
                for j in sp:
                    sp_str += j + ' '
                write_list.append(sp_str)
                virtual_f.write(':'.join(write_list) + '\n')

