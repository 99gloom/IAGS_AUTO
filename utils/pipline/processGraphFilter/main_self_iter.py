from utils.pipline.processGraphFilter.graphAlgorithm.BaseSyntenyClass import BaseSyntenyClass
from utils.pipline.processGraphFilter.graphAlgorithm.SelfSyntenyGraph import SelfSyntenyGraph
from utils.pipline.processGraphFilter.graphAlgorithm.MergedSyntenyGraph import MergedSyntenyGraph

if __name__ == '__main__':
    '''
    输入格式  数据需要额外生成
    '''
    sp_list = [
        'Brachy',
        'Sorghum',
        'Maize',
        'Rice',
        'Telongatum'
    ]
    sp_copy_number = {
        'Brachy': 2,
        'Sorghum': 2,
        'Maize': 4,
        'Rice': 2,
        'Telongatum': 2,
    }
    order = {
        'ANC3': ['Maize', 'Sorghum'],
        'ANC2': ['Brachy', 'Telongatum'],
        'ANC1': ['ANC2', 'Rice'],
        'root': ['ANC1', 'ANC3']
    }
    dir = './example/drimmBlocks/'


    '''
    匹配算法
    '''
    # 自匹配
    syn = {}
    for i in sp_list:
        g = SelfSyntenyGraph(i, dir, sp_copy_number[i])
        syn[i] = g


    # 自底向上匹配
    for k, v in order.items():
        g = MergedSyntenyGraph(syn[v[0]], syn[v[1]])
        syn[k] = g


    '''
    输出文件
    '''
    # 输出
    root_name = 'root'
    virtual_name = 'virtual'
    res_dir = './graph_process_res/'
    syn['root'].out_res(res_dir, root_name, virtual_name)