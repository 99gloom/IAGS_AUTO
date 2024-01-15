import copy
import os
import random

from ..IAGS_tool.models import MultiGMPmodel, GGHPmodel, MultiGGHPmodel, GMPmodel
from ..IAGS_tool.util import calculatedCRBrateAndEstimationAccuracy, chromosomeRearrangementPainting, calculateFissionsAndFusions
from utils.tools.makeDirFunc import make_dir
from ..processTREE.EvolutionaryTree import Evolutionary_tree


class ProcessIAGS:
    def __init__(self, tree_dir: str, final_block_dir: str, output_dir: str, evo_tree: Evolutionary_tree):
        self.computing_queue = self.__read_computing_queue(os.path.join(tree_dir, 'model_and_outgroup.txt'))
        self.nodes_copynumber = self.__read_nodes_copynumber(os.path.join(tree_dir, 'all.ratio'))
        self.used_node = self.__get_used_node()

        self.final_block_dir = final_block_dir
        self.output_dir = output_dir
        self.leave_nodes_id = evo_tree.leaves_nodes_id
        self.evo_tree = evo_tree


        self.multiplied_file_dir = os.path.join(output_dir, 'temp_block_file_for_multiplied')
        make_dir(self.multiplied_file_dir)
        self.painting_dir = os.path.join(output_dir, 'painting')
        make_dir(self.painting_dir)


    def __read_computing_queue(self, path):
        computing_queue = []
        with open(path, 'r') as f:
            for row in f:
                row_info = row.rstrip('\n').rstrip(' ').split(':')
                computing_queue.append(row_info)
        return computing_queue

    def __read_nodes_copynumber(self, path):
        nodes_copynumber = {}
        with open(path, 'r') as f:
            for row in f:
                row_info = row.rstrip('\n').rstrip(' ').split('\t')
                nodes_copynumber[row_info[0]] = int(row_info[1])
        return nodes_copynumber

    def __get_used_node(self):
        used_node = []
        for line in self.computing_queue:
            each_node = []
            anc_name = line[0]
            children_name_list = line[2:-1]
            outgroup_name = line[-1].split('*')[0]
            each_node.append(anc_name)
            for i in children_name_list:
                each_node.append(i)
            each_node.append(outgroup_name)
            used_node.append(each_node)
        return used_node

    def process_ancestor_block_and_evaluate(self):
        for line in self.computing_queue:
            anc_name = line[0]
            mode_type = line[1]
            children_name_list = line[2:-1]
            outgroup_name = line[-1].split('*')[0]
            times = 0
            outgroup_multiplied_flag = False

            # 创造祖先基因组文件夹
            anc_out_dir = os.path.join(self.output_dir, anc_name) + '/'
            make_dir(anc_out_dir)

            # 判断outgroup文件是否需要加倍,MultiCopyGMP需要加倍
            if 'GMP' in mode_type and self.nodes_copynumber[anc_name] / self.nodes_copynumber[outgroup_name] > 1:
                times = int(self.nodes_copynumber[anc_name] / self.nodes_copynumber[outgroup_name])
                outgroup_multiplied_flag = True

            # 获取所有文件的路径
            children_path_list = []
            for i in children_name_list:
                children_path_list.append(self.__get_node_path(i))
            outgroup_file_path = self.__get_node_path(outgroup_name, outgroup_multiplied_flag, times)


            # 根据mode选择模型
            if 'GMP' == mode_type:
                self.__processGMP(children_path_list, children_name_list, outgroup_file_path, outgroup_name,
                                  anc_out_dir, anc_name)
            elif 'MultiCopyGMP' == mode_type:
                self.__processMultiCopyGMP(children_path_list, children_name_list, outgroup_file_path, outgroup_name,
                                           anc_out_dir, anc_name, self.nodes_copynumber[anc_name])
            elif 'GGHP' == mode_type:
                self.__processGGHP(children_path_list, children_name_list, outgroup_file_path, outgroup_name,
                                   anc_out_dir, anc_name, self.nodes_copynumber[children_name_list[0]],
                                   self.nodes_copynumber[outgroup_name])
            elif 'MultiCopyGGHP' == mode_type:
                self.__processMultiCopyGGHP(children_path_list, children_name_list, outgroup_file_path, outgroup_name,
                                            anc_out_dir, anc_name, self.nodes_copynumber[children_name_list[0]],
                                            self.nodes_copynumber[outgroup_name], self.nodes_copynumber[anc_name])
            else:
                print('File error! please check file or rerun')
                exit()



    def process_painting(self):
        color_list = ['#DF1159','#1E93C9','#26AF67','#D5A1C5','#EBCA6D','#94B51E','#000000','#A9A9A9','#62C1BD',
                      '#FF8C00','#4169E1','#FF00FF','#9932CC','00FFFF','00FF00']
        block_length_file = os.path.join(self.final_block_dir,'blockindex.genenumber')
        target_species_block_file = self.__get_node_path(self.evo_tree.painting_anc_node)
        target_species_name = self.evo_tree.painting_anc_node
        target_species_copy_number = self.nodes_copynumber[self.evo_tree.painting_anc_node]

        # 处理颜色列表不够的情况
        HEX_LIST = ['0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F']
        chr = []
        with open(target_species_block_file,'r') as f:
            for row in f:
               chr.append(row.rstrip('\n').rstrip(' '))
        while (len(color_list) < len(chr)):
            color = '#'
            for j in range(6):
                color += HEX_LIST[random.randint(0,15)]
            if color not in color_list:
                color_list.append(color)

        painting_list = copy.deepcopy(self.evo_tree.all_nodes_id_without_root)
        painting_list.remove(self.evo_tree.painting_anc_node)
        for i in painting_list:
            rearranged_species_block_file = self.__get_node_path(i)
            rearranged_species_name = i
            rearranged_species_copy_number = self.nodes_copynumber[i]
            chromosomeRearrangementPainting.plotChrsRearrangement(block_length_file,
                                  rearranged_species_block_file, rearranged_species_name,
                                  rearranged_species_copy_number,target_species_block_file,
                                  target_species_name, target_species_copy_number,
                                  color_list, self.painting_dir + '/')

    def process_Calculating_Fissions_and_Fusions(self):
        Fissions_and_Fusions_path = os.path.join(self.output_dir, 'shufflingEvents.txt')
        outfile =  open(Fissions_and_Fusions_path,'w')
        tree = self.evo_tree.tree
        node_id_list = []
        node_id_list.append('root')
        while node_id_list:
            if node_id_list[0] == 'root':
                node_id_list.pop(0)
                if not tree.children('root'):
                    continue

                root_children_list = []
                for i in tree.children('root'):
                    root_children_list.append(i.identifier)
                    node_id_list.append(i.identifier)
                root_children_list.remove(self.evo_tree.painting_anc_node)
                ancestor_file = self.__get_node_path(self.evo_tree.painting_anc_node)
                ancestor_copy_number = self.nodes_copynumber[self.evo_tree.painting_anc_node]
                ancestor_block_sequence_dir = os.path.join(self.output_dir, self.evo_tree.painting_anc_node) + '/'

                for j in root_children_list:
                    descendant_file = self.__get_node_path(j)
                    descendant_copy_number = self.nodes_copynumber[j]
                    fissions, fusions = calculateFissionsAndFusions.calculateFissionAndFussions(descendant_file, ancestor_file,
                                                                    descendant_copy_number, ancestor_copy_number,
                                                                    ancestor_block_sequence_dir)
                    outfile.write(self.evo_tree.painting_anc_node + '->' +  j + '\n' +
                                  'fissions: ' + str(fissions) + '\t' + 'fusions: ' + str(fusions) + '\n')
            else:
                ancestor_name = node_id_list.pop(0)
                if not tree.children('root'):
                    continue

                root_children_list = []
                for i in tree.children(ancestor_name):
                    root_children_list.append(i.identifier)
                    node_id_list.append(i.identifier)
                ancestor_file = self.__get_node_path(ancestor_name)
                ancestor_copy_number = self.nodes_copynumber[ancestor_name]
                ancestor_block_sequence_dir = os.path.join(self.output_dir, ancestor_name) + '/'
                for j in root_children_list:
                    descendant_file = self.__get_node_path(j)
                    descendant_copy_number = self.nodes_copynumber[j]
                    fissions, fusions = calculateFissionsAndFusions.calculateFissionAndFussions(descendant_file, ancestor_file,
                                                                    descendant_copy_number, ancestor_copy_number,
                                                                    ancestor_block_sequence_dir)
                    outfile.write(ancestor_name + '->' +  j + '\n' +
                                  'fissions: ' + str(fissions) + '\t' + 'fusions: ' + str(fusions) + '\n')

        outfile.close()


    def __processGMP(self, children_path_list, children_name_list, outgroup_path, outgroup_name, anc_out_dir, anc_name):
        # 计算节点
        species_file_list = []
        for i in children_path_list:
            species_file_list.append(i)
        species_file_list.append(outgroup_path)
        GMPmodel.GMPmodel(species_file_list, anc_out_dir, anc_name)

        # 评估节点
        all_calculate_nodes_name = [outgroup_name] + children_name_list
        all_calculate_nodes_path = [outgroup_path] + children_path_list
        matching_target_file_name = None
        matching_target_file_path = None
        for i in range(len(all_calculate_nodes_name)):
            if all_calculate_nodes_name[i] in self.leave_nodes_id:
                matching_target_file_name = all_calculate_nodes_name[i]
                matching_target_file_path = all_calculate_nodes_path[i]
                break
        if not matching_target_file_name:
            matching_target_file_name = outgroup_name
        if not matching_target_file_path:
            matching_target_file_path = outgroup_path
        matching_target_copy_number = 1
        speciesAndCopyList = []
        for i in range(len(all_calculate_nodes_name)):
            species_temp_info = []
            species_temp_info.append(all_calculate_nodes_path[i])
            species_temp_info.append(1)
            species_temp_info.append(all_calculate_nodes_name[i])
            speciesAndCopyList.append(species_temp_info)
        calculatedCRBrateAndEstimationAccuracy.calculatedCRBrateAndEstimationAccuracy(matching_target_file_path,
                                                                                      matching_target_copy_number,
                                                                                      matching_target_file_name,
                                                                                      speciesAndCopyList,
                                                                                      anc_out_dir,
                                                                                      'GMP')
        return

    def __processGGHP(self, children_path_list, children_name_list, outgroup_path, outgroup_name, anc_out_dir, anc_name,
                      dup_copy_number, out_copy_number):

        GGHPmodel.GGHPmodel(children_path_list[0], outgroup_path, anc_out_dir, anc_name, dup_copy_number,
                            out_copy_number)

        child_elem = [children_path_list[0], dup_copy_number, children_name_list[0]]
        out_elem = [outgroup_path, out_copy_number, outgroup_name]
        speciesAndCopyList = [child_elem, out_elem]

        calculatedCRBrateAndEstimationAccuracy.calculatedCRBrateAndEstimationAccuracy(outgroup_path,
                                                                                      out_copy_number,
                                                                                      outgroup_name,
                                                                                      speciesAndCopyList,
                                                                                      anc_out_dir,
                                                                                      'GGHP')

    def __processMultiCopyGGHP(self, children_path_list, children_name_list, outgroup_path, outgroup_name, anc_out_dir,
                               anc_name, dup_copy_number, out_copy_number, ancestor_target_copy_number):

        MultiGGHPmodel.MultiCopyGGHPmodel(children_path_list[0], outgroup_path, anc_out_dir, anc_name, dup_copy_number,
                                          out_copy_number, ancestor_target_copy_number)

        speciesAndCopyList = []
        child_elem = [children_path_list[0], dup_copy_number, children_name_list[0]]
        out_elem = [outgroup_path, out_copy_number, outgroup_name]
        speciesAndCopyList = [child_elem, out_elem]

        calculatedCRBrateAndEstimationAccuracy.calculatedCRBrateAndEstimationAccuracy(outgroup_path,
                                                                                      out_copy_number,
                                                                                      outgroup_name,
                                                                                      speciesAndCopyList,
                                                                                      anc_out_dir,
                                                                                      'GGHP')

    def __processMultiCopyGMP(self, children_path_list, children_name_list, outgroup_path, outgroup_name, anc_out_dir,
                              anc_name, ancestor_target_copy_number=2):
        # 计算
        species_file_list = []
        for i in children_path_list:
            species_file_list.append(i)
        species_file_list.append(outgroup_path)
        MultiGMPmodel.MultiCopyGMPmodel(species_file_list, anc_out_dir, outgroup_path, anc_name,
                                        ancestor_target_copy_number)

        # 评估
        all_calculate_nodes_name = [outgroup_name] + children_name_list
        all_calculate_nodes_path = [outgroup_path] + children_path_list
        matching_target_file_name = None
        matching_target_file_path = None
        for i in range(len(all_calculate_nodes_name)):
            if all_calculate_nodes_name[i] in self.leave_nodes_id:
                matching_target_file_name = all_calculate_nodes_name[i]
                matching_target_file_path = all_calculate_nodes_path[i]
                break
        if not matching_target_file_name:
            matching_target_file_name = outgroup_name
        if not matching_target_file_path:
            matching_target_file_path = outgroup_path
        matching_target_copy_number = ancestor_target_copy_number
        speciesAndCopyList = []
        for i in range(len(all_calculate_nodes_name)):
            species_temp_info = []
            species_temp_info.append(all_calculate_nodes_path[i])
            species_temp_info.append(ancestor_target_copy_number)
            species_temp_info.append(all_calculate_nodes_name[i])
            speciesAndCopyList.append(species_temp_info)
        calculatedCRBrateAndEstimationAccuracy.calculatedCRBrateAndEstimationAccuracy(matching_target_file_path,
                                                                                      matching_target_copy_number,
                                                                                      matching_target_file_name,
                                                                                      speciesAndCopyList,
                                                                                      anc_out_dir,
                                                                                      'MultiCopyGMP')
        return

    def __multiply_file(self, infile, outfile, times):

        outfile = open(outfile, 'w')
        sequence = []
        with open(infile, 'r') as f:
            while True:
                line = f.readline()
                if not line:
                    break
                sequence.append(line)
        for j in range(times):
            for i in sequence:
                outfile.write(i)

    def __get_node_path(self, file_name, outgroup_multiplied_flag=False, times=1) -> str:
        if outgroup_multiplied_flag:
            multiplied_outgroup_path = os.path.join(self.multiplied_file_dir,
                                                    file_name + '_' + str(times) + '.multiplied.block')
            self.__multiply_file(self.__get_node_path(file_name), multiplied_outgroup_path, times)
            return multiplied_outgroup_path

        else:
            if file_name.startswith('ANC') or file_name.startswith('PRE'):
                return os.path.join(self.output_dir, file_name + '/' + file_name + '.block')
            else:
                return os.path.join(self.final_block_dir, file_name + '.final.block')
