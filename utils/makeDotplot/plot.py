import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
import pandas as pd



def get_gff(gff_path):
    # 获取gff中每个genename和[start，end]的对应
    gene_range = {}
    with open(gff_path, 'r') as gffFile:
        for row in gffFile:
            gene_info = row.rstrip('\n').rstrip(' ').split('\t')
            gene_range[gene_info[1]] = [int(gene_info[2]), int(gene_info[3])]
    return gene_range

def get_chr_range(gff_path):
        dataFrame = pd.read_csv(gff_path, header=None, sep='\t')
        name = dataFrame.at[0, 0]
        number_index = 0
        for j in range(len(name) - 1, -1, -1):
            if name[j].isdigit():
                continue
            else:
                number_index = j + 1
                break
        dataFrame[0] = dataFrame[0].map(lambda x: int(x[number_index:]))
        dataFrame = dataFrame.sort_values(by=[0, 2], ascending=[True, True])
        dataFrame[0] = dataFrame[0].map(lambda x: name[:number_index] + str(x))
        dataFrame[1] = dataFrame[1].map(lambda x: x.split(';')[0])

        chr_name = dataFrame[0].unique()
        chr_len = []
        for i in chr_name:
            row = dataFrame[dataFrame[0] == i][-1:]
            length = row[3]
            chr_len.append(int(length))
        return chr_len



def get_name_and_number(sequence_path, sequence_name_path, gene_range):
    sequence_id = []
    sequence_name = []
    with open(sequence_path, 'r') as sequenceFile:
        # 读取sequence文件
        for row in sequenceFile:
            row_info = row.rstrip('\n').rstrip(' ').split(' ')
            sequence_id.append(row_info)

    with open(sequence_name_path, 'r') as nameFile:
        # 读取sequence name文件
        for row in nameFile:
            row_info = row.rstrip('\n').rstrip(' ').split(' ')
            sequence_name.append(row_info)

    for i in range(len(sequence_id)):
        for j in range(len(sequence_id[i])):
            id = sequence_id[i][j]
            name = sequence_name[i][j]
            if name in gene_range.keys():
                gene_range[name].append(int(id))
    return sequence_name, gene_range


def abs_plot(speciesA_name, speciesB_name, sequence_dir, gff_dir, save_dir):
    speciesA_sequence_path = os.path.join(sequence_dir, speciesA_name) + '.all.sequence'
    speciesA_sequence_name_path = os.path.join(sequence_dir, speciesA_name) + '.all.sequence.genename'
    speciesA_gff_path = os.path.join(gff_dir, speciesA_name) + '.gff'

    speciesB_sequence_path = os.path.join(sequence_dir, speciesB_name) + '.all.sequence'
    speciesB_sequence_name_path = os.path.join(sequence_dir, speciesB_name) + '.all.sequence.genename'
    speciesB_gff_path = os.path.join(gff_dir, speciesB_name) + '.gff'

    save_dir = os.path.join(save_dir, speciesA_name + '-' + speciesB_name + '.all.pdf')
    # print(speciesA_name + '-' + speciesB_name + '-plot')

    gene_range_A = get_gff(speciesA_gff_path)
    gene_range_B = get_gff(speciesB_gff_path)
    chr_len_A = get_chr_range(speciesA_gff_path)
    chr_len_B = get_chr_range(speciesB_gff_path)

    sequence_name_A, gene_range_A = get_name_and_number(speciesA_sequence_path, speciesA_sequence_name_path,
                                                        gene_range_A)

    sequence_name_B, gene_range_B = get_name_and_number(speciesB_sequence_path, speciesB_sequence_name_path,
                                                        gene_range_B)

    sum_chr_len_A = sum(chr_len_A)
    sum_chr_len_B = sum(chr_len_B)

    fig = plt.figure(figsize=(40, 40))
    ax = fig.add_subplot(1, 1, 1)
    ax.axis('equal')
    ax.axis('off')
    start_plot_x = 0
    start_plot_y = 0

    for n in range(len(sequence_name_A)):
        for i in range(len(sequence_name_B)):
            for k in sequence_name_A[n]:
                for j in sequence_name_B[i]:
                    if gene_range_A[k][-1] == gene_range_B[j][-1]:
                        ax.add_patch(
                            patches.Rectangle(
                                (start_plot_x + gene_range_B[j][0], start_plot_y + gene_range_A[k][0]),
                                500000, 500000,
                                # 500K可视性更高
                                # gene_range_B[j][1] - gene_range_B[j][0],
                                # gene_range_A[k][1] - gene_range_A[k][0],
                                facecolor='red'
                            )

                        )
            start_plot_x += chr_len_B[i]
        start_plot_y += chr_len_A[n]
        start_plot_x = 0
        print('chr_y:' + str(start_plot_y))

    ax.hlines(0, 0, sum_chr_len_B, color="black")
    ax.vlines(0, 0, sum_chr_len_A, color="black")

    line_position_x = 0
    for i in range(len(chr_len_B)):
        line_position_x += chr_len_B[i]
        ax.vlines(line_position_x, 0, sum_chr_len_A, color='black')

    line_position_y = 0
    for i in range(len(chr_len_A)):
        line_position_y += chr_len_A[i]
        ax.hlines(line_position_y, 0, sum_chr_len_B, color='black')



    plt.text(x=-sum_chr_len_B / 6, y=sum_chr_len_A / 2, fontsize=50, verticalalignment='center',
             s=speciesA_name)
    plt.text(x=sum_chr_len_B / 2, y=-sum_chr_len_A / 9, fontsize=50, verticalalignment='center',
             s=speciesB_name)

    fig.savefig(save_dir)

def process_plot(species, sequence_dir, gff_dir, save_dir):
    for i in range(len(species)):
        for j in range(i+1, len(species)):
            abs_plot(species[i], species[j], sequence_dir, gff_dir, save_dir)



