from ..inferringAncestorGenomeStructure.BlockMatchingOptimization import BlockMatchingOptimization
from ..inferringAncestorGenomeStructure.EndpointMatchingOptimization import EndpointMatchingOptimization
from ..inferringAncestorGenomeStructure.GGHP import GGHP
from ..util.transformToAdjacency import transformToAdjacency


def MultiCopyGGHPmodel(dup_child_file, outgroup_file, outdir,
                       ancestor_name, dup_copy_number, out_copy_number, ancestor_target_copy_number):
    """
    Multi-copy GGHP model takes into duplicated and outgroup species block sequences.
    Ancestor block copy number can be more than one.
    IAGS transforms both block sequences into block adjacencies which is same with GGHP.
    But GGHP integer programming formulations can just obtain ancestral block adjacencies.
    Ancestral block adjacencies are multi-copy.
    IAGS followed child guide strategy to transform multi-copy ancestral block adjacencies to sequences.
    IAGS first used self-BMO integer programming formulation
    to remove the influence of WGD in child species and then used EMO integer programming formulation.

    :param dup_child_file: block sequence file for duplicated species
    :param outgroup_file: block sequence file for outgroup species
    :param outdir: output directory
    :param ancestor_name: ancestor name
    :param dup_copy_number: target copy number of duplicated species
    :param out_copy_number: target copy number of outgroup species
    :param ancestor_target_copy_number: target copy number of ancestor species
    """
    adj_file = outdir + ancestor_name + '.adj'
    filelist = [dup_child_file,
                outgroup_file]
    # transform to adjacencies
    transformToAdjacency(filelist, adj_file)
    # GGHP integer programming, ancestor target copy number is not one
    ggap = GGHP(adj_file, target_copy_number=ancestor_target_copy_number,
                duptype=dup_copy_number / ancestor_target_copy_number,
                duptype_out=out_copy_number / ancestor_target_copy_number)
    ggap.optimization()
    adjacency_matrix = ggap.ancestor_adjacency_matrix()
    adjacency_matrix.output(outdir + ancestor_name + '.matrix.xls')
    # self-BMO for reduce matching dimension
    mo = BlockMatchingOptimization(dup_child_file, dup_child_file,
                                   matching_dim1=dup_copy_number,
                                   matching_dim2=dup_copy_number,
                                   relation1=1,
                                   relation2=1,
                                   self_matching=True)
    mo.optimization()
    mo.matching_relation()
    outmatching = outdir + 'dup_self_matching.table.txt'
    mo.output_matching_relation(outmatching)
    output_sequence_file_list = [outdir + 'dup_self_matching.block.relabel']
    mo.out_relabel_sequence(output_sequence_file_list)

    # endpoint mathing with self matched guided child species to tansform adjacencies to sequences

    mo = EndpointMatchingOptimization(outdir + ancestor_name + '.matrix.xls',
                                      outdir + 'dup_self_matching.block.relabel',
                                      matching_dim1=ancestor_target_copy_number,
                                      matching_dim2=ancestor_target_copy_number,
                                      relation1=1,
                                      relation2=1,
                                      relabel=False)

    mo.optimization()
    mo.matching_relation()
    adjacency_matrix = mo.build_adjacency_matrix()

    # block sequence
    adjacency_matrix.assemble()
    adjacency_matrix.out_assembly(outdir + ancestor_name + '.block', remove_bar=True)
