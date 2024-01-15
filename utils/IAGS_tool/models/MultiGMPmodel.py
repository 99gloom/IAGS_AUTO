from ..inferringAncestorGenomeStructure.EndpointMatchingOptimization import EndpointMatchingOptimization
from ..inferringAncestorGenomeStructure.GMP import GMP
from ..util.transformToAdjacency import transformToAdjacency


def MultiCopyGMPmodel(species_file_list, outdir, guided_species_for_matching,
                      ancestor_name, ancestor_target_copy_number):
    """
    Multi-copy GMP model takes into some species block sequence files and
    transforms block sequence into block adjacencies which is same with GMP.
    But GMP integer programming formulations can just obtain ancestral block adjacencies.
    Ancestral block adjacencies are multi-copy.
    IAGS followed child guide strategy to transform multi-copy ancestral block adjacencies
    to sequences using EMO integer programming formulations.

    :param species_file_list: input species block sequence file list.
    :param outdir: output directory
    :param guided_species_for_matching: a guided child species block sequence file
    :param ancestor_name: ancestor name
    :param ancestor_target_copy_number: target copy number of ancestor species
    """
    adj_file = outdir + ancestor_name + '.adj'
    # transform to adjacencies
    transformToAdjacency(species_file_list, adj_file)
    output_matrix_file = outdir + ancestor_name + '.matrix.xls'
    # GMP integer programming, ancestor target copy number is not one
    gmp = GMP(adj_file,
              target_copy_number=ancestor_target_copy_number)
    gmp.optimization()
    adjacency_matrix = gmp.ancestor_adjacency_matrix()
    adjacency_matrix.output(output_matrix_file)
    # endpoint mathing with guided child species to tansform adjacencies to sequences
    mo = EndpointMatchingOptimization(output_matrix_file, guided_species_for_matching,
                                      matching_dim1=ancestor_target_copy_number,
                                      matching_dim2=ancestor_target_copy_number,
                                      relation1=1,
                                      relation2=1, relabel=True)
    mo.optimization()
    mo.matching_relation()
    adjacency_matrix = mo.build_adjacency_matrix()
    # block sequence
    adjacency_matrix.assemble()
    adjacency_matrix.out_assembly(outdir + ancestor_name + '.block', remove_bar=True)
