from ..inferringAncestorGenomeStructure.GMP import GMP
from ..util.transformToAdjacency import transformToAdjacency


def GMPmodel(species_file_list, outdir, ancestor_name):
    """
    GMP model takes into some species block sequence files
    and transforms block sequence into block adjacencies.
    IAGS uses GMP integer programming formulations
    based on these block adjacencies to get ancestral block adjacencies
    and then directly transforms to block sequence.

    :param species_file_list: input species block sequence file list.
    :param outdir: output directory
    :param ancestor_name: ancestor name
    """
    adj_file = outdir + ancestor_name + '.adj'
    # transform to adjacencies
    transformToAdjacency(species_file_list, adj_file)
    output_matrix_file = outdir + ancestor_name + '.matrix.xls'
    # GMP integer programming
    gmp = GMP(adj_file, target_copy_number=1)
    gmp.optimization()
    adjacency_matrix = gmp.ancestor_adjacency_matrix()
    adjacency_matrix.output(output_matrix_file)
    # block sequence
    adjacency_matrix.assemble()
    adjacency_matrix.out_assembly(outdir + ancestor_name + '.block', remove_bar=False)
