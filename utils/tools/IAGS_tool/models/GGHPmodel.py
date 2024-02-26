from ..inferringAncestorGenomeStructure.GGHP import GGHP
from ..util.transformToAdjacency import transformToAdjacency


def GGHPmodel(dup_child_file, outgroup_file, outdir,
              ancestor_name, dup_copy_number, out_copy_number):
    """
    GGHP model takes into duplicated and outgroup species block sequences.
    Ancestor block copy number should be only one.
    IAGS transforms both block sequences into block adjacencies.
    IAGS uses GGHP integer programming formulations
    based on block adjacencies to get ancestral block adjacencies
    and then directly transforms to block sequence.
    For basic GGHP, target copy number of duplicated species is two and outgroup species is one.
    IAGS allow multiple species as input
    which duplicated species block sequences and outgroup species block sequences
    should be merged together, respectively and
    the input target block copy number should be summed, respectively.

    :param dup_child_file: block sequence file for duplicated species
    :param outgroup_file: block sequence file for outgroup species
    :param outdir: output directory
    :param ancestor_name: ancestor name
    :param dup_copy_number: target copy number of duplicated species
    :param out_copy_number: target copy number of outgroup species
    """
    filelist = [dup_child_file, outgroup_file]
    adj_file = outdir + ancestor_name + '.adj'
    output_matrix_file = outdir + ancestor_name + '.matrix.xls'
    # transform to adjacencies
    transformToAdjacency(filelist, adj_file)
    # GGHP integer programming, ancestor target copy number is only one for GGHP model
    ggap = GGHP(adj_file, target_copy_number=1, duptype=dup_copy_number / 1, duptype_out=out_copy_number / 1)
    ggap.optimization()
    adjacency_matrix = ggap.ancestor_adjacency_matrix()
    adjacency_matrix.output(output_matrix_file)
    # block sequence
    adjacency_matrix.assemble()
    adjacency_matrix.out_assembly(outdir + ancestor_name + '.block', remove_bar=False)
