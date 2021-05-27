#import random

#import numpy as np
#import pandas as pd
import scipy as sp
from scipy.cluster.hierarchy import average, dendrogram, to_tree
#import skbio
from skbio import Sequence, DNA, TabularMSA, TreeNode, DistanceMatrix
#from skbio.alignment import local_pairwise_align_ssw
from skbio.alignment import global_pairwise_align_nucleotide
from functools import partial

gpa = partial(global_pairwise_align_nucleotide, penalize_terminal_gaps=True)


def kmer_distance(sequence1, sequence2, k=3, overlap=True):
    """Compute the kmer distance between a pair of sequences
    Parameters
    ----------
    sequence1 : skbio.Sequence
    sequence2 : skbio.Sequence
    k : int, optional
        The word length.
    overlapping : bool, optional
        Defines whether the k-words should be overlapping or not
        overlapping.
    Returns
    -------
    float
        Fraction of the set of k-mers from both sequence1 and
        sequence2 that are unique to either sequence1 or
        sequence2.
    Raises
    ------
    ValueError
        If k < 1.
    Notes
    -----
    k-mer counts are not incorporated in this distance metric.
    """
    sequence1_kmers = set(map(str, sequence1.iter_kmers(k, overlap)))
    sequence2_kmers = set(map(str, sequence2.iter_kmers(k, overlap)))
    all_kmers = sequence1_kmers | sequence2_kmers
    shared_kmers = sequence1_kmers & sequence2_kmers
    number_unique = len(all_kmers) - len(shared_kmers)
    fraction_unique = number_unique / len(all_kmers)
    return fraction_unique


def guide_tree_from_sequences(sequences,
                              metric=kmer_distance,
                              display_tree = False):
    """ Build a UPGMA tree by applying metric to sequences
    Parameters
    ----------
    sequences : list of skbio.Sequence objects (or subclasses)
      The sequences to be represented in the resulting guide tree.
    metric : function
      Function that returns a single distance value when given a pair of
      skbio.Sequence objects.
    display_tree : bool, optional
      Print the tree before returning.
    Returns
    -------
    skbio.TreeNode
    """
    guide_dm = DistanceMatrix.from_iterable(
                    sequences, metric=metric, key='id')
    guide_lm = average(guide_dm.condensed_form())
    guide_tree = to_tree(guide_lm)
    if display_tree:
        guide_d = dendrogram(guide_lm, labels=guide_dm.ids, orientation='right',
               link_color_func=lambda x: 'black')
    return guide_tree

def progressive_msa(sequences, pairwise_aligner = gpa, guide_tree=None,
                    metric=kmer_distance):
    """ Perform progressive msa of sequences
    Parameters
    ----------
    sequences : skbio.SequenceCollection
        The sequences to be aligned.
    metric : function, optional
      Function that returns a single distance value when given a pair of
      skbio.Sequence objects. This will be used to build a guide tree if one
      is not provided.
    guide_tree : skbio.TreeNode, optional
        The tree that should be used to guide the alignment process.
    pairwise_aligner : function
        Function that should be used to perform the pairwise alignments,
        for example skbio.alignment.global_pairwise_align_nucleotide. Must
        support skbio.Sequence objects or skbio.TabularMSA objects
        as input.
    Returns
    -------
    skbio.TabularMSA
    """

    if guide_tree is None:
        guide_dm = DistanceMatrix.from_iterable(
                        sequences, metric=metric, key='id')
        guide_lm = average(guide_dm.condensed_form())
        guide_tree = TreeNode.from_linkage_matrix(guide_lm, guide_dm.ids)

    seq_lookup = {s.metadata['id']: s for i, s in enumerate(sequences)}
    c1, c2 = guide_tree.children
    if c1.is_tip():
        c1_aln = seq_lookup[c1.name]
    else:
        c1_aln = progressive_msa(sequences, pairwise_aligner, c1)

    if c2.is_tip():
        c2_aln = seq_lookup[c2.name]
    else:
        c2_aln = progressive_msa(sequences, pairwise_aligner, c2)

    alignment, _, _ = pairwise_aligner(c1_aln, c2_aln)
    # this is a temporary hack as the aligners in skbio 0.4.1 are dropping
    # metadata - this makes sure that the right metadata is associated with
    # the sequence after alignment
    if isinstance(c1_aln, Sequence):
        alignment[0].metadata = c1_aln.metadata
        len_c1_aln = 1
    else:
        for i in range(len(c1_aln)):
            alignment[i].metadata = c1_aln[i].metadata
        len_c1_aln = len(c1_aln)
    if isinstance(c2_aln, Sequence):
        alignment[1].metadata = c2_aln.metadata
    else:
        for i in range(len(c2_aln)):
            alignment[len_c1_aln + i].metadata = c2_aln[i].metadata

    return alignment

def progressive_msa_and_tree(sequences,
                             pairwise_aligner = gpa,
                             metric=kmer_distance,
                             guide_tree=None,
                             display_aln=False,
                             display_tree=False):
    """ Perform progressive msa of sequences and build a UPGMA tree
    Parameters
    ----------
    sequences : skbio.SequenceCollection
        The sequences to be aligned.
    pairwise_aligner : function
        Function that should be used to perform the pairwise alignments,
        for example skbio.alignment.global_pairwise_align_nucleotide. Must
        support skbio.Sequence objects or skbio.TabularMSA objects
        as input.
    metric : function, optional
      Function that returns a single distance value when given a pair of
      skbio.Sequence objects. This will be used to build a guide tree if one
      is not provided.
    guide_tree : skbio.TreeNode, optional
        The tree that should be used to guide the alignment process.
    display_aln : bool, optional
        Print the alignment before returning.
    display_tree : bool, optional
        Print the tree before returning.
    Returns
    -------
    skbio.alignment
    skbio.TreeNode
    """
    msa = progressive_msa(sequences, pairwise_aligner=pairwise_aligner,
                          guide_tree=guide_tree)

    if display_aln:
        print(msa)

    msa_dm = DistanceMatrix.from_iterable(msa, metric=metric, key='id')
    msa_lm = sp.cluster.hierarchy.average(msa_dm.condensed_form())
    msa_tree = TreeNode.from_linkage_matrix(msa_lm, msa_dm.ids)
    if display_tree:
        d = sp.cluster.hierarchy.dendrogram(msa_lm, labels=msa_dm.ids, orientation='right',
                   link_color_func=lambda x: 'black')
    return msa, msa_tree

def iterative_msa_and_tree(sequences,
                           num_iterations,
                           pairwise_aligner = gpa,
                           metric=kmer_distance,
                           display_aln=False,
                           display_tree=False):
    """ Perform progressive msa of sequences and build a UPGMA tree
    Parameters
    ----------
    sequences : skbio.SequenceCollection
       The sequences to be aligned.
    num_iterations : int
       The number of iterations of progressive multiple sequence alignment to
       perform. Must be greater than zero and less than five.
    pairwise_aligner : function
       Function that should be used to perform the pairwise alignments,
       for example skbio.alignment.global_pairwise_align_nucleotide. Must
       support skbio.Sequence objects or skbio.TabularMSA objects
       as input.
    metric : function, optional
      Function that returns a single distance value when given a pair of
      skbio.Sequence objects. This will be used to build a guide tree if one
      is not provided.
    display_aln : bool, optional
       Print the alignment before returning.
    display_tree : bool, optional
       Print the tree before returning.
    Returns
    -------
    skbio.alignment
    skbio.TreeNode
   """
    if num_iterations > 5:
        raise ValueError("A maximum of five iterations is allowed."
                         "You requested %d." % num_iterations)
    previous_iter_tree = None
    for i in range(num_iterations):
        if i == (num_iterations - 1):
            # only display the last iteration
            display = True
        else:
            display = False
        previous_iter_msa, previous_iter_tree = \
            progressive_msa_and_tree(sequences,
             pairwise_aligner=pairwise_aligner,
             metric=metric,
             guide_tree=previous_iter_tree,
             display_aln=display_aln and display,
             display_tree=display_tree and display)

    return previous_iter_msa, previous_iter_tree

