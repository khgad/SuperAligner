
from skbio import DNA
import matplotlib.pyplot as plt

from msa_algo.algo import *

class MSA:
    seq = []
    msa ,tree = None, None
    
    def __init__(self, sequences, progressive_alignment = True, iterative_alignment = False):
        self.seq = self.creat_DNAseq(sequences)
        if progressive_alignment:
            self.msa, self.tree  = progressive_msa_and_tree(self.seq, display_tree=True)
        else:
            self.msa, self.tree = iterative_msa_and_tree(sequences = self.seq, num_iterations=3, display_tree=True)


    def creat_DNAseq(self, sequences):
        DNA_seq = []
        id = 0
        for s in sequences:
            DNA_seq.append(DNA(s, {"id":"%s%s"%('s', id)}))
            id += 1
        return DNA_seq


if __name__ == '__main__':

    seq = ["ACCGGTGACCAGTTGACCAGT", "ATCGGTACCGGTAGAAGT", "GGTACCAAATAGAA", "GGCACCAAACAGAA", "GGCCCACTGAT"]
    #seq = ['ACTC', 'AAGC', 'AATGC']
    #seq = ['AGGCTATCACCTGACCTCCA', 'TAGCTATCACGACCGC', 'TAGCTGACCGC', 'TCACGACCGACA']
    ds = MSA(seq)
    print(ds.msa)
    #plt.show()
