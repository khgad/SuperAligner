import numpy as np

class pairwise:
    
    match=1
    mismatch=-1
    gap=-1
    seq1=""
    seq2=""
    mtx=None
    align1=None
    align2=None
    
    def match_score(self, c1,c2):
        if c1==c2:
            return self.match
        elif c1=='-' or c2 == '-':
            return self.gap
        else:
            return self.mismatch