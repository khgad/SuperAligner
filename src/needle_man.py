import numpy as np

class global_pairwise_align:
    
    match=1
    mismatch=-1
    gap=-1
    seq1=""
    seq2=""
    mtx=None
    align1=None
    align2=None

    def __init__(self,seq1,seq2):
        self.seq1=seq1
        self.seq2=seq2
        self.mtx=self.matrix_score(seq1, seq2)
        self.align1,self.align2=self.trace_back(self.mtx, seq1, seq2)

    def match_score(self,c1,c2):
        if c1==c2:
            return self.match
        elif c1=="-" or c2=="-":
            return self.gap
        else:
            return self.mismatch

    def matrix_score(self, seq1, seq2):
        x=0
        m,n=len(seq1),len(seq2)
        matrix_score=np.zeros((m+1,n+1),dtype=int)
        for i in range(0,m+1):
            matrix_score[i][0]=self.gap*i
        for j in range(0,n+1):
            matrix_score[0][j]=self.gap*j

        for i in range(1,m+1):
            for j in range(1,n+1):
                diagonal=matrix_score[i-1][j-1]+self.match_score(seq1[i-1],seq2[j-1])
                up=matrix_score[i-1][j]+self.gap
                left=matrix_score[i][j-1]+self.gap
                matrix_score[i][j]=max(diagonal,up,left)
        return matrix_score 
                
    def trace_back(self,matrix_score, seq1, seq2):
        align1,align2='',''
        i,j=len(seq1),len(seq2)

        while i>0 or j>0:
            score_current=matrix_score[i][j]
            score_diagonal=matrix_score[i-1][j-1]
            score_left=matrix_score[i][j-1]
            score_up=matrix_score[i-1][j]

            if score_current==score_diagonal+self.match_score(seq1[i-1],seq2[j-1]):
                a1,a2=seq1[i-1],seq2[j-1]
                i,j=i-1,j-1
            elif score_current==score_up+self.gap:
                a1,a2=seq1[i-1],"-"
                i-=1
            elif score_current==score_left+self.gap:
                a1,a2="-",seq2[j-1]
                j-=1

            align1+=a1
            align2+=a2
        align1=align1[::-1]
        align2=align2[::-1]
        return align1,align2

    
    