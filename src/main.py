from arg import commandline
from needle_man import global_pairwise_align
from smith_waterman import local_pairwise_align
from MSA import MSA
from input import Input
from trans import translator
from charts import MGM_chart, plot_fastq_qualities

if __name__ == '__main__':
    
    args = commandline()
    
    seq, qul = [], []
    file_path = ''
    
    
    if args.file:
        file_path = args.file
        _in = Input(file_path)
        seq, qul = _in.sequences, _in.quality
        
    elif args.sequence:
        seq = args.sequence
    
    if args.quality:
        qul = args.quality
       
    
    if args.local:
        #
        local_align = local_pairwise_align(seq[0], seq[1])
        print(local_align.align1)
        print(local_align.align2)
        
    elif args.Global:
        # seq[0] = AGAT --- seq[1] = ATGT
        global_align = global_pairwise_align(seq[0], seq[1])
        print(global_align.align1)
        print(global_align.align2)
        
    elif args.progressive:
        pro_align = MSA(seq)
        print(pro_align.msa)
        
    elif args.iterative:
        ite_align = MSA(seq, False, True)
        print(ite_align.msa)
    
    if args.quality_chart and len(qul):
        plot_fastq_qualities(qul)
        
    if args.score_chart:
        MGM_chart(17,4,10)
            
    if args.translator:
        print(translator(seq))