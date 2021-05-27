import argparse


def commandline():

    parser = argparse.ArgumentParser(description='python script for pairwise and multiple sequence alignment')
    parser.add_argument('-Q', '--quality_chart',
                        action='store_true',
                        help='display sequences quality score chart'
                        )
    parser.add_argument('-G', '--score_chart',
                        action='store_true', 
                        help='display chart of rate of gaps, mismatches, and matches'
                        )
    parser.add_argument('-T', '--translator',
                        action='store_true',
                        help='translate DNA sequence into protien'
                        )

    alignment_type = parser.add_mutually_exclusive_group(required=True)
    alignment_type.add_argument('-S', '--local', 
                                action='store_true', 
                                help='smith-waterman algorithm'
                                )
    alignment_type.add_argument('-N', '--Global',
                                action='store_true',
                                help='needleman wunsch algorithm'
                                )
    alignment_type.add_argument('-P', '--progressive', 
                                action='store_true', 
                                help='progressive algorithm'
                                )
    alignment_type.add_argument('-I', '--iterative', 
                                action='store_true', 
                                help='iterative algorithm'
                                )
    
    input_type = parser.add_mutually_exclusive_group(required = True)
    input_type.add_argument('-f', '--file',
                            metavar='path',
                            help='path of vaild file [ fasta | fastq ]')
    input_type.add_argument('-s', '--sequence', 
                            nargs='+', 
                            metavar='seq', 
                            help='sequences to be aligned', 
                            )
    
    parser.add_argument('-q', '--quality', 
                        nargs='+', 
                        metavar='qul', 
                        help='quality of the sequences'
                        )
    
    '''
    alignment_type = parser.add_subparsers(help='type of the alignment method')

    pairwise = alignment_type.add_parser('pw', help='pairwise sequence alignment')
    m = pairwise.add_mutually_exclusive_group(required = True)
    m.add_argument('-S', '--local', 
                action='store_true', 
                help='smith-waterman algorithm'
                )
    m.add_argument('-N', '--Global',
                action='store_true',
                help='needleman wunsch algorithm'
                )
    input_type = pairwise.add_subparsers(help='type of the input [ file | enter text ]')
    file = input_type.add_parser('f', help = 'provide vaild files [ fasta | fastq ]')
    file.add_argument('file', help='path of the file', default=None)
    seq = input_type.add_parser('s', help='give sequence by hand')
    seq.add_argument('-s', '--sequence', 
                    nargs=2, 
                    metavar='seq', 
                    help='sequences to be aligned', 
                    required=True
                    )
    seq.add_argument('-q', '--quality', 
                    nargs=2, 
                    metavar='qul', 
                    help='quality of the sequences'
                    )

    multiple = alignment_type.add_parser('mp', help='multiple sequence alignment')
    m = multiple.add_mutually_exclusive_group(required = True)
    m.add_argument('-P', '--progressive', 
                action='store_true', 
                help='progressive algorithm'
                )
    m.add_argument('-I', '--iterative', 
                action='store_true', 
                help='iterative algorithm'
                )
    input_type = multiple.add_subparsers(help='type of the input [ file | enter text ]')
    file = input_type.add_parser('f', help = 'provide vaild files [ fasta | fastq ]')
    file.add_argument('file', help='path of the file')
    seq = input_type.add_parser('s', help='give sequence by hand')
    seq.add_argument('-s', '--sequence', 
                    nargs='+', 
                    metavar='seq', 
                    help='sequences to be aligned', 
                    required=True
                    )
    seq.add_argument('-q', '--quality', 
                    nargs='+', 
                    metavar='qul', 
                    help='quality of the sequences'
                    )
    '''
    return parser.parse_args()

if __name__ == '__main__':
    
    commandline()
    