class Input:

    file_name = ""
    sequences = []
    quality = []

    def __init__(self, file_name):
        self.file_name = file_name
        self.sequences, self.quality = self.readFile()


    def readFastq(self):
        sequences = []
        qualities = []
        with open(self.file_name) as fq:
            while True:
                fq.readline()
                read = fq.readline().strip()
                fq.readline()
                ql = fq.readline().strip()
                if len(read) == 0:
                    break
                sequences.append(read)
                qualities.append(ql)
        return sequences, qualities

    def readFasta(self):
        seq = ''
        sequences = []
        with open(self.file_name) as fa:
            while True:
                read = fa.readline().strip()
                if len(read) == 0:
                    break
                if read[0] == '>':
                    if seq != '':
                        sequences.append(seq)
                    seq = ''
                    continue
                seq += read
            sequences.append(seq)
        return sequences, []


    # Main function to use in reading fasta or fastq files.

    def readFile(self):
        fileID = ''
        with open(self.file_name) as fl:
            fileID = fl.readline()[0]
        if fileID == '@':
            return self.readFastq()
        elif fileID == '>':
            return self.readFasta()
        else:
            return -1, -1



if __name__ == '__main__':
    
    seq = Input('../tests/c.fastq').sequences
    print(seq)
