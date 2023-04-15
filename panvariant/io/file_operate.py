class FastaIO:
    def __init__(self, fasta_file):
        self.__fasta_file = fasta_file
        self.fasta_db = {}

    def read_fasta(self):
        with open(self.__fasta_file, 'r') as fin:
            gid = ""
            seq = ""
            for line in fin:
                if line[0] == '>':
                    if seq != "":
                        self.fasta_db[gid] = seq
                    gid = line.strip().split()[0][1:]
                    seq = ""
                else:
                    seq += line.strip().upper()
        self.fasta_db[id] = seq

    def read_aln(self):
        with open(self.__fasta_file, 'r') as fin:
            gid = ""
            seq = ""
            for line in fin:
                if line[0] == '>':
                    if seq != "":
                        self.fasta_db[gid] = seq
                    gid = line.strip().split()[0][1:].replace('_R_', '')
                    seq = ""
                else:
                    seq += line.strip().upper()
            self.fasta_db[gid] = seq
