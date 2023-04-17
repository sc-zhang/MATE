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
        self.fasta_db[gid] = seq

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


class VariantIO:
    def __init__(self):
        self.samples = None
        # self.variants is 2-D list, like:
        # [[pos, ref, alt(dict), [0, 1, 1, 0, ...](sample genotype)
        # ...
        # ]
        self.variants = []

    def read_var(self, variant_file):
        with open(variant_file, 'r') as fin:
            for line in fin:
                data = line.strip().split()
                if line[0] == '#':
                    self.samples = data[3:]
                else:
                    pos = int(data[0])
                    ref = data[1]
                    alt_sites = data[2].split(',')
                    alt = {_+1: alt_sites[_] for _ in range(len(alt_sites))}
                    self.variants.append([pos, ref, alt, data[3:]])

    @staticmethod
    def write_file(variant_file, samples, variants, mode='raw'):
        with open(variant_file, 'w') as fout:
            fout.write("#POS\tREF\tALT\t%s\n" % ('\t'.join(sorted(samples))))
            for info in variants:
                pos, ref, alt, geno = info
                if mode == 'raw':
                    alt_list = [_ for _ in sorted(alt, key=lambda x: alt[x])]
                else:
                    alt_list = [alt[_] for _ in sorted(alt)]
                fout.write("%d\t%s\t%s\t%s\n" % (pos, ref, ','.join(alt_list), '\t'.join(map(str, geno))))
