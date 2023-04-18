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
                    self.variants.append([pos, ref, alt, list(map(int, data[3:]))])

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


class PhenoIO:
    def __init__(self, pheno_file):
        self.__pheno_file = pheno_file
        self.pheno_db = {}

    def read_pheno(self):
        with open(self.__pheno_file, 'r') as fin:
            for line in fin:
                data = line.strip().split()
                self.pheno_db[data[0]] = float(data[1])


class AssociateIO:
    def __init__(self):
        pass

    @staticmethod
    def write_file(asc_file, samples, full_info):
        sample_list = sorted(samples)
        with open(asc_file, 'w') as fout:
            fout.write("#Gene\tPOS\t%s\tLevene_pvalue\tTtest_pvalue\t"
                       "Best_info(Count|Avg|Std)\tSecond_best_info(Count|Avg|Std)\tREF\tALT\n"
                       % ('\t'.join(sample_list)))
            for gene in sorted(full_info):
                cur_samples = full_info[gene]['samples']
                for idx in range(len(full_info[gene]['variants'])):
                    pos, ref, alt, geno = full_info[gene]['variants'][idx]
                    cur_sample_geno = {cur_samples[_]: geno[_] for _ in range(len(cur_samples))}
                    cur_sample_set = set(cur_samples)
                    full_geno = ['-' for _ in range(len(sample_list))]
                    for sample_idx in range(len(sample_list)):
                        sample = sample_list[sample_idx]
                        if sample in cur_sample_set:
                            full_geno[sample_idx] = str(cur_sample_geno[sample])
                    levene_pvalue, ttest_pvalue, best_info, sec_best_info = full_info[gene]['stats'][idx]
                    alt_list = [alt[_] for _ in sorted(alt)]
                    fout.write("%s\t%d\t%s\t%.2f\t%.2f\t%s\t%s\t%s\t%s\n"
                               % (gene, pos, '\t'.join(full_geno), levene_pvalue, ttest_pvalue,
                                  best_info, sec_best_info, ref, ','.join(alt_list)))
