class FastaIO:
    def __init__(self, fasta_file):
        self.__fasta_file = fasta_file
        self.fasta_db = {}
        self.seq_len_db = {}

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

        for gid in self.fasta_db:
            self.seq_len_db[gid] = len(self.fasta_db[gid])

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
        self.samples = []
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
                    alt = {_ + 1: alt_sites[_] for _ in range(len(alt_sites))}
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


class MatrixIO:
    def __init__(self):
        pass

    @staticmethod
    def read_merge_mat(mat_file):
        allele_db = {}
        with open(mat_file, 'r') as fin:
            start_read = False
            for line in fin:
                if line.strip() == "":
                    continue
                if line[0] == '#':
                    if line.startswith("#Allele Sequences"):
                        start_read = True
                        continue
                    if start_read:
                        data = line.strip().split()
                        allele_id = data[1]
                        allele_seq = data[2]
                        if allele_seq == 'Absence':
                            continue
                        gene_id = '-'.join(allele_id.split('-')[:-1])
                        if gene_id not in allele_db:
                            allele_db[gene_id] = {}
                        allele_db[gene_id][allele_id] = allele_seq
        return allele_db

    @staticmethod
    def write_merge_mat(merge_file, converted_var_db, var_cnt, cluster_sample_db, pheno_db):
        # converted_var_db is a dict:
        # gene_id:string => variant string 1:string => variant index 1:int
        # var_cnt is a int list
        # [gene1 variant type count, gene2 variant type count, ..., geneN variant type count]
        # cluster_sample_db is a dict
        # sample:string => gene_id:string => variant index:int
        # pheno_db is a dict:
        # sample:string => phenotype:float
        with open(merge_file, 'w') as fout:
            seq_db = {}
            fout.write("#TypeInfo\t%s\n" % ('\t'.join(map(str, var_cnt))))
            fout.write("#Sample")
            for gid in sorted(converted_var_db):
                for idx in range(len(converted_var_db[gid])):
                    fout.write("\t%s-Allele%d" % (gid, idx + 1))
                for seq in converted_var_db[gid]:
                    idx = converted_var_db[gid][seq]
                    seq_id = "%s-Allele%d" % (gid, idx)
                    seq_db[seq_id] = seq

            fout.write("\tPhenotype\n")

            for sample in sorted(cluster_sample_db):
                fout.write("%s" % sample)
                for gid in sorted(cluster_sample_db[sample]):
                    var_idx = cluster_sample_db[sample][gid]
                    var_info = [0 for _ in range(len(converted_var_db[gid]))]
                    var_info[var_idx - 1] = 1
                    fout.write("\t%s" % ('\t'.join(map(str, var_info))))
                fout.write("\t%f\n" % pheno_db[sample])

            fout.write("#\n#\n#Allele Sequences\n")
            for seq_id in sorted(seq_db):
                fout.write("# %s\t%s\n" % (seq_id, seq_db[seq_id] if seq_db[seq_id] else "Absence"))


class AlignIO:
    def __init__(self):
        pass

    @staticmethod
    def write_file(aln_file, aln_db):
        with open(aln_file, 'w') as fout:
            for gid in sorted(aln_db):
                fout.write(">%s\n%s\n" % (gid, aln_db[gid]))


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

    @staticmethod
    def read_asc(asc_file):
        var_db = {}
        with open(asc_file, 'r') as fin:
            for line in fin:
                if line[0] == '#':
                    continue
                data = line.strip().split()
                gid = data[0]
                pos = int(data[1])
                ref = data[-2]
                alt = data[-1].split(',')
                best_geno = data[-4].split('(')[0]
                if gid not in var_db:
                    var_db[gid] = {'pos': [pos], 'geno': [[_] for _ in data[2: -6]], 'ref': [ref], 'alt': [alt],
                                   'best_geno': [best_geno]}
                else:
                    var_db[gid]['pos'].append(pos)
                    # read genotype of each sample
                    for _ in range(len(data[2:-6])):
                        geno_idx = _ + 2
                        var_db[gid]['geno'][_].append(data[geno_idx])
                    var_db[gid]['ref'].append(ref)
                    var_db[gid]['alt'].append(alt)
                    var_db[gid]['best_geno'].append(best_geno)

        converted_var_db = {}
        for gid in var_db:
            converted_var_db[gid] = {'pos': var_db[gid]['pos'], 'geno': [], 'ref': var_db[gid]['ref'],
                                     'alt': var_db[gid]['alt'], 'best_geno': var_db[gid]['best_geno']}
            geno_set = set()
            # drop missing samples
            for _ in range(len(var_db[gid]['geno'])):
                if list(set(var_db[gid]['geno'][_])) == ['-']:
                    continue
                geno_set.add(tuple(var_db[gid]['geno'][_]))
            converted_var_db[gid]['geno'] = sorted(geno_set)
        return converted_var_db

    @staticmethod
    def read_asc_for_merge(asc_file):
        var_db = {}
        with open(asc_file, 'r') as fin:
            for line in fin:
                data = line.strip().split()
                if line[0] == '#':
                    samples = data[2:-6]
                else:
                    gid = data[0]
                    if gid not in var_db:
                        var_db[gid] = {}
                    for idx in range(len(samples)):
                        var = data[idx + 2]
                        sample = samples[idx]
                        if sample not in var_db[gid]:
                            var_db[gid][sample] = []
                        var_db[gid][sample].append(var)
        return var_db
