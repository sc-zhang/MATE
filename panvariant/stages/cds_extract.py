from os import path, listdir
from panvariant.io.file_operate import FastaIO
from panvariant.io.message import Message as Msg
from pathos.multiprocessing import Pool


class CDSExtract:
    def __init__(self, gmap_out_gff3_dir, genome_dir, out_cds_dir, thread):
        self.__gmap_out_gff3_dir = path.abspath(gmap_out_gff3_dir)
        self.__genome_dir = path.abspath(genome_dir)
        self.__out_cds_dir = path.abspath(out_cds_dir)
        self.__thread = thread

    @staticmethod
    def __reverse_seq(seq):
        rev_seq = ""
        base_db = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        for base in seq[::-1]:
            if base in base_db:
                rev_seq += base_db[base]
            else:
                rev_seq += base
        return rev_seq

    def __extract_cds(self, in_gff3, in_fa, out_cds_file):
        Msg.info("\tLoading genome %s" % in_fa)
        fasta_io = FastaIO(in_fa)
        fasta_io.read_fasta()

        Msg.info("\tExtracting cds")
        with open(in_gff3, 'r') as fin:
            with open(out_cds_file, 'w') as fout:
                region_db = {}
                for line in fin:
                    if line.strip() == '' or line[0] == '#':
                        continue
                    data = line.strip().split()
                    if data[2] == 'gene':
                        for info in data[8].split(';'):
                            if info.startswith("Name"):
                                gid = info.split('=')[1]
                        chrn = data[0]
                        if gid not in region_db:
                            region_db[gid] = {'chrn': chrn, 'rna': [], 'best_score': 0, 'direct': "", 'cds': []}

                    elif data[2] == 'mRNA':
                        score = 1
                        for info in data[8].split(';'):
                            if info.startswith('cov') or info.startswith('iden'):
                                score *= float(info.split('=')[1])
                        if score > region_db[gid]['best_score']:
                            region_db[gid]['chrn'] = data[0]
                            region_db[gid]['best_score'] = score
                            region_db[gid]['rna'] = [int(data[3]), int(data[4])]
                            region_db[gid]['direct'] = data[6]
                            region_db[gid]['cds'] = []
                    elif data[2] == 'CDS':
                        if score == region_db[gid]['best_score']:
                            region_db[gid]['cds'].append([int(data[3]), int(data[4])])

                for gid in sorted(region_db):
                    chrn = region_db[gid]['chrn']
                    if region_db[gid]['direct'] == '+':
                        cds = ""
                        for sp, ep in region_db[gid]['cds']:
                            cds += fasta_io.fasta_db[chrn][sp - 1: ep]
                        fout.write('>%s\n%s\n' % (gid, cds))
                    else:
                        cds = ""
                        for sp, ep in region_db[gid]['cds']:
                            cds += self.__reverse_seq(fasta_io.fasta_db[chrn][sp - 1: ep])
                        fout.write('>%s\n%s\n' % (gid, cds))

        Msg.info("\tCDS extracted")

    def extract(self):
        Msg.info("Extracting")
        pool = Pool(processes=self.__thread)

        res = []
        for fn in listdir(self.__genome_dir):
            sample_id = '.'.join(fn.split('.')[:-1])
            genome_file = path.join(self.__genome_dir, fn)
            gff3_file = path.join(self.__gmap_out_gff3_dir, "%s.gff3" % sample_id)
            out_cds = path.join(self.__out_cds_dir, "%s.cds" % sample_id)
            res.append([fn, pool.apply_async(self.__extract_cds, (gff3_file, genome_file, out_cds,))])
        pool.close()
        pool.join()

        # If subprocess failed, the error will be caught.
        for genome_fn, r in res:
            try:
                r.get()
            except Exception as e:
                Msg.warn("\tException caught with {}: {}".format(genome_fn, e))

        Msg.info("Extract finished")
