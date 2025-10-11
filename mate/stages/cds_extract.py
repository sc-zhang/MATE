from os import path, listdir
from mate.io.file_operate import FastaIO, Gff3IO
from mate.io.message import Message as Msg
from pathos.multiprocessing import Pool


class CDSExtract:
    def __init__(self, gmap_out_gff3_dir, genome_dir, out_cds_dir, thread):
        self.__gmap_out_gff3_dir = path.abspath(gmap_out_gff3_dir)
        self.__genome_dir = path.abspath(genome_dir)
        self.__out_cds_dir = path.abspath(out_cds_dir)
        self.__thread = thread

    @staticmethod
    def __extract_cds(in_gff3, in_fa, out_cds_file):
        Msg.info("\tLoading genome %s" % in_fa)
        fasta_io = FastaIO(in_fa)
        fasta_io.read_fasta()

        Msg.info("\tExtracting cds")
        gff3_io = Gff3IO(in_gff3)
        gff3_io.read_gff3()
        gff3_io.extract_cds(fasta_io.fasta_db)
        gff3_io.write_cds(out_cds_file)

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
