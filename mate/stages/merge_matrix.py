from mate.io.file_operate import PhenoIO, FastaIO, MatrixIO
from mate.io.message import Message as Msg
from pathos.multiprocessing import Pool
from os import listdir, path


def __merge_with_single_pheno(pheno_file, cleanup_aln_dir, merge_file):
    Msg.info("\tLoading phenotypes")
    pheno = PhenoIO(pheno_file)
    pheno.read_pheno()

    Msg.info("Converting data")

    converted_allele_db = {}
    cluster_sample_db = {}
    for fn in listdir(cleanup_aln_dir):
        cleanup_aln_file = path.join(cleanup_aln_dir, fn)
        gid = fn.split('.')[0]
        fasta_io = FastaIO(cleanup_aln_file)
        fasta_io.read_aln()
        converted_allele_db[gid] = {}
        allele_idx = 1
        for smp in sorted(pheno.pheno_db):
            if smp not in fasta_io.fasta_db:
                continue
            else:
                seq = fasta_io.fasta_db[smp]
                if seq not in converted_allele_db[gid]:
                    converted_allele_db[gid][seq] = allele_idx
                    allele_idx += 1

        # for absence, the allele idx set to last allele
        converted_allele_db[gid][""] = allele_idx
        for smp in sorted(pheno.pheno_db):
            if smp not in cluster_sample_db:
                cluster_sample_db[smp] = {}
            if smp in fasta_io.fasta_db:
                seq = fasta_io.fasta_db[smp]
            else:
                seq = ""
            cluster_sample_db[smp][gid] = converted_allele_db[gid][seq]

    allele_cnt = []
    for gid in sorted(converted_allele_db):
        allele_cnt.append(len(converted_allele_db[gid]))

    Msg.info("\tWriting merged matrix")
    allele_io = MatrixIO()
    allele_io.write_merge_mat(merge_file, converted_allele_db, allele_cnt, cluster_sample_db, pheno.pheno_db)
    Msg.info("Finished")


def merge_variant_matrix(pheno_dir, cleanup_aln_dir, merge_dir, thread):
    pool = Pool(processes=thread)
    Msg.info("Merging data")

    res = []
    for pheno_fn in listdir(pheno_dir):
        Msg.info("\tMerging with %s" % pheno_fn)
        pheno_file = path.join(pheno_dir, pheno_fn)
        merge_file = path.join(merge_dir, '.'.join(pheno_fn.split('.')[:-1]) + '.mat')
        res.append([pheno_fn, pool.apply_async(__merge_with_single_pheno, (pheno_file, cleanup_aln_dir, merge_file, ))])
    pool.close()
    pool.join()

    # If subprocess failed, the error will be caught.
    for pheno_fn, r in res:
        try:
            r.get()
        except Exception as e:
            Msg.warn("\tException caught with {}: {}".format(pheno_fn, e))

    Msg.info("Association finished")
