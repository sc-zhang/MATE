import copy

from mate.io.file_operate import PhenoIO, FastaIO, MatrixIO, AssociateIO
from mate.io.message import Message as Msg
from pathos.multiprocessing import Pool
from os import listdir, path


def __get_sig_vars(seq, gene_var_sites_db):
    sig_vars = ['=' for _ in range(len(seq))]
    for i in range(len(gene_var_sites_db['pos'])):
        pos_start = gene_var_sites_db['pos'][i]
        var_len = len(gene_var_sites_db['ref'][i])
        for j in range(pos_start, pos_start+var_len):
            sig_vars[j] = seq[j]
    return ''.join(sig_vars)


def __merge_with_single_pheno(pheno_file, cleanup_aln_dir, asc_file, merge_cleanup_file, merge_sig_file):
    Msg.info("\tLoading phenotypes")
    pheno = PhenoIO(pheno_file)
    pheno.read_pheno()

    Msg.info("Converting data")

    converted_cleanup_allele_db = {}
    converted_sig_allele_db = {}
    cluster_sample_cleanup_db = {}
    cluster_sample_sig_db = {}
    cleanup_allele_idx = 1
    sig_allele_idx = 1

    for fn in listdir(cleanup_aln_dir):
        is_sig = True
        cleanup_aln_file = path.join(cleanup_aln_dir, fn)
        gid = '.'.join(fn.split('.')[:-1])
        if not path.exists(asc_file):
            var_sites_db = {}
        else:
            var_sites_db = AssociateIO.read_asc(asc_file)
        if gid not in var_sites_db:
            is_sig = False
        fasta_io = FastaIO(cleanup_aln_file)
        fasta_io.read_aln()
        converted_cleanup_allele_db[gid] = {}
        cleanup_allele_idx = 1
        if is_sig:
            converted_sig_allele_db[gid] = {}
            sig_allele_idx = 1

        # for one allele which supported samples less than 5% or more than 95% samples, mark it as absence
        smp_cnt = len(pheno.pheno_db)
        lower_threshold = smp_cnt*.05
        upper_threshold = smp_cnt*.95
        support_cnt = {}
        sig_support_cnt = {}
        for smp in sorted(pheno.pheno_db):
            if smp not in fasta_io.fasta_db:
                continue
            else:
                seq = fasta_io.fasta_db[smp]
            if seq not in support_cnt:
                support_cnt[seq] = 0
            support_cnt[seq] += 1
            if is_sig:
                sig_vars = __get_sig_vars(seq, var_sites_db[gid])
                if sig_vars not in sig_support_cnt:
                    sig_support_cnt[sig_vars] = 0
                sig_support_cnt[sig_vars] += 1

        for smp in sorted(pheno.pheno_db):
            if smp not in fasta_io.fasta_db:
                continue
            else:
                seq = fasta_io.fasta_db[smp]
                if not (lower_threshold < support_cnt[seq] < upper_threshold):
                    continue
                if seq not in converted_cleanup_allele_db[gid]:
                    converted_cleanup_allele_db[gid][seq] = cleanup_allele_idx
                    cleanup_allele_idx += 1
                if is_sig:
                    sig_vars = __get_sig_vars(seq, var_sites_db[gid])
                    if not (lower_threshold < sig_support_cnt[sig_vars] < upper_threshold):
                        continue
                    if sig_vars not in converted_sig_allele_db[gid]:
                        converted_sig_allele_db[gid][sig_vars] = sig_allele_idx
                        sig_allele_idx += 1

        # for absence, the allele idx set to last allele
        converted_cleanup_allele_db[gid][""] = cleanup_allele_idx
        if is_sig:
            converted_sig_allele_db[gid][""] = sig_allele_idx

        for smp in sorted(pheno.pheno_db):
            if smp not in cluster_sample_cleanup_db:
                cluster_sample_cleanup_db[smp] = {}
            if smp in fasta_io.fasta_db:
                seq = fasta_io.fasta_db[smp]
            else:
                seq = ""
            if seq not in converted_cleanup_allele_db[gid]:
                seq = ""
            cluster_sample_cleanup_db[smp][gid] = converted_cleanup_allele_db[gid][seq]
            if is_sig:
                if smp not in cluster_sample_sig_db:
                    cluster_sample_sig_db[smp] = {}
                if seq:
                    sig_vars = __get_sig_vars(seq, var_sites_db[gid])
                else:
                    sig_vars = ""
                if sig_vars not in converted_sig_allele_db[gid]:
                    sig_vars = ""
                cluster_sample_sig_db[smp][gid] = converted_sig_allele_db[gid][sig_vars]

        missing_cnt = 0
        for smp in cluster_sample_cleanup_db:
            if cluster_sample_cleanup_db[smp][gid] == converted_cleanup_allele_db[gid][""]:
                missing_cnt += 1
        if missing_cnt >= upper_threshold:
            converted_cleanup_allele_db.pop(gid)
            for smp in cluster_sample_cleanup_db:
                cluster_sample_cleanup_db[smp].pop(gid)

        missing_cnt = 0
        if is_sig:
            for smp in cluster_sample_sig_db:
                if cluster_sample_sig_db[smp][gid] == converted_sig_allele_db[gid][""]:
                    missing_cnt += 1
            if missing_cnt >= len(cluster_sample_sig_db)*.95:
                converted_sig_allele_db.pop(gid)
            for smp in cluster_sample_sig_db:
                cluster_sample_sig_db[smp].pop(gid)

    Msg.info("\tWriting merged matrix")
    allele_io = MatrixIO()
    allele_cnt = [len(converted_cleanup_allele_db[gid]) for gid in sorted(converted_cleanup_allele_db)]
    allele_io.write_merge_mat(merge_cleanup_file, converted_cleanup_allele_db,
                              allele_cnt, cluster_sample_cleanup_db, pheno.pheno_db)
    allele_cnt = [len(converted_sig_allele_db[gid]) for gid in sorted(converted_sig_allele_db)]
    allele_io.write_merge_mat(merge_sig_file, converted_sig_allele_db,
                              allele_cnt, cluster_sample_sig_db, pheno.pheno_db)
    Msg.info("Finished")


def merge_variant_matrix(pheno_dir, cleanup_aln_dir, asc_dir, merge_cleanup_dir, merge_sig_dir, thread):
    pool = Pool(processes=thread)
    Msg.info("Merging data")

    res = []
    for pheno_fn in listdir(pheno_dir):
        Msg.info("\tMerging with %s" % pheno_fn)
        pheno_file = path.join(pheno_dir, pheno_fn)
        asc_file = path.join(asc_dir, '.'.join(pheno_fn.split('.')[:-1]) + '.asc')
        merge_cleanup_file = path.join(merge_cleanup_dir, '.'.join(pheno_fn.split('.')[:-1]) + '.mat')
        merge_sig_file = path.join(merge_sig_dir, '.'.join(pheno_fn.split('.')[:-1]) + '.mat')
        res.append([pheno_fn, pool.apply_async(__merge_with_single_pheno,
                                               (pheno_file, cleanup_aln_dir, asc_file,
                                                merge_cleanup_file, merge_sig_file, ))])
    pool.close()
    pool.join()

    # If subprocess failed, the error will be caught.
    for pheno_fn, r in res:
        try:
            r.get()
        except Exception as e:
            Msg.warn("\tException caught with {}: {}".format(pheno_fn, e))

    Msg.info("Association finished")
