from panvariant.io.file_operate import VariantIO, PhenoIO, AssociateIO
from panvariant.io.message import Message as Msg
from pathos.multiprocessing import Pool
from os import listdir, path


def __merge_with_single_pheno(pheno_file, asc_file, merge_file):
    Msg.info("\tLoading phenotypes")
    pheno = PhenoIO(pheno_file)
    pheno.read_pheno()

    Msg.info("\tLoading association file")
    asc = AssociateIO()
    var_db = asc.read_asc_for_merge(asc_file)

    Msg.info("Converting data")
    converted_var_db = {}
    cluster_sample_db = {}
    for gid in var_db:
        converted_var_db[gid] = {}
        var_idx = 1
        for sample in var_db[gid]:
            if sample not in pheno.pheno_db:
                continue
            var = ''.join(var_db[gid][sample])
            if var not in converted_var_db[gid]:
                converted_var_db[gid][var] = var_idx
                var_idx += 1

            if sample not in cluster_sample_db:
                cluster_sample_db[sample] = {}
            cluster_sample_db[sample][gid] = converted_var_db[gid][var]

    var_cnt = []
    for gid in sorted(converted_var_db):
        var_cnt.append(len(converted_var_db[gid]))

    Msg.info("\tWriting merged matrix")
    var_io = VariantIO()
    var_io.write_merge_mat(merge_file, converted_var_db, var_cnt, cluster_sample_db, pheno.pheno_db)
    Msg.info("Finished")


def final_merge_variants(pheno_dir, asc_dir, merge_dir, thread):
    pool = Pool(processes=thread)
    Msg.info("Merging data")

    res = []
    for pheno_fn in listdir(pheno_dir):
        Msg.info("\tMerging with %s" % pheno_fn)
        pheno_file = path.join(pheno_dir, pheno_fn)
        asc_file = path.join(asc_dir, '.'.join(pheno_fn.split('.')[:-1]) + '.asc')
        merge_file = path.join(merge_dir, '.'.join(pheno_fn.split('.')[:-1]) + '.mat')
        res.append([pheno_fn, pool.apply_async(__merge_with_single_pheno, (pheno_file, asc_file, merge_file, ))])
    pool.close()
    pool.join()

    # If subprocess failed, the error will be caught.
    for pheno_fn, r in res:
        try:
            r.get()
        except Exception as e:
            Msg.warn("\tException caught with {}: {}".format(pheno_fn, e))

    Msg.info("Association finished")
