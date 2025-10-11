from mate.io.file_operate import VariantIO, PhenoIO, AssociateIO
from mate.io.message import Message as Msg
from pathos.multiprocessing import Pool
from os import listdir, path
from scipy.stats import levene, ttest_ind
from outliers import smirnov_grubbs as grubbs
from numpy import array, average, std
import warnings

warnings.filterwarnings("error")


def __associate_with_single_pheno(pheno_file, cla_dir, asc_file, is_lower_better):
    Msg.info("\tLoading phenotypes")
    pheno = PhenoIO(pheno_file)
    pheno.read_pheno()

    # full_info is a dictionary like below:
    # gene => 'samples' => sample_list
    #      => 'variants' => variants
    #      => 'stats' => stat_info
    full_info = {}
    samples_set = set()

    is_empty = True
    for fn in sorted(listdir(cla_dir)):
        gene = '.'.join(fn.split('.')[:-1])
        cla_file = path.join(cla_dir, fn)
        Msg.info("\tLoading %s" % cla_file)
        var_io = VariantIO()
        var_io.read_var(cla_file)
        if not var_io.variants:
            Msg.info("\tNo variant found, skipping...")
            continue

        full_info[gene] = {'samples': var_io.samples}
        for s in var_io.samples:
            samples_set.add(s)
        Msg.info("\tComparing best and second best variant type")

        variants = []
        stat_info = []
        for i in range(len(var_io.variants)):
            pos, ref, alt, geno = var_io.variants[i]
            cur_site_pheno_db = {}
            for geno_idx in range(len(geno)):
                if var_io.samples[geno_idx] not in pheno.pheno_db:
                    continue
                if geno[geno_idx] not in cur_site_pheno_db:
                    cur_site_pheno_db[geno[geno_idx]] = []
                cur_site_pheno_db[geno[geno_idx]].append(pheno.pheno_db[var_io.samples[geno_idx]])

            # # remove outliers with grubbs test and get best and second best pheno list for comparison
            cur_site_pheno_list = []
            for pheno_idx in cur_site_pheno_db:
                try:
                    clean_pheno = list(grubbs.test(array(cur_site_pheno_db[pheno_idx]), alpha=0.05))

                except RuntimeWarning:
                    clean_pheno = cur_site_pheno_db[pheno_idx]
                cur_site_pheno_list.append([average(clean_pheno), pheno_idx, clean_pheno])
            if len(cur_site_pheno_list) < 2:
                continue
            cur_site_pheno_list = sorted(cur_site_pheno_list, reverse=True if not is_lower_better else False)
            _, best_type, best_pheno = cur_site_pheno_list[0]
            _, sec_best_type, sec_best_pheno = cur_site_pheno_list[1]

            if len(best_pheno) <= 1 or len(sec_best_pheno) <= 1:
                continue
            try:
                _, levene_pvalue = levene(best_pheno, sec_best_pheno)
            except RuntimeWarning:
                continue
            if levene_pvalue > 0.05:
                equal_var = True
            else:
                equal_var = False
            try:
                t_pvalue = ttest_ind(best_pheno, sec_best_pheno, equal_var=equal_var).pvalue
            except RuntimeWarning:
                continue
            if t_pvalue <= 0.05:
                variants.append(var_io.variants[i])
                stat_info.append([levene_pvalue, t_pvalue,
                                  "%d(%d|%.2f|%.2f)" % (best_type, len(best_pheno),
                                                        average(best_pheno), std(best_pheno)),
                                  "%d(%d|%.2f|%.2f)" % (sec_best_type, len(sec_best_pheno),
                                                        average(sec_best_pheno), std(sec_best_pheno))])
        if variants:
            is_empty = False
        full_info[gene]['variants'] = variants
        full_info[gene]['stats'] = stat_info

    if not is_empty:
        asc_io = AssociateIO()
        Msg.info("\tWriting %s" % asc_file)
        asc_io.write_file(asc_file, samples_set, full_info)
    else:
        Msg.info("\tNo association sites found, skipping...")
    Msg.info("\tFinished")


def associate_with_pheno(pheno_dir, cla_dir, asc_dir, thread):
    pool = Pool(processes=thread)
    Msg.info("Associating with phenotypes")

    res = []
    for pheno_fn in listdir(pheno_dir):
        Msg.info("\tAssociating with %s" % pheno_fn)
        if pheno_fn.startswith("LOW-"):
            is_lower_better = True
        else:
            is_lower_better = False

        pheno_file = path.join(pheno_dir, pheno_fn)
        asc_file = path.join(asc_dir, '.'.join(pheno_fn.split('.')[:-1]) + '.asc')
        res.append([pheno_fn, pool.apply_async(__associate_with_single_pheno,
                                               (pheno_file, cla_dir, asc_file, is_lower_better,))])
    pool.close()
    pool.join()

    # If subprocess failed, the error will be caught.
    for pheno_fn, r in res:
        try:
            r.get()
        except Exception as e:
            Msg.warn("\tException caught with {}: {}".format(pheno_fn, e))

    Msg.info("Association finished")
