from mate.io.file_operate import VariantIO
from mate.io.message import Message as Msg
from pathos.multiprocessing import Pool
from os import listdir, path


def __variant_classifier_for_single_file(var_file, classified_file):
    Msg.info("\tLoading %s" % var_file)
    var_io = VariantIO()
    var_io.read_var(var_file)
    # var_io.variants is 2-D list, like:
    # [[pos, ref, alt(dict), [0, 1, 1, 0, ...](sample genotype)
    # ...
    # ]
    last_pos_range = [-1, -1]
    last_ref = ""
    last_alt = {}
    last_genotype = []
    full_info = []

    Msg.info("\tMerging")
    for _ in range(len(var_io.variants)):
        if _ == 0:
            last_pos, last_ref, last_alt, last_genotype = var_io.variants[_]
            last_pos_range[0] = last_pos
            last_pos_range[1] = last_pos
            continue
        cur_pos, cur_ref, cur_alt, cur_genotype = var_io.variants[_]
        if cur_genotype == last_genotype and cur_pos == last_pos_range[1]+1:
            last_ref += cur_ref
            for _ in last_alt:
                alt = last_alt[_] + cur_alt[_]
                last_alt[_] = alt
            last_pos_range[1] = cur_pos
        else:
            full_info.append([last_pos_range[0], last_ref, last_alt, last_genotype])
            last_pos_range = [cur_pos, cur_pos]
            last_ref = cur_ref
            last_alt = cur_alt
            last_genotype = cur_genotype
    if last_genotype:
        full_info.append([last_pos_range[0], last_ref, last_alt, last_genotype])

    var_io.write_file(classified_file, var_io.samples, full_info, mode='cla')

    Msg.info("\tFinished")


def variant_classifier(var_dir, cla_dir, thread):
    pool = Pool(processes=thread)
    Msg.info("Variant classifying")

    res = []
    for fn in listdir(var_dir):
        Msg.info("\tClassifying %s" % fn)
        var_file = path.join(var_dir, fn)
        cla_file = path.join(cla_dir, fn.replace('.var', '.cla'))
        res.append([fn, pool.apply_async(__variant_classifier_for_single_file, (var_file, cla_file, ))])
    pool.close()
    pool.join()

    # If subprocess failed, the error will be caught.
    for var_fn, r in res:
        try:
            r.get()
        except Exception as e:
            Msg.warn("\tException caught with {}: {}".format(var_fn, e))

    Msg.info("Variant classified")
