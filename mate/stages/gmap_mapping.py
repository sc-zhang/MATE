from os import listdir, path, getcwd, chdir
from mate.base.dep_check import DepCheck
from mate.base.runner import Runner
from mate.io.message import Message as Msg
from pathos.multiprocessing import Pool


def __para_gmap(ref_cds, genome_dir, genome_fn, gmap_out_gff3, ploidy, thread):
    genome_file = path.join(genome_dir, genome_fn)
    runner = Runner()
    ferr = open(gmap_out_gff3 + '.log', 'w')
    Msg.info("\tBuilding GMAP database %s" % genome_file)
    runner.set_cmd("gmap_build -D . -d %s_DB -t %d %s" % (genome_fn, thread, genome_file))
    runner.run()
    ferr.write("%s\n" % runner.get_err())

    Msg.info("\tMapping %s" % genome_file)
    fsize = path.getsize(genome_file)
    if fsize > 2 ** 32:
        gmap_program = "gmapl"
    else:
        gmap_program = "gmap"
    runner.set_cmd("%s -D . -d %s_DB -f 2 -n %d -t %d %s > %s" %
                   (gmap_program, genome_fn, ploidy, thread, ref_cds, gmap_out_gff3))
    runner.run()
    ferr.write("%s\n" % runner.get_err())
    ferr.close()


def mapping(ref_cds, genome_dir, gmap_out_gff3_dir, ploidy, thread):
    Msg.info("Checking dependencies")
    if (not DepCheck.check("gmap_build")) or (not DepCheck.check("gmap")):
        Msg.error("GMAP not found, please install it and add to the PATH environment variable")
        exit(-1)

    Msg.info("Pass")
    Msg.info("Running GMAP on pan-genome")
    ref_cds = path.abspath(ref_cds)
    genome_dir = path.abspath(genome_dir)
    cur_dir = path.abspath(getcwd())

    Msg.info("Entering %s" % gmap_out_gff3_dir)
    chdir(gmap_out_gff3_dir)
    genome_cnt = len(listdir(genome_dir))
    main_thread = min(genome_cnt, thread)
    sub_thread = max(thread // main_thread, 1)
    pool = Pool(processes=main_thread)
    res = []
    for fn in listdir(genome_dir):
        gmap_out_gff3 = path.join(gmap_out_gff3_dir, '.'.join(fn.split('.')[:-1]) + ".gff3")
        res.append([fn, pool.apply_async(__para_gmap,
                                         (ref_cds, genome_dir, fn, gmap_out_gff3, ploidy, sub_thread,))])
    pool.close()
    pool.join()

    for fn, r in res:
        try:
            r.get()
        except Exception as e:
            Msg.warn("\tException caught with {}: {}".format(fn, e))

    Msg.info("GMAP finished")
    chdir(cur_dir)
