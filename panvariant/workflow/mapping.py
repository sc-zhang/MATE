from os import listdir, path, getcwd, chdir
from panvariant.operate.dep_check import DepCheck
from panvariant.operate.runner import Runner
from panvariant.io.message import Message as Msg


def mapping(ref_cds, genome_dir, gmap_out_gff3_dir, ploidy, thread):
    Msg.info("Checking dependencies")
    if (not DepCheck.check("gmap_build")) or (not DepCheck.check("gmap")):
        Msg.error("GMAP not found, please install it and add to the PATH environment variable")
        exit(-1)

    Msg.info("Pass")

    runner = Runner()
    Msg.info("Running GMAP on pan-genome")
    ref_cds = path.abspath(ref_cds)
    genome_dir = path.abspath(genome_dir)
    cur_dir = path.abspath(getcwd())
    
    Msg.info("Entering %s" % gmap_out_gff3_dir)
    chdir(gmap_out_gff3_dir)
    for fn in listdir(genome_dir):
        genome_file = path.join(genome_dir, fn)
        gff3_file = path.join(gmap_out_gff3_dir, '.'.join(fn.split('.')[:-1])+".gff3")
        Msg.info("\tBuild %s" % fn)
        runner.set_cmd("gmap_build -D . -d %s_DB -t %d %s" % (fn, thread, genome_file))
        runner.run()
        Msg.info("\tGmap %s" % fn)
        runner.set_cmd("gmap -D . -d %s_DB -f %d -n 1 -t %d %s > %s" % (fn, ploidy, thread, ref_cds, gff3_file))
        runner.run()
    Msg.info("GMAP finished")
    chdir(cur_dir)
