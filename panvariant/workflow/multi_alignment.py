from panvariant.operate.dep_check import DepCheck
from panvariant.operate.runner import Runner
from panvariant.io.message import Message as Msg
from os import path, makedirs


def mafft_alignment(fasta_file_dir_for_mafft, aln_dir, thread):
    Msg.info("Checking mafft")
    if not DepCheck.check("mafft"):
        Msg.error("mafft not found, please install it and add to the PATH environment variable")
        exit(-1)
    Msg.info("Pass")

    Msg.info("Running mafft")
    if not path.exists(aln_dir):
        makedirs(aln_dir)
    runner = Runner()
    for fn in fasta_file_dir_for_mafft:
        in_fa = path.join(fasta_file_dir_for_mafft, fn)
        out_aln = path.join(aln_dir, '.'.join(fn.split('.')[:-1])+'.aln')
        Msg.info("\tmafft %s" % fn)
        runner.set_cmd("mafft --adjustdirection --thread %d %s > %s" % (thread, in_fa, out_aln))
        runner.run()
    Msg.info("Mafft finished")
