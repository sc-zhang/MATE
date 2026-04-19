from mate.base.dep_check import DepCheck
from mate.base.runner import Runner
from mate.io.message import Message as Msg
from mate.io.file_operate import Converter
from os import path, makedirs, listdir


def mafft_alignment(fasta_file_dir_for_mafft, aln_dir, thread):
    Msg.info("Checking mafft")
    if not DepCheck.check("which mafft"):
        Msg.error("mafft not found, please install it and add to the PATH environment variable")
        exit(-1)
    Msg.info("Pass")

    Msg.info("Running mafft")
    if not path.exists(aln_dir):
        makedirs(aln_dir)
    runner = Runner()
    converter = Converter()
    for fn in listdir(fasta_file_dir_for_mafft):
        if not fn.endswith(".pep"):
            continue
        in_fa = path.join(fasta_file_dir_for_mafft, fn)
        out_aln = path.join(aln_dir, '.'.join(fn.split('.')[:-1]) + '.aln_pep')
        Msg.info("\tmafft %s" % fn)
        runner.set_cmd("mafft --adjustdirection --thread %d %s > %s" % (thread, in_fa, out_aln))
        runner.run()
        out_cds_aln = out_aln.replace(".aln_pep", ".aln")
        converter.restore_pep_aln_to_cds_aln(out_aln, in_fa.replace(".pep", ".cds"), out_cds_aln)

    Msg.info("Mafft finished")
