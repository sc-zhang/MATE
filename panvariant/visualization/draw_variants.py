import matplotlib as mpl
import matplotlib.pyplot as plt
from bioplotz import multialign
from panvariant.io.file_operate import FastaIO
from panvariant.io.message import Message as Msg
from os import listdir, path
from pathos.multiprocessing import Pool


mpl.use('Agg')


def __draw_genes_in_association(cleanup_aln_file, pic_file):
    fasta_io = FastaIO(cleanup_aln_file)
    fasta_io.read_aln()
    aln_db = {}
    seq_set = set()
    seq_len = 0
    for smp in sorted(fasta_io.fasta_db):
        seq = fasta_io.fasta_db[smp]
        seq_len = len(seq)
        if seq not in seq_set:
            aln_db[smp] = seq
            seq_set.add(seq)
    if len(aln_db) <= 1:
        Msg.info("\tToo less samples, Aborting")
        return

    Msg.info("\tPlotting %s" % cleanup_aln_file)
    fig_height = max(int(seq_len/200.*len(aln_db)/5), 2)
    plt.figure(figsize=(20, fig_height), dpi=100)
    multialign(aln_db, base_per_line=200)
    plt.rcParams['font.sans-serif'] = 'Courier New'
    plt.savefig(pic_file, bbox_inches='tight')
    Msg.info("\tFinished")


def draw_variant_sites_in_association(cleanup_aln_dir, pic_dir, thread):
    pool = Pool(processes=thread)
    Msg.info("Visualizing allele genes with phenotypes")

    res = []
    for fn in listdir(cleanup_aln_dir):
        cleanup_aln_file = path.join(cleanup_aln_dir, fn)
        pic_file = path.join(pic_dir, fn.replace('.aln', '.pdf'))
        res.append([fn, pool.apply_async(__draw_genes_in_association, (cleanup_aln_file, pic_file,))])
    pool.close()
    pool.join()

    for asc_fn, r in res:
        try:
            r.get()
        except Exception as e:
            Msg.warn("\tException caught with {}: {}". format(asc_fn, e))

    Msg.info("Visualization finished")
