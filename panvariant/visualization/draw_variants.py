import matplotlib as mpl
import matplotlib.pyplot as plt
from bioplotz import multialign
from panvariant.io.file_operate import FastaIO, AssociateIO
from panvariant.io.message import Message as Msg
from os import listdir, path, makedirs
from pathos.multiprocessing import Pool


mpl.use('Agg')


def __draw_genes_in_association(cleanup_aln_dir, asc_file, pic_dir):
    if pic_dir and not path.exists(pic_dir):
        makedirs(pic_dir)

    var_sites_db = AssociateIO.read_asc(asc_file)

    for gene in var_sites_db:
        cleanup_aln_file = path.join(cleanup_aln_dir, "%s.aln" % gene)
        if not path.exists(cleanup_aln_file):
            Msg.warn("Gene not found, Aborting")
            return

        fasta_io = FastaIO(cleanup_aln_file)
        fasta_io.read_aln()
        aln_db = {}
        seq_set = set()
        seq_len = 0
        pic_file = path.join(pic_dir, "%s.pdf" % gene)

        # remove dup sequences
        allele_type = 1
        for smp in sorted(fasta_io.fasta_db):
            seq = fasta_io.fasta_db[smp]
            seq_len = len(seq)

            if seq not in seq_set:
                allele_name = "Allele %d" % allele_type
                allele_type += 1
                aln_db[allele_name] = seq
                seq_set.add(seq)

        if len(aln_db) <= 1:
            Msg.warn("\tToo less samples, Aborting")
            return

        highlight_pos = []
        for i in range(len(var_sites_db[gene]['pos'])):
            pos_start = var_sites_db[gene]['pos'][i]
            var_len = len(var_sites_db[gene]['ref'][i])
            for j in range(pos_start, pos_start+var_len):
                highlight_pos.append(j)

        Msg.info("\tPlotting %s" % cleanup_aln_file)
        fig_height = max(int(seq_len/200.*len(aln_db)/4), 5)
        plt.figure(figsize=(20, fig_height), dpi=100)
        plt.rcParams['font.sans-serif'] = 'Courier New'
        multialign(aln_db, base_per_line=200, highlight_positions=highlight_pos, highlight_color='blue',
                   match_color='lightgrey', mismatch_color='lightsalmon')
        plt.savefig(pic_file, bbox_inches='tight')
        plt.close('all')
        Msg.info("\tFinished")


def draw_variant_sites_in_association(cleanup_aln_dir, asc_dir, pic_dir, thread):
    pool = Pool(processes=thread)
    Msg.info("Visualizing allele genes with phenotypes")

    res = []
    for fn in listdir(asc_dir):
        asc_file = path.join(asc_dir, fn)
        sub_pic_dir = path.join(pic_dir, fn.replace('.asc', ''))
        res.append([fn, pool.apply_async(__draw_genes_in_association, (cleanup_aln_dir, asc_file, sub_pic_dir,))])
    pool.close()
    pool.join()

    for asc_fn, r in res:
        try:
            r.get()
        except Exception as e:
            Msg.warn("\tException caught with {}: {}". format(asc_fn, e))

    Msg.info("Visualization finished")
