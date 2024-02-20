import matplotlib as mpl
import matplotlib.pyplot as plt
from bioplotz import multialign
from mate.io.file_operate import MatrixIO
from mate.io.message import Message as Msg
from os import listdir, path, makedirs
from pathos.multiprocessing import Pool


mpl.use('Agg')


def __draw_sig_alleles(mat_file, pic_dir):
    if pic_dir and not path.exists(pic_dir):
        makedirs(pic_dir)

    allele_io = MatrixIO()
    allele_db = allele_io.read_merge_mat(mat_file)
    for gene_id in allele_db:
        Msg.info("\tPlotting %s" % mat_file)
        seq_len = 0
        for allele_id in allele_db[gene_id]:
            seq_len = len(allele_db[gene_id][allele_id])
            break

        fig_height = max(int(seq_len/200.*len(allele_db[gene_id])/4), 5)
        plt.figure(figsize=(20, fig_height), dpi=100)
        plt.rcParams['font.sans-serif'] = 'Courier New'
        multialign(allele_db[gene_id],
                   base_per_line=200,
                   match_color='lightgrey',
                   mismatch_color='blue',
                   mismatch_background_color='lightskyblue')
        pic_file = path.join(pic_dir, "%s.pdf" % gene_id)
        plt.savefig(pic_file, bbox_inches='tight')
        plt.close('all')
        Msg.info("\tFinished")


def draw_sig_alleles(mat_dir, pic_dir, thread):
    pool = Pool(processes=thread)
    Msg.info("Visualizing allele genes with phenotypes")

    res = []
    for fn in listdir(mat_dir):
        mat_file = path.join(mat_dir, fn)
        sub_pic_dir = path.join(pic_dir, fn.replace('.mat', ''))
        res.append([fn, pool.apply_async(__draw_sig_alleles, (mat_file, sub_pic_dir,))])
    pool.close()
    pool.join()

    for asc_fn, r in res:
        try:
            r.get()
        except Exception as e:
            Msg.warn("\tException caught with {}: {}". format(asc_fn, e))

    Msg.info("Visualization finished")


def __draw_cleanup_alleles(gene_id, gene_allele_db, pic_dir):
    Msg.info("\tPlotting %s" % gene_id)
    seq_len = 0
    for allele_id in gene_allele_db:
        seq_len = len(gene_allele_db[allele_id])
        break

    fig_height = max(int(seq_len / 200. * len(gene_allele_db) / 4), 5)
    plt.figure(figsize=(20, fig_height), dpi=100)
    plt.rcParams['font.sans-serif'] = 'Courier New'
    multialign(gene_allele_db,
               base_per_line=200,
               match_color='lightgrey',
               mismatch_color='blue',
               mismatch_background_color='lightskyblue')
    pic_file = path.join(pic_dir, "%s.pdf" % gene_id)
    plt.savefig(pic_file, bbox_inches='tight')
    plt.close('all')
    Msg.info("\tFinished")

def draw_cleanup_alleles(mat_dir, pic_dir, thread):
    pool = Pool(processes=thread)
    Msg.info("Visualizing cleanup allele genes")

    for fn in listdir(mat_dir):
        mat_file = path.join(mat_dir, fn)
        break

    res = []
    allele_io = MatrixIO()
    allele_db = allele_io.read_merge_mat(mat_file)
    for gene_id in allele_db:
        res.append([gene_id, pool.apply_async(__draw_cleanup_alleles, (gene_id, allele_db[gene_id], pic_dir,))])
    pool.close()
    pool.join()

    for asc_fn, r in res:
        try:
            r.get()
        except Exception as e:
            Msg.warn("\tException caught with {}: {}". format(asc_fn, e))

    Msg.info("Visualization finished")
