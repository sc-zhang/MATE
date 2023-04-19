import matplotlib as mpl
import matplotlib.pyplot as plt
from panvariant.io.file_operate import FastaIO, AssociateIO
from panvariant.io.message import Message as Msg
from os import listdir, path, makedirs
from pathos.multiprocessing import Pool


mpl.use('Agg')


def __draw_genes_in_association(ref_cds, asc_file, pic_dir):
    if pic_dir and not path.exists(pic_dir):
        makedirs(pic_dir)
    fast_io = FastaIO(ref_cds)
    fast_io.read_fasta()
    ref_cds_len_db = fast_io.seq_len_db
    var_sites_db = AssociateIO.read_asc(asc_file)

    variant_color_db = {'Missing': 'black',
                        '-': 'crimson',
                        'A': 'darkgreen',
                        'T': 'deeppink',
                        'C': 'darkorange',
                        'G': 'gold'}
    for gene in var_sites_db:
        Msg.info("\tPlotting %s in %s" % (gene, pic_dir))
        plt.figure(figsize=(20, 8), dpi=300)
        best_cnt = {}
        for sample_idx in range(len(var_sites_db[gene]['geno'])):
            best_cnt[sample_idx] = 0
            plt.plot([0, ref_cds_len_db[gene]], [sample_idx, sample_idx], color='grey', linewidth=1)
            for geno_idx in range(len(var_sites_db[gene]['geno'][sample_idx])):
                geno = var_sites_db[gene]['geno'][sample_idx][geno_idx]
                pos = var_sites_db[gene]['pos'][geno_idx]
                ref = var_sites_db[gene]['ref'][geno_idx]
                alt = var_sites_db[gene]['alt'][geno_idx]
                best_geno = var_sites_db[gene]['best_geno'][geno_idx]

                if geno == best_geno:
                    best_cnt[sample_idx] += 1

                if geno == '-':
                    plt.scatter(pos, sample_idx, s=50, marker='x', color='black')
                elif geno == '0':
                    for site in range(len(ref)):
                        ref_geno = ref[site]
                        plt.scatter(pos+site, sample_idx, s=50,
                                    marker="p" if ref_geno == '-' else 's',
                                    color=variant_color_db[ref_geno],
                                    zorder=99)
                else:
                    for site in range(len(alt[int(geno)-1])):
                        alt_geno = alt[int(geno)-1][site]
                        plt.scatter(pos + site, sample_idx, s=50,
                                    marker="p" if alt_geno == '-' else 's',
                                    color=variant_color_db[alt_geno],
                                    zorder=99)

        for variant_type in variant_color_db:
            plt.scatter(0, -1, s=50, marker="p" if variant_type == '-' else 's',
                        color=variant_color_db[variant_type],
                        label=variant_type)
        plt.legend(bbox_to_anchor=(1.01, 0.5), loc="center left", fontsize=20, frameon=False)
        plt.ylim(-.5, len(var_sites_db[gene]['geno'])-.5)
        best_idx = 0
        best_value = 0
        for sample_idx in best_cnt:
            if best_cnt[sample_idx] > best_value:
                best_value = best_cnt[sample_idx]
                best_idx = sample_idx
        plt.text(-25, best_idx, "Best Type", va='center', ha='right', fontsize=15)
        plt.title("Alleles of gene %s" % gene, fontsize=25)
        plt.axis('off')
        pic_file = path.join(pic_dir, "%s.pdf" % gene)
        plt.savefig(pic_file, bbox_inches='tight')
        pic_file = path.join(pic_dir, "%s.png" % gene)
        plt.savefig(pic_file, bbox_inches='tight')
        plt.close('all')
    Msg.info("\tFinished")


def draw_variant_sites_in_association(ref_cds, asc_dir, pic_dir, thread):
    pool = Pool(processes=thread)
    Msg.info("Visualizing allele genes with phenotypes")

    res = []
    for fn in listdir(asc_dir):
        asc_file = path.join(asc_dir, fn)
        pic_sub_dir = path.join(pic_dir, fn.replace('.asc', ''))
        res.append([fn, pool.apply_async(__draw_genes_in_association, (ref_cds, asc_file, pic_sub_dir,))])
    pool.close()
    pool.join()

    for asc_fn, r in res:
        try:
            r.get()
        except Exception as e:
            Msg.warn("\tException caught with {}: {}". format(asc_fn, e))

    Msg.info("Visualization finished")
