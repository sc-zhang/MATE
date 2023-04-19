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
    fasta_io = FastaIO(ref_cds)
    fasta_io.read_fasta()
    var_sites_db = AssociateIO.read_asc(asc_file)

    variant_color_db = {'-': 'crimson',
                        'A': 'darkgreen',
                        'T': 'deeppink',
                        'C': 'darkorange',
                        'G': 'gold'}
    for gene in var_sites_db:
        Msg.info("\tPlotting %s in %s" % (gene, pic_dir))
        plt.figure(figsize=(20, 2*len(var_sites_db[gene]['geno'])), dpi=300)
        best_cnt = {}
        min_pos = fasta_io.seq_len_db[gene]
        max_pos = 0
        for sample_idx in range(len(var_sites_db[gene]['geno'])):
            best_cnt[sample_idx] = 0
            for geno_idx in range(len(var_sites_db[gene]['geno'][sample_idx])):
                geno = var_sites_db[gene]['geno'][sample_idx][geno_idx]
                pos = var_sites_db[gene]['pos'][geno_idx]
                ref = var_sites_db[gene]['ref'][geno_idx]
                alt = var_sites_db[gene]['alt'][geno_idx]
                best_geno = var_sites_db[gene]['best_geno'][geno_idx]

                if pos+len(ref) > max_pos:
                    max_pos = pos+len(ref)
                if pos < min_pos:
                    min_pos = pos

                if geno == best_geno:
                    best_cnt[sample_idx] += 1

                if geno == '-':
                    Msg.warn("\tUnexpected variant found")
                elif geno == '0':
                    for site in range(len(ref)):
                        ref_geno = ref[site]
                        plt.scatter(pos+site, sample_idx, s=100,
                                    marker="X" if ref_geno == '-' else 's',
                                    color=variant_color_db[ref_geno],
                                    alpha=0.75,
                                    zorder=99)
                else:
                    for site in range(len(alt[int(geno)-1])):
                        alt_geno = alt[int(geno)-1][site]
                        plt.scatter(pos + site, sample_idx, s=100,
                                    marker="X" if alt_geno == '-' else 's',
                                    color=variant_color_db[alt_geno],
                                    alpha=0.75,
                                    zorder=99)

        for variant_type in variant_color_db:
            plt.scatter(0, -1, s=100, marker="X" if variant_type == '-' else 's',
                        color=variant_color_db[variant_type],
                        label=variant_type if variant_type != '-' else "Deletion")
        plt.legend(bbox_to_anchor=(1.01, 0.5), loc="center left", fontsize=20, frameon=False)

        left = min_pos - fasta_io.seq_len_db[gene] * .1
        right = max_pos + fasta_io.seq_len_db[gene] * .1

        # avoid plot range out of gene range
        left = max(0, left)
        right = min(right, fasta_io.seq_len_db[gene])
        for sample_idx in range(len(var_sites_db[gene]['geno'])):
            plt.plot([left, right], [sample_idx, sample_idx], color='grey', linewidth=1, zorder=1)

        plt.xlim(left-10 if left == 0 else left, right+10 if right == fasta_io.seq_len_db[gene] else right)
        plt.xticks([])
        plt.xlabel("%s (%d~%d:%d)" % (gene, int(left) + 1, int(right), fasta_io.seq_len_db[gene]), fontsize=20)

        plt.ylim(-.1, len(var_sites_db[gene]['geno'])-.9)
        best_idx = 0
        best_value = 0
        for sample_idx in best_cnt:
            if best_cnt[sample_idx] > best_value:
                best_value = best_cnt[sample_idx]
                best_idx = sample_idx
        plt.yticks([best_idx], ["Best"], fontsize=20)
        ax = plt.gca()
        ax.tick_params("y", length=0)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        pic_file = path.join(pic_dir, "%s.pdf" % gene)
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
