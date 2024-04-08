from os import path, makedirs, getcwd, listdir, chdir
from mate.io.message import Message as Msg
from mate.stages import *
from mate.visualization.draw_variants import draw_variant_sites_in_association
from mate.visualization.draw_alleles import draw_alleles
from time import time


def get_sample_set(in_dir, filetype="genome"):
    sample_set = set()
    for fn in listdir(in_dir):
        if filetype == "bam" and (not fn.endswith(".bam")):
            continue
        sample_set.add('.'.join(fn.split('.')[:-1]))
    return sample_set

def pipeline(args):
    start_time = time()
    ref_cds = path.abspath(args.ref)
    if args.genome:
        genome_dir = path.abspath(args.genome)
    else:
        genome_dir = ""
    if args.bam:
        bam_dir = path.abspath(args.bam)
    else:
        bam_dir = ""
    out_dir = path.abspath(args.output)
    ploidy = args.ploidy
    pheno_dir = path.abspath(args.pheno)
    kmer_length = args.kmer
    thread = args.thread
    save_pdf = args.show

    cur_dir = path.abspath(getcwd())
    if not path.exists(out_dir):
        makedirs(out_dir)
    Msg.info("Entering %s" % out_dir)
    chdir(out_dir)

    cur_stage = 1
    if genome_dir:
        sample_set = get_sample_set(genome_dir)
        Msg.info("Step%d: GMAP" % cur_stage)
        gmap_out_gff3_dir = path.join(getcwd(), "%02d.GMAP" % cur_stage)
        cur_stage += 1
        is_finished = True
        if not path.exists(gmap_out_gff3_dir):
            makedirs(gmap_out_gff3_dir)
            is_finished = False
        else:
            if not listdir(gmap_out_gff3_dir):
                is_finished = False
            gff3_set = set()
            for fn in listdir(gmap_out_gff3_dir):
                if not fn.endswith('.gff3'):
                    continue
                gff3_set.add(fn.replace('.gff3', ''))
            if gff3_set != sample_set:
                is_finished = False
        if is_finished:
            Msg.info("Gmap results are found, skipping...")
        else:
            mapping(ref_cds, genome_dir, gmap_out_gff3_dir, ploidy, thread)

        Msg.info("Step%d: Extracting CDS" % cur_stage)
        extract_cds_dir = path.join(getcwd(), "%02d.CDS" % cur_stage)
        cur_stage += 1
        is_finished = True
        if not path.exists(extract_cds_dir):
            makedirs(extract_cds_dir)
            is_finished = False
        else:
            if not listdir(extract_cds_dir):
                is_finished = False
            cds_set = set()
            for fn in listdir(extract_cds_dir):
                cds_set.add(fn.replace('.cds', ''))
            if cds_set != sample_set:
                is_finished = False
        if is_finished:
            Msg.info("CDS are extracted, skipping...")
        else:
            cds_extractor = CDSExtract(gmap_out_gff3_dir, genome_dir, extract_cds_dir, thread)
            cds_extractor.extract()
    else:
        sample_set = get_sample_set(bam_dir, filetype="bam")
        Msg.info("Step%d: Converting Bam to CDS" % cur_stage)
        extract_cds_dir = path.join(getcwd(), "%02d.CDS" % cur_stage)
        cur_stage += 1
        is_finished = True
        if not path.exists(extract_cds_dir):
            makedirs(extract_cds_dir)
            is_finished = False
        else:
            if not listdir(extract_cds_dir):
                is_finished = False
            cds_set = set()
            for fn in listdir(extract_cds_dir):
                cds_set.add(fn.replace('.cds', ''))
            if cds_set != sample_set:
                is_finished = False
        if is_finished:
            Msg.info("CDS are extracted, skipping...")
        else:
            bam_convertor = BAM2CDS(bam_dir, ref_cds, extract_cds_dir, thread)
            bam_convertor.convert()

    Msg.info("Step%d: Converting cds" % cur_stage)
    converted_cds_dir = path.join(getcwd(), "%02d.Convert" % cur_stage)
    cur_stage += 1
    is_finished = True
    if not path.exists(converted_cds_dir):
        makedirs(converted_cds_dir)
        is_finished = False
    else:
        if not listdir(converted_cds_dir):
            is_finished = False
    if is_finished:
        Msg.info("CDS are converted, skipping...")
    else:
        convert_cds_files_for_mafft(extract_cds_dir, converted_cds_dir)

    Msg.info("Step%d: MAFFT" % cur_stage)
    out_mafft_dir = path.join(getcwd(), "%02d.MAFFT" % cur_stage)
    cur_stage += 1
    is_finished = True
    if not path.exists(out_mafft_dir):
        makedirs(out_mafft_dir)
        is_finished = False
    else:
        if not listdir(out_mafft_dir):
            is_finished = False
    if is_finished:
        Msg.info("MAFFT results found, skipping...")
    else:
        mafft_alignment(converted_cds_dir, out_mafft_dir, thread)

    Msg.info("Step%d: Variant calling" % cur_stage)
    out_aln_dir = path.join(getcwd(), "%02d.Variants" % cur_stage, "01.CleanupAlign")
    out_var_dir = path.join(getcwd(), "%02d.Variants" % cur_stage, "02.Variant")
    cur_stage += 1
    is_finished = True
    if not path.exists(out_aln_dir):
        makedirs(out_aln_dir)
        is_finished = False
    else:
        if not listdir(out_aln_dir):
            is_finished = False

    if not path.exists(out_var_dir):
        makedirs(out_var_dir)
        is_finished = False
    else:
        if not listdir(out_var_dir):
            is_finished = False

    if is_finished:
        Msg.info("Variant results found, skipping...")
    else:
        variant_caller(out_mafft_dir, out_var_dir, out_aln_dir, kmer_length, thread)

    Msg.info("Step%d: Variant classifying" % cur_stage)
    out_cla_dir = path.join(getcwd(), "%02d.ClassifiedVariants" % cur_stage)
    cur_stage += 1
    is_finished = True
    if not path.exists(out_cla_dir):
        makedirs(out_cla_dir)
        is_finished = False
    else:
        if not listdir(out_cla_dir):
            is_finished = False
    if is_finished:
        Msg.info("Variant classified results found, skipping...")
    else:
        variant_classifier(out_var_dir, out_cla_dir, thread)

    Msg.info("Step%d: Associating with phenotypes" % cur_stage)
    out_asc_dir = path.join(getcwd(), "%02d.Association" % cur_stage)
    cur_stage += 1
    is_finished = True
    if not path.exists(out_asc_dir):
        makedirs(out_asc_dir)
        is_finished = False
    else:
        if not listdir(out_asc_dir):
            is_finished = False
    if is_finished:
        Msg.info("Association results found, skipping...")
    else:
        associate_with_pheno(pheno_dir, out_cla_dir, out_asc_dir, thread)

    Msg.info("Step%d: Generating variant matrix" % cur_stage)
    out_merge_cleanup_dir = path.join(getcwd(), "%02d.VariantMatrix" % cur_stage, "01.CleanupAlleles")
    out_merge_sig_dir = path.join(getcwd(), "%02d.VariantMatrix" % cur_stage, "02.SignificantAlleles")
    cur_stage += 1
    is_finished = True
    if not path.exists(out_merge_cleanup_dir):
        makedirs(out_merge_cleanup_dir)
        is_finished = False
    else:
        if not listdir(out_merge_cleanup_dir):
            is_finished = False

    if not path.exists(out_merge_sig_dir):
        makedirs(out_merge_sig_dir)
        is_finished = False
    else:
        if not listdir(out_merge_sig_dir):
            is_finished = False
    if is_finished:
        Msg.info("Variant matrix found, skipping...")
    else:
        merge_variant_matrix(pheno_dir, out_aln_dir, out_asc_dir, out_merge_cleanup_dir, out_merge_sig_dir, thread)

    if save_pdf:
        Msg.info("Step%d: Visualizing variants" % cur_stage)
        out_vis_var_dir = path.join(getcwd(), "%02d.Visualization" % cur_stage, "01.Variants")
        out_vis_cleanup_allele_dir = path.join(getcwd(), "%02d.Visualization" % cur_stage, "02.Alleles", "01.CleanupAlleles")
        out_vis_sig_allele_dir = path.join(getcwd(), "%02d.Visualization" % cur_stage, "02.Alleles", "02.SignificantAlleles")
        cur_stage += 1
        is_finished = True
        if not path.exists(out_vis_var_dir):
            makedirs(out_vis_var_dir)
            is_finished = False
        else:
            if not listdir(out_vis_var_dir):
                is_finished = False
        if is_finished:
            Msg.info("Visualization of variants found, skipping")
        else:
            draw_variant_sites_in_association(out_aln_dir, out_asc_dir, out_vis_var_dir, thread)

        is_finished = True
        if not path.exists(out_vis_cleanup_allele_dir):
            makedirs(out_vis_cleanup_allele_dir)
            is_finished = False
        else:
            if not listdir(out_vis_cleanup_allele_dir):
                is_finished = False
        if is_finished:
            Msg.info("Visualization of cleanup alleles found, skipping")
        else:
            draw_alleles(out_merge_cleanup_dir, out_vis_cleanup_allele_dir, thread)

        is_finished = True
        if not path.exists(out_vis_sig_allele_dir):
            makedirs(out_vis_sig_allele_dir)
            is_finished = False
        else:
            if not listdir(out_vis_sig_allele_dir):
                is_finished = False
        if is_finished:
            Msg.info("Visualization of significant alleles found, skipping")
        else:
            draw_alleles(out_merge_sig_dir, out_vis_sig_allele_dir, thread)

    Msg.info("Return %s" % cur_dir)
    chdir(cur_dir)
    end_time = time()
    Msg.info("All done, cost: %.2f sec." % (end_time-start_time))
