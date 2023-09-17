from os import path, makedirs, getcwd, listdir, chdir
from mate.io.message import Message as Msg
from mate.stages import *
from mate.visualization.draw_variants import draw_variant_sites_in_association
from mate.visualization.draw_alleles import draw_alleles
from time import time


def get_sample_set(in_dir):
    sample_set = set()
    for fn in listdir(in_dir):
        sample_set.add('.'.join(fn.split('.')[:-1]))
    return sample_set


def pipeline(args):
    start_time = time()
    ref_cds = path.abspath(args.ref)
    genome_dir = path.abspath(args.genome)
    out_dir = path.abspath(args.output)
    ploidy = args.ploidy
    pheno_dir = path.abspath(args.pheno)
    kmer_length = args.kmer
    thread = args.thread

    cur_dir = path.abspath(getcwd())
    if not path.exists(out_dir):
        makedirs(out_dir)
    Msg.info("Entering %s" % out_dir)
    chdir(out_dir)

    sample_set = get_sample_set(genome_dir)
    Msg.info("Step1: GMAP")
    gmap_out_gff3_dir = path.join(getcwd(), "01.GMAP")
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

    Msg.info("Step2: Extracting CDS")
    extract_cds_dir = path.join(getcwd(), "02.CDS")
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

    Msg.info("Step3: Converting cds")
    converted_cds_dir = path.join(getcwd(), "03.Convert")
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

    Msg.info("Step4: MAFFT")
    out_mafft_dir = path.join(getcwd(), "04.MAFFT")
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

    Msg.info("Step5: Variant calling")
    out_aln_dir = path.join(getcwd(), "05.Variants", "01.CleanupAlign")
    out_var_dir = path.join(getcwd(), "05.Variants", "02.Variant")
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

    Msg.info("Step6: Variant classifying")
    out_cla_dir = path.join(getcwd(), "06.ClassifiedVariants")
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

    Msg.info("Step7: Associating with phenotypes")
    out_asc_dir = path.join(getcwd(), "07.Association")
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

    Msg.info("Step8: Generating variant matrix")
    out_merge_dir = path.join(getcwd(), "08.VariantMatrix")
    is_finished = True
    if not path.exists(out_merge_dir):
        makedirs(out_merge_dir)
        is_finished = False
    else:
        if not listdir(out_merge_dir):
            is_finished = False
    if is_finished:
        Msg.info("Variant matrix found, skipping...")
    else:
        merge_variant_matrix(pheno_dir, out_aln_dir, out_merge_dir, thread)

    Msg.info("Step9: Visualizing variants")
    out_vis_var_dir = path.join(getcwd(), "09.Visualization", "01.Variants")
    out_vis_allele_dir = path.join(getcwd(), "09.Visualization", "02.Alleles")
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
    if not path.exists(out_vis_allele_dir):
        makedirs(out_vis_allele_dir)
        is_finished = False
    else:
        if not listdir(out_vis_allele_dir):
            is_finished = False
    if is_finished:
        Msg.info("Visualization of alleles found, skipping")
    else:
        draw_alleles(out_merge_dir, out_vis_allele_dir, thread)

    Msg.info("Return %s" % cur_dir)
    chdir(cur_dir)
    end_time = time()
    Msg.info("All done, cost: %.2f sec." % (end_time-start_time))
