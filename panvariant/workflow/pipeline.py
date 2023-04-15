from os import path, makedirs, getcwd, listdir, chdir
from panvariant.io.message import Message as Msg
from panvariant.workflow.mapping import mapping
from panvariant.workflow.cds_extract import CDSExtract
from panvariant.workflow.convert import convert_cds_files_for_mafft
from panvariant.workflow.multi_alignment import mafft_alignment
from panvariant.workflow.variant_caller import variant_caller_for_single_file


def get_sample_set(in_dir):
    sample_set = set()
    for fn in listdir(in_dir):
        sample_set.add('.'.join(fn.split('.')[:-1]))
    return sample_set


def pipeline(args):
    ref_cds = path.abspath(args.ref)
    genome_dir = path.abspath(args.genome)
    out_dir = path.abspath(args.output)
    ploidy = args.ploidy
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
        for fn in listdir(gmap_out_gff3_dir):
            if not fn.endswith('.gff3'):
                continue
            if fn.replace('.gff3', '') not in sample_set:
                is_finished = False
                break
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
        for fn in listdir(extract_cds_dir):
            if fn.replace('.cds', '') not in sample_set:
                is_finished = False
                break
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
    if is_finished:
        Msg.info("MAFFT results found, skipping...")
    else:
        mafft_alignment(converted_cds_dir, out_mafft_dir, thread)

    Msg.info("Step5: Variant calling")
    for fn in out_mafft_dir:
        variant_caller_for_single_file(path.join(out_mafft_dir, fn))
        break

    chdir(cur_dir)
