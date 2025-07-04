#!/usr/bin/env python3
import argparse
from mate.workflow.pipeline import pipeline


def get_opts():
    groups = argparse.ArgumentParser()

    mut_input_group = groups.add_mutually_exclusive_group(required=True)
    mut_input_group.add_argument('-g', '--genome', help="Directory contain all genomes")
    mut_input_group.add_argument('-b', '--bam', help="Directory contain all bam files by mapping "
                                                     "Reseq reads to reference cds")
    # groups.add_argument('-g', '--genome', help="Directory contain all genomes", required=True)
    mut_ref_group = groups.add_mutually_exclusive_group(required=True)
    mut_ref_group.add_argument('--cds', help="Reference cds file of candidate genes, can be set with "
                                             "-g/--genome and -b/--bam")
    mut_ref_group.add_argument('--bed', help="Reference bed file of candidate genes, only effect with "
                                             "-b/--bam")

    groups.add_argument('-l', '--ploidy', help="Ploidy of genomes, only effect with -g, default=1",
                        type=int, default=1)
    groups.add_argument('-p', '--pheno', help="Directory contain phenotypes for association, if the "
                                              "filename of phenotype starts with \"LOW-\", means lower value is better"
                        , required=True)
    groups.add_argument('--variant_filter', help="Threshold string for variant caller, "
                                                 "\"kmer_length:lower_threshold:missing_threshold\", "
                                                 "kmer_length means the size of kmer; "
                                                 "lower_threshold means if one position contain a kmer with less than"
                                                 "this ratio of samples support, drop it; "
                                                 "missing_threshold means for one position if more than this ratio of "
                                                 "samples are \"-\" at  drop it; "
                                                 "default=5:0.05:0.9", default="5:0.05:0.9")
    groups.add_argument('--allele_filter', help="Threshold string for final allele construction, "
                                                "\"lower_threshold:upper_threshold:missing_threshold:min_allele\", "
                                                "lower_threshold means if one allele with less than this ratio of "
                                                "samples supported, it would be dropped; "
                                                "upper_threshold means if one allele with more than this ratio of "
                                                "samples supported, it would be dropped; "
                                                "missing_threshold means if one gene with more than this ratio of "
                                                "samples marked as absence, it would be dropped; "
                                                "min_allele means if one gene with less than this count of alleles "
                                                "(ignore absence), it would be dropped; "
                                                "default=0:1:1:1 (no filter)", default="0:1:1:1")
    groups.add_argument('-o', '--output', help="Output directory", required=True)
    groups.add_argument('-s', '--show', help="The multi-alignment of variants would be stored if this "
                                             "parameter is set",
                        action="store_true")
    groups.add_argument('--format',
                        help="The picture format of multi-alignment, only effect when parameter 'show' is setting, "
                             "default=pdf", default="pdf")
    groups.add_argument('-t', '--thread', help="Thread number, default=10", type=int, default=10)

    return groups.parse_args()


if __name__ == "__main__":
    opts = get_opts()
    pipeline(opts)
