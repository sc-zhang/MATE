#!/usr/bin/env python3
import argparse
from mate.workflow.pipeline import pipeline
from mate.__version__ import __version__


def get_opts():
    groups = argparse.ArgumentParser()

    groups.add_argument(
        "-q", "--query", help="Directory contain all query files", required=True
    )
    groups.add_argument("--query_type", help="Type of query file, could be genome/bam/cds, "
                                             "default=genome",
                        choices=["genome", "bam", "cds"], default="genome")
    groups.add_argument("-r", "--reference",
                        help="Reference file of candidate genes",
                        required=True)
    groups.add_argument("--ref_type",
                        help="Type of reference file, could be cds/gff3, when query_type is genome, "
                             "it only could be set to cds, default=cds",
                        choices=["cds", "gff3"], default="cds")
    groups.add_argument(
        "-l",
        "--ploidy",
        help="Ploidy of genomes, only effect with -g, default=1",
        type=int,
        default=1,
    )
    groups.add_argument(
        "-p",
        "--pheno",
        help="Directory contain phenotypes for association, if the "
             'filename of phenotype starts with "LOW-", means lower value is better',
        required=True,
    )
    groups.add_argument("--cds_align",
                        help="If set this parameter, the multi sequences alignment would be applied on cds, "
                             "otherwise it would be converted to pep first",
                        action="store_true")
    groups.add_argument(
        "--variant_filter",
        help="Threshold string for variant caller, "
             '"kmer_length:lower_threshold:missing_threshold", '
             "kmer_length means the size of kmer; "
             "lower_threshold means if one position contain a kmer with less than"
             "this ratio of samples support, drop it; "
             "missing_threshold means for one position if more than this ratio of "
             'samples are "-" at  drop it; '
             "default=5:0.05:0.9",
        default="5:0.05:0.9",
    )
    groups.add_argument(
        "--allele_filter",
        help="Threshold string for final allele construction, "
             '"lower_threshold:upper_threshold:missing_threshold:min_allele", '
             "lower_threshold means if one allele with less than this ratio of "
             "samples supported, it would be dropped; "
             "upper_threshold means if one allele with more than this ratio of "
             "samples supported, it would be dropped; "
             "missing_threshold means if one gene with more than this ratio of "
             "samples marked as absence, it would be dropped; "
             "min_allele means if one gene with less than this count of alleles "
             "(ignore absence), it would be dropped; "
             "default=0:1:1:1 (no filter)",
        default="0:1:1:1",
    )
    groups.add_argument("-o", "--output", help="Output directory", required=True)
    groups.add_argument(
        "-s",
        "--show",
        help="The multi-alignment of variants would be stored if this "
             "parameter is set",
        action="store_true",
    )
    groups.add_argument(
        "--format",
        help="The picture format of multi-alignment, only effect when parameter 'show' is setting, "
             "default=pdf",
        default="pdf",
    )
    groups.add_argument(
        "-t", "--thread", help="Thread number, default=10", type=int, default=10
    )
    groups.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s {version}".format(version=__version__),
    )
    return groups.parse_args()


if __name__ == "__main__":
    opts = get_opts()
    pipeline(opts)
