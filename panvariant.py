#!/usr/bin/env python3
import argparse
from panvariant.workflow.pipeline import pipeline


def get_opts():
    groups = argparse.ArgumentParser()

    groups.add_argument('-r', '--ref', help='Reference cds file', required=True)
    groups.add_argument('-g', '--genome', help="Directory contain all genomes", required=True)
    groups.add_argument('-l', '--ploidy', help="Ploidy of genomes, default=2", type=int, default=2)
    groups.add_argument('-p', '--pheno', help="Directory contain phenotypes for association", required=True)
    groups.add_argument('-o', '--output', help="Output directory", required=True)
    groups.add_argument('-t', '--thread', help="Thread number, default=10", type=int, default=10)

    return groups.parse_args()


if __name__ == "__main__":
    opts = get_opts()
    pipeline(opts)
