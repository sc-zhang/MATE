## Introduction
This software is a tool for identifying variants associate with phenotypes from pan-genome.

## Dependencies
### Software
 - Python >= 3.7
 - GMAP >= 2021
 - MAFFT
### Python Modules
 - pathos==0.3.0
 - outlier_utils
 - scipy
 - bioplotz
 - pysam

## Installation
```bash
cd /path/to/install
git clone https://github.com/sc-zhang/MATE.git
pip install -r requirements.txt
chmod +x MATE/mate.py
echo 'export PATH=/path/to/install/MATE:$PATH' >> ~/.bash_profile
source ~/.bash_profile
```

## Usage
```bash                                                                                                                                                                                                                                                                                             ─╯
usage: mate.py [-h] -r REF (-g GENOME | -b BAM) [-l PLOIDY] -p PHENO [-k KMER] [--filter FILTER] -o OUTPUT [-s] [-t THREAD]

options:
  -h, --help            show this help message and exit
  -r REF, --ref REF     Reference cds file
  -g GENOME, --genome GENOME
                        Directory contain all genomes
  -b BAM, --bam BAM     Directory contain all bam files by mapping Reseq reads to reference cds
  -l PLOIDY, --ploidy PLOIDY
                        Ploidy of genomes, only effect with -g default=2
  -p PHENO, --pheno PHENO
                        Directory contain phenotypes for association, if the filename of phenotype starts with "LOW-", means lower value is better
  -k KMER, --kmer KMER  kmer length for cleanup mafft result, default=5
  --filter FILTER       Threshold string, "lower_threshold:upper_threshold:missing_threshold:min_allele", lower_threshold means if one allele with less than this ratio of samples supported, it would be dropped; upper_threshold means if one allele with more than this ratio of samples supported, it would
                        be dropped; missing_threshold means if one gene with more than this ratio of samples marked as absence, it would be dropped; min_allele means if one gene with less than this count of alleles (ignore absence), it would be dropped; default=0.05:1:0.25:1
  -o OUTPUT, --output OUTPUT
                        Output directory
  -s, --show            The multi-alignment of variants would be stored as pdf file if this parameter is set
  -t THREAD, --thread THREAD
                        Thread number, default=10
```
**Notice** the id of genes in reference cds file must not contain invalid characters that cannot use in path, like '/', 
'\', '?', et al.  
**Notice** if one phenotype file in pheno directory starts with "LOW-" (all three letters in up case) means the 
lower value of this phenotype is better than higher.  
**Notice** font "Courier New" is required for visualizing, user can copy ttf file of "Courier New" to 
~/.local/share/fonts and use command below to make cache.
```bash
fc-cache -f -v
rm -rf ~/.cache/matplotlib/
```

## Results
- **07.Association**: asc files in 07.Association contain all significant variants. The significant variants were 
identified by divided samples into groups with same variant for each variant site, then, calculated the average value 
of phenotypes in the same group, after that, a t-test was taken between the top 2 groups with the best average 
values, the site with p-value less than 0.05 were saved as significant variant sites.
- **08.VariantMatrix**: for each phenotype, we classified genes into several alleles with cleaned alignments 
(01.CleanupAlleles) and significant variant site(02.SignificantAlleles) in samples, then for each sample, 
we set allele in this sample as '1', otherwise '0'.
- **09.Visualization**: for important variants, we created a subdirectory "01.Variants", and for each phenotype, 
we created a subdirectory for it in "01.Variants", and drawn the variants on genes, and for each gene, a pdf file 
with the same name of gene was saved. The bases with different colors have different means; for different alleles we 
did same thing as variants, but created a subdirectory "02.Alleles".
```bash
lightgrey:    Match
lightskyblue: Mismatch
red:          Key site
```