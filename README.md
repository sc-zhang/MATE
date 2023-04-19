## Dependencies
### Software
 - GMAP
 - MAFFT
### Python Modules
 - pathos
 - outlier_utils
 - scipy.stats

## Installation
```bash
cd /path/to/install
git clone https://github.com/sc-zhang/PanVariant.git
chmod +x PanVariant/panvariant.py
echo 'export PATH=/path/to/install/PanVariant:$PATH' >> ~/.bash_profile
source ~/.bash_profile
```

## Usage
```bash
usage: panvariant.py [-h] -r REF -g GENOME [-l PLOIDY] -p PHENO -o OUTPUT [-t THREAD]

options:
  -h, --help            show this help message and exit
  -r REF, --ref REF     Reference cds file
  -g GENOME, --genome GENOME
                        Directory contain all genomes
  -l PLOIDY, --ploidy PLOIDY
                        Ploidy of genomes, default=2
  -p PHENO, --pheno PHENO
                        Directory contain phenotypes for association
  -o OUTPUT, --output OUTPUT
                        Output directory
  -t THREAD, --thread THREAD
                        Thread number, default=10
```
**Notice** the id of genes in reference cds file must not contain invalid characters that cannot use in path, like '/', 
'\', '?', et al.

## Results
- **07.Association**: asc files in 07.Association contain all significant variants. The significant variants were 
identified by divided samples into groups with same variant for each variant site, then, calculated the average value 
of phenotypes in the same group, after that, a t-test was taken between the top 2 groups with the highest average 
values, the site with p-value less than 0.05 were saved as significant variant sites.
- **08.Visualization**: for each phenotype, we created a subdirectory for it, and drawn the variants on genes, and for
each gene, a pdf file with the same name of gene was saved.