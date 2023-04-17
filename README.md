## Dependencies
### Software
 - GMAP
 - MAFFT
### Python Modules
 - pathos
 - outlier-utils

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
                        Phenotype values of all genomes
  -o OUTPUT, --output OUTPUT
                        Output directory
  -t THREAD, --thread THREAD
                        Thread number, default=10
```
**Notice** the id of genes in reference cds file must not contain invalid characters that cannot use in path, like '/', 
'\', '?', et al.