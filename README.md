# Hiding-in-plain-sight---unannotated-structural-variants-in-public-genomic-data-sets

## Getting Started
1. [Download miniconda](https://docs.conda.io/en/latest/miniconda.html) and choose Python 3.7 64-bit installer.
```
Make sure to open a new terminal and test with `which conda`.
```
For example:
```
MDMBASTYLIANOU:~ astylianou900045$ which conda
/Users/astylianou900045/miniconda3/bin/conda
```
2. Add conda channels - do this to get all the packages needed for this project:
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```
Then check that everything installed correctly with: conda config --show channels

For example:
```
MDMBASTYLIANOU:~ astylianou900045$ conda config --show channels
channels:
  - conda-forge
  - bioconda
  - defaults
```
Now that conda is installed, clone the repo

3. `git clone https://github.com/NCBI-Hackathons/Hiding-in-plain-sight-unannotated-structural-variants-in-public-genomic-data-sets.git`

Cloning will make a new directory. After it is cloned, enter that directory with 'cd' 

For example:
4. `cd` into `Hiding-in-plain-sight-unannotated-structural-variants-in-public-genomic-data-sets` repo

Time to activate the environment we created for this project:

5. ```conda env create -f environment.yml```

6. `conda activate sv_env`

You will know that it was activated because sv_env will show in your command prompt, for example:
```
(sv_env) MDMBASTYLIANOU:~ astylianou900045$
```
7. Run test data in MEI-only mode. Should create a bed file with MEI calls.
```
python /path/to/repo/vcfToBed.py \
-inFile /path/to/repo/test/test.vcf \
-outFile out.fa -chr 21 \
-dir /path/to/repo/test -meiOnly True
```
## Running on other SV types (deletions, inversions, and tandem duplications)
1. [Download human reference sequences for each chromosome for hg19](ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/README.txt)

2. Make BLAST database for each chromosome
```
for f in *.fa.gz; do gunzip -k $f; done
for f in *.fa; do makeblastdb -in $f -dbtype nucl; done
```
3. Make sure your MEI BLAST database and chromosome BLAST databases are in the same directory

4. Run test data
```
python /path/to/repo/vcfToBed.py \
-inFile /path/to/repo/test/test.vcf \
-outFile out.fa -chr 21 \
-dir /path/to/blast_databases/test
```
