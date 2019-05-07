# Hiding-in-plain-sight---unannotated-structural-variants-in-public-genomic-data-sets

## Getting Started
1. [Download miniconda](https://docs.conda.io/en/latest/miniconda.html) and install
Make sure to open a new terminal and test with `which conda`

2. Add conda channels:
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

3. `git clone git@github.com:NCBI-Hackathons/Hiding-in-plain-sight-unannotated-structural-variants-in-public-genomic-data-sets.git`

4. `cd` into `Hiding-in-plain-sight-unannotated-structural-variants-in-public-genomic-data-sets` repo

5. ```conda env create -f environment.yml```

6. `conda activate sv_env`

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
