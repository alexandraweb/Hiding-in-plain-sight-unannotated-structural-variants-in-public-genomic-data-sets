# Hiding-in-plain-sight---unannotated-structural-variants-in-public-genomic-data-sets
1. For front end development -  download R and Rstudio - more info [here](https://www.ics.uci.edu/~sternh/courses/210/InstallingRandRStudio.pdf)

Once Rstudio is installed, import:
```
shiny
dplyr
shinyWidgets
devtools
```
To host the shiny app we will use:
```
https://www.shinyapps.io
```

2. For back end development - follow the steps below to create your environment and install an editor.
I prefer the free version of Pycharm, found [here](https://www.jetbrains.com/pycharm/)

## Getting Started
1. [Download miniconda](https://docs.conda.io/en/latest/miniconda.html) and choose Python 3.7 64-bit installer.
```
After downloading it make sure to open a new terminal and test with `which conda`.
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
```
3. `git clone https://github.com/NCBI-Hackathons/Hiding-in-plain-sight-unannotated-structural-variants-in-public-genomic-data-sets.git`
```
Cloning will make a new directory. After it is cloned, enter that directory with 'cd' 

For example:
```
4. `cd Hiding-in-plain-sight-unannotated-structural-variants-in-public-genomic-data-sets`
```
Time to create and activate the environment for this project:

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
