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
-outFile /tmp/out.fa -chr 21 \
-dir /path/to/repo/test -meiOnly True
```
