#!/bin/bash -x
set -o errexit


source enviorments/project1_env/bin/activate
#python Part1_RIT.py
#python parseBlastFile.py
#python Part1_RIT_RT.py
#python Part1_RIT.py
#python vcf_to_fastaAS.py -i /mnt/hnas/bioinfo/projects/VariantInfoDownloads/ExAC/gnomAD/Version2.1/genomes/gnomad.genomes.r2.1.sites.vcf.gz
#vcf_to_fastaAS.py
python vcfToFasta.py -inFile /mnt/hnas/bioinfo/projects/VariantInfoDownloads/ExAC/gnomAD/Version2.1/genomes/gnomad.genomes.r2.1.sites.vcf.gz -outFile chr21FastaOOP.fa -chr 21
