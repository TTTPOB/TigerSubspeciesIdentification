# Tiger Subspecies Identifier
This is a snakemake pipepline designed for tiger subspecies identifying. It can speed up the process by using subset of the whole genome.

## Input files
This workflow requires；
 - A file descibes which scaffolds should be used (txt)
 - The orignal, high quality snp set (vcf.gz)
 - fastq.gz files to be identified
 - The whole genome