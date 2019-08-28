This is custom code for CRISPR adaptation analysis in KD403 strain.
The analysis is performed in R, ShortRead package will be required. Run R as administrator and install the package by running code:
source("https://bioconductor.org/biocLite.R")
biocLite("ShortRead")

Make sure that a folder with fastq reads also contains a reference genome sequence "kd403-real.fasta"
