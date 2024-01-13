# Chapter 3: Stacks of SNPs at the disco: Comparing two SNP discovery approaches using _G. calcis_ genetic data

Here are the data processing methods I used for SNP discovery in "South Canterbury" _G. calcis_ samples using Stacks. Demultiplexing, adapter removal and cutside trimming was performed by the Elshire Group using their standard approaches.

## 1: Exclusion of short reads
After assessing quality of reads using [fastQC](https://github.com/s-andrews/FastQC) and [MultiQC](https://github.com/MultiQC/MultiQC) I used TrimGalore 0.6.7 to remove reads of less than 100BP. This enabled better SNP discovery in Stacks by removing short sequences that might be difficult to form into loci. This also had the advantage of removing a small amount of adapter contamination that was still present.

## 2: SNP discovery in Stacks
Stacks is [well-known](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12775) for needing extensive parameter optimisation prior to final SNP discovery. This script was designed to trial the combinations M/n 1-9 for each value of 3, 6 and 9 for m using a handful of representative samples from South Canterbury populations. 
