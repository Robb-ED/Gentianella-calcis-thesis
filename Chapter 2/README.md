# This is a future readme file for chapter 2
Here are the methods I used to explore patterns of population genetic diversity and connectivity in all sampled _Gentianella calcis_ and _G. astonii_ populations. 

## 1 & 2: Preparation for SNP filtering
[vcftools](https://vcftools.github.io/index.html) 0.1.13 was used to gather information on per-sample SNP missingness, minor allele frequency, etc, using modified code from the amazing [Jana Wold](https://github.com/janawold1). These outputs were visualised in R and used to inform SNP filtering performed using the scripts provided by the authors of DiscoSNP-RAD. 

## 3: Analyses without prior assignment
Once filtered, initial analyses that didn't need prior population assignment (e.g., PCA, DAPC) were performed using R. SNP error was calculated using a modified version of the code from [Mastretta-Yanes et al. (2015)](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12291). The main modifications to this were mostly to enable replicate detection using an alternate naming method, removal of allele error calculation, and little tidy-ups to comments, etc. 

## 4: Other analyses
Various approaches were then used to explore patterns of genetic diversity and population connectivity in the SNP data. 
