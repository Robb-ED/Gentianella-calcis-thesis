# This is a future readme file for chapter 2
Here are the methods I used to explore patterns of population genetic diversity and connectivity in all sampled _Gentianella calcis_ and _G. astonii_ populations. 

## 1 & 2: Preparation for SNP filtering
[vcftools](https://vcftools.github.io/index.html) 0.1.13 was used to gather information on per-sample SNP missingness, minor allele frequency, etc, using modified code from the amazing [Jana Wold](https://github.com/janawold1). These outputs were visualised in R and used to inform SNP filtering performed using the scripts provided by the authors of DiscoSNP-RAD. 

## 3: Analyses without prior assignment
Once filtered, initial analyses that didn't need prior population assignment (e.g., PCA, DAPC) were performed using R.

## 4: Other analyses
Various approaches were then used to explore patterns of genetic diversity and population connectivity in the SNP data. 
