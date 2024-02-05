# _Genotyping-By-Sequencing Reveals High Population Structure in a Group of Threatened Plants from New Zealand Limestone Ecosystems_
Here are the methods I used to explore patterns of population genetic diversity and connectivity in all sampled _Gentianella calcis_ and _G. astonii_ populations. 

## 1 & 2: Preparation for SNP filtering
[vcftools](https://vcftools.github.io/index.html) 0.1.13 was used to gather information on per-sample SNP missingness, minor allele frequency, etc, using modified code from the amazing [Jana Wold](https://github.com/janawold1). These outputs were visualised in R and used to inform SNP filtering performed using the scripts provided by the authors of DiscoSNP-RAD. 

## 3: Analyses without prior assignment
Once filtered, initial analyses that didn't need prior population assignment (e.g., PCA, DAPC) were performed using R. SNP error was calculated using a modified version of the code from [Mastretta-Yanes et al. (2015)](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12291). The main modifications to this were mostly to enable replicate detection using an alternate naming method, removal of allele error calculation, and little tidy-ups to comments, etc. 

## 4: Other analyses
Various approaches were then used to explore patterns of genetic diversity and population connectivity in the SNP data. 

## 5: *TreeMix* for *OptM*
This script randomly subsamples data from the *TreeMix* input file created previously to create new input files that contain random 80% subsets of the original. These were then used to run [TreeMix](https://bitbucket.org/nygcresearch/treemix/wiki/Home). Once complete, these output files were zipped and used to calculate the optimal number of migration edges in [OptM](https://rfitak.shinyapps.io/OptM/).  

## 6: *TreeMix*
Once the optimal number of migration edges had been calculated, TreeMix was re-run for migration edges 1-10 using 100% of the SNP data. Visualisation was completed in R using the R scripts provided with TreeMix.
