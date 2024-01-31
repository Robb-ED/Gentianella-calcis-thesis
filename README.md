# Lonely on Limestone: A Conservation Genomics Study of _Gentianella calcis_

Presented here are used for formal data analysis and processing in Chapters 2 and 3 of my Master's thesis. In each directory I have provided a brief overview of the chapter contents, and details about any tools used (when necessary). To ensure anonymity of landowners and sampling sites, any sample IDs or metadata have been removed from these scripts.

## Chapter 2: _Genotyping-By-Sequencing Reveals High Population Structure in a Group of Threatened Plants from New Zealand Limestone Ecosystems_
SNPs from Genotyping-By-Sequencing data processed by the Elshire Group were used to explore patterns of genetic diversity and population connectivity in _Gentianella calcis_ and _G. astonii_ populations. [DiscoSNP-RAD](https://github.com/GATB/DiscoSnp) was used for SNP discovery.

## Chapter 3: _Stacks of SNPs at the disco: Comparing Two SNP Discovery Approaches Using_ Gentianella calcis _Genetic Data_  
[Stacks 2.64](https://catchenlab.life.illinois.edu/stacks/) was used for discovery SNPs in a subset of the trimmed and demultiplexed reads from the previous chapter. Using this new dataset and a subset of the DiscoSNP-RAD dataset, the two different SNP discovery approaches were compared in terms of data error, RAD loci assembled and population genetic statistics (e.g., genetic diversity, population structure). Read length adjustment and extra trimming used [TrimGalore 0.6.7](https://github.com/FelixKrueger/TrimGalore).  
