#!/bin/bash

#Note: These outputs were only used to estimate the optimal number of migration edges and not for the final run of TreeMix and visualisation.

for m in $(seq 1 10) 
do
    for i in $(seq 1 10) 
    do
        # Generate random seed
        s=$RANDOM
        # Generate bootstrapped input file with ~80% of the SNP loci
        gunzip -c treemix_input.gz | awk 'BEGIN {srand()} { if (NR==1) {print $0} else if (rand() <= .8) print $0}' | gzip > treemix${i}.${m}.treemix.gz
        # Run treemix on bootstrapped input file
        treemix -i treemix${i}.${m}.treemix.gz -o out${i}.${m} -noRoot -m ${m} -seed ${s} &
    done &
    done
    
