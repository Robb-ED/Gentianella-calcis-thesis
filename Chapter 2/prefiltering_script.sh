#!/bin/bash

file=YOURFILENAMEHERE.vcf
printf "\nNOW BEGINNING TO CALCULATE STATS...\n"
for vcf in ${file}
  do
      echo "Calculating stats"
      vcftools --gzvcf ${vcf} --out total --freq2 &
      vcftools --gzvcf ${vcf} --out total --missing-site &
      vcftools --gzvcf ${vcf} --out total --missing-indv &
      vcftools --gzvcf ${vcf} --out total --het
      echo "Now done calculating stats"
done