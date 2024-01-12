###########Script for testing filtering and related parameters##############
setwd("~/R/R_stuff/Gentianella/Gentianella_project/Whole_Data/Trials/Filtering/")

library(tidyverse)
library(ggplot2)

###Read in the data####
var_miss <- read_delim("total.lmiss", delim = "\t", 
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
var_freq <- read_delim("total.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
var_ind_miss  <- read_delim("total.imiss", delim = "\t",
                            col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
var_ind_het <- read_delim("het.het", delim = "\t",
                          col_names = c("ind","ho_obs", "ho_exp", "nsites", "f"), skip = 1)
#Add my custom colour
col <- "orchid3"

####Make prettiful plots####
# Variant missingness
missdens <- ggplot(var_miss, aes(fmiss)) + 
  geom_density(colour = col, fill =col, lwd=1, alpha = 0.7) +
  labs(title = "Variant Missingness Distribution") +
  theme_light() +
  ylab("Density")+
  xlab("Proportion of samples missing SNP")
missdens


#INDV missing
indvmiss <- ggplot(var_ind_miss, aes(x=ind, y=nmiss)) + 
  geom_col(fill = col, colour = "black", alpha=0.7) +
  labs(title = "Number Variants Missing Per Sample") +
  xlab("Sample") +
  ylab("Number of Variants Missing") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
indvmiss


# Minor allele frequency
var_freq$maf <- var_freq %>% 
  select(a1, a2) %>%
  apply(1, function(z) min(z))

maf_plot <- ggplot(var_freq, aes(maf)) + 
  geom_density(fill = col, colour = "black", alpha = 0.7) +
  labs(title = "MAF distribution") +
  theme_light()

maf_plot


#Heterozygosity
het_plot<- ggplot(var_ind_het, aes(x=ind, y=ho_obs)) + 
  geom_col(fill = col, colour = "black", alpha = 0.7) +
  labs(title = "Homozygosity per individual") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("Number of Homozygous sites")+
  xlab("Sample")
het_plot


#Calculate proportion of heterozygous loci
var_ind_het$he <- (var_ind_het$nsites - var_ind_het$ho_obs)

#Heterozygosity (not homozygosity)
het_plot2<- ggplot(var_ind_het, aes(x=ind, y=he)) + 
  geom_col(fill = col, colour = "black", alpha = 0.7) +
  labs(title = "Heterozygosity per individual") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("Number of Heterozygous sites")+
  xlab("Sample")
het_plot2


# Plot Variant Mean Depth
#Note: VCFtools counted NA as 0, leading to unrealistically low depth estimates. 
#to get around this I've used another strategy which uses a modified version of a function 
#employed by "SNPfiltR" to do it a way that works better. Because this uses the VCF file (which when unfiltered can be large)
#I recommend only using this once you've filtered your data
#SNPfiltR strategy
library(vcfR)

calcis <- read.vcfR("populations.snps.vcf") #your VCF file that is being filtered

dp.matrix <- extract.gt(calcis, element = "DP", as.numeric = TRUE)
snpdepth <- rowSums(dp.matrix, na.rm = TRUE)/rowSums(is.na(dp.matrix) == FALSE)
snpdepth <- as.data.frame(snpdepth[!snpdepth>200])
colnames(snpdepth)[1]="SNP_Depth"

SnpDepthplot <- ggplot(snpdepth, aes(SNP_Depth)) + 
  geom_density(lwd=1, alpha = 0.7) +
  labs(title = "SNP Depth") +
  theme_light() +
  ylab("Density")+
  xlab("Average Depth of SNP") +
  scale_colour_manual(values=col) +
  scale_fill_manual(values=col) +
  scale_x_continuous(n.breaks=20)
SnpDepthplot



