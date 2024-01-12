####Pairwise Fst Matrix Script######
library(vcfR)
library(adegenet)
library(StAMPP)
library(reshape2)

setwd("~/R/R_stuff/Gentianella/Gentianella_project/Whole_Data/50_20miss_norep")

#Read in files 
#read in vcf file
calcis <- read.vcfR("wholedata_50rem_0.05maf_20miss_6depth_noreps.vcf")

#read in pop data 
pop.data <- read.table("popmap_full_altsorting_50percentrem.txt", sep = "\t", header = T)

#Create genlight and add information to it
v_genlight <- vcfR2genlight(calcis)
pop(v_genlight) <- pop.data$Population
indNames(v_genlight) <- pop.data$AccessID

#Fst with P values from the StAMMP package
pvaluefst <- stamppFst(v_genlight, nboots = 5000, percent = 95.000, nclusters = 3)
fstmat <- pvaluefst$Fsts
##If you use too few bootstraps, your p values will be small (i.e., 2dp)
##To increase this, do more bootstraps and increase the number of cores (clusters)
## to give the calculations more beef

#Make a nice triangle plot 
#Use the matrix of pairwise Fsts to extract their scores
pvaluematrixrounded <- round(pvaluefst$Fsts, digits = 4)
melted_cormat <- melt(pvaluematrixrounded, na.rm = TRUE)

#Heatmap
heatmap <- ggplot(data = melted_cormat, aes(Var2, Var1, label = value, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "yellow", 
                       midpoint = 0, limit = c(0,0.75), space = "Lab", #adjust limits depending on Fst  
                       name="Fst Value") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed() +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_blank())

#add the fst values and prevent them from being in scientific notation
heatmap + 
  geom_text(aes(Var2, Var1, label = format(value, scientific = FALSE)), 
            color = "black", size = 3)

#As the p values are 0 even with large bootstrap values, we need to make them slightly greater in order to apply a BF correction
#Based off trials using lower Fst values, I can be certain everything is at least 0.00001
pvals <- pvaluefst$Pvalues
pvals <- format((pvals[1-12,1-12]+0.00001), scientific = FALSE)

#apply BF correction for multiple comparisons
adjusted <- p.adjust(pvals, method = "bonferroni")
adjusted