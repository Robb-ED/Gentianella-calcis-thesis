##My attempt at Amova, using https://grunwaldlab.github.io/Population_Genetics_in_R/AMOVA.html
#Start by setting working directory, packages and reading in the data
setwd("~/R/R_stuff/Gentianella/Gentianella_project/Whole_Data/50_20miss_norep/")

#Libraries
library(vcfR)
library(ade4)
library(poppr)


##Read in the VCF
calcis <- read.vcfR("wholedata_50rem_0.05maf_20miss_6depth_noreps.vcf")

#read in popdata (you should have four columns if you want to do it the way I did)
pop.data <- read.table("popmap_full_altsorting_50percentrem.txt", sep = "\t", header = TRUE)

v_genlight <- vcfR2genlight(calcis)
pop(v_genlight) <- pop.data$Taxon

#Add the AMOVA strata file
strata(v_genlight) <- (pop.data[,2-4])
nameStrata(v_genlight) <- ~ID/Species/Subsp/Population

#Perform an AMOVA that tests at the population level
amova_obj <- poppr.amova(v_genlight, 
            hier = ~Population,
            threads = 3,
            method = "ade4",
            within = FALSE)

#Check your results
amova_obj
amova_obj$varcomp/sum(amova_obj$varcomp)

#Test if pops significantly different using randomization test
set.seed(2847)
v_sig <- randtest(amova_obj, nrepet = 2000)


#Visualise the results - black line shouldn't fall within observed 
plot(v_sig)
v_sig

#Test to see what happens when we randomize the population structure
v_shuffled <- v_genlight
strata(v_shuffled)

#Set seed for reproducbility 
set.seed(6451)

#Shuffle the population assignments
strata(v_shuffled)[sample(nInd(v_shuffled)),3:4]
strata(v_shuffled) <-strata(v_genlight)[sample(nInd(v_genlight)),3:4]

#amova
V_new <- poppr.amova(v_shuffled, hier = ~Population, method="ade4", within = FALSE)
V_new

#Test
V_new_test <- randtest(V_new, nrepet = 2000)
V_new_test

#See what the result is
plot(V_new_test)

