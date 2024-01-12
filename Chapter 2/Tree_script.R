################Create Intermediate files########################
#Use this script for getting input files for things like treemix etc
#set working directory
library(dartR)
library(vcfR)
library(StAMPP)
library(adegenet)

#Import file
calcis <- read.vcfR("ecology_dp8_sspc_paralog_x0.75_nobadSites_maxdp60_TCR.vcf")

#Make sure that your 'population labels' are in the ecology_population_data.txt file. 
#You need three columns with the headers 'AccessID', 'Country', 'Site' 
pop.data <- read.table("ecology_population_data.txt", sep = "\t", header = TRUE)

v_genlight <- vcfR2genlight(calcis)

pop(genlight) <- pop.data$Population


#Write a treemix file
gl2treemix(noreps,
           outfile = "treemix_input.gz",
           outpath = getwd())


#Input for splitstree
#Create a matrix of pairwise Nei's distances between samples
gen_dist <- stamppNeisD(v_genlight, pop = FALSE)

#Write this out as a .phy file to the working directory
stamppPhylip(gen_dist, file = "input_for_splitstree")

#UPGMA
calcis_colours <- c("#C97E8C","deepskyblue","#EAB64D","orchid4","yellowgreen","#DA8051")
library(poppr)
library(ape)

#Create UPGMA tree
tree <- aboot(v_genlight, tree = "upgma", distance = bitwise.dist, sample = 100, 
              showtree = F, cutoff = 50, quiet = T)

#Show UPGMA tree with populations labelled: make sure to list the population names alphabetical below in the 'legend' line!
cols <- brewer.pal(n = nPop(v_genlight), name = "Set3")

plot.phylo(tree, cex = 0.5, font = 2, adj = 0, tip.color =  calcis_colours[pop(v_genlight)])
nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 3, xpd = TRUE)
legend("topright", inset = 0, legend = c("arduana","astonii","calcis","manahune","taiko","waipara"), 
       fill = calcis_colours, border = FALSE, bty = "n", cex = 0.5, yjust = 2)
axis(side = 1)
title(xlab = "Genetic distance")


