##################Analyses without prior assignment###################
#In this script we will perform all necessary population genetic analyses that don't need prior assignment
#Start by setting working directory, packages and reading in the data

setwd("~/R/R_stuff/Gentianella/Gentianella_project/Southcant_subset_discoSNP/")

#Libraries
library(ggplot2)
library(vcfR)
library(adegenet)
library(RColorBrewer)
library(pals)
library(ggpubr)
library(reshape2)

#Choose the colour palette for the dataset you are working with
calcis_colours <- c("#0BB24E","#5c6f9d","#857E61","#90C8AD","orchid4","yellowgreen")

calcis_colours <- c("#C97E8C","deepskyblue","#EAB64D","orchid4","yellowgreen","#DA8051")

calcis_colours <- c("#EAB64D","#AFC494","#0BB24E","#5c6f9d","deepskyblue","#DA8051","#857E61"
                    ,"#90C8AD","orchid4","yellowgreen","#C97E8C","#CBD4C2")

#read in file
calcis <- read.vcfR("disco_southcant_sspc_maf0.05_cov_6_miss20_noreps_noGHTA1_4.vcf")

#read in pop data - this will allow us to see if groups match up with a prior expectations
pop.data <- read.table("popmap_southcant_noreps_nobadGHTA.txt", sep = "\t", header = T) #header to F if below is true

#colnames(pop.data) <- (c("AccessID","Population")) #if popfile doesn't have this

#Create Genlight and add information to it
v_genlight <- vcfR2genlight(calcis)
pop(v_genlight) <- pop.data$Population
indNames(v_genlight) <- pop.data$AccessID


##############Analyses############
#PCA using adgenet and ggplot2
v_PCA <- glPca(v_genlight, center = TRUE, scale = TRUE, nf = 5, loadings = TRUE,
               alleleAsUnit = FALSE, useC = TRUE, parallel = FALSE,
               n.cores = NULL, returnDotProd=FALSE, matDotProd=NULL)

barplot(100*v_PCA$eig/sum(v_PCA$eig), col = heat.colors(90), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)

#Calculate variation explained by axes
eig <- round((100*v_PCA$eig/sum(v_PCA$eig)), digits = 1)

#Extract PCA scores and merge with the population data
pca_scores <- as.data.frame(v_PCA$scores)
pca_scores["AccessID"] <- rownames(pca_scores)
row.names(pca_scores) <- NULL
pca_scores <- merge(pca_scores, pop.data, by ="AccessID")

#PLOT the PCA - the sprintf function allows us to use a variable in a string 
ggplot(pca_scores, aes(x=PC1, y = PC2)) + #change axes here if needed
  geom_point(aes(col = Population), size =2) + #change symbol here if needed
  stat_ellipse(aes(group = Population, color = Population)) + #change groupings if needed
  xlab(sprintf("PC1 %s %%", eig[1])) + #Change axes
  ylab(sprintf("PC2 %s %%", eig[2]))+ #Change axes
  ggtitle("PCA Stacks Dataset")+ #title
  scale_color_manual(values = calcis_colours)

####DAPC####
#DAPC using groups found via K-means clustering
#https://grunwaldlab.github.io/Population_Genetics_in_R/gbs_analysis.html based off this tutorial
set.seed(258)

maxK <- 10 #set the max value you want
myMat <- matrix(nrow=10, ncol=maxK)
colnames(myMat) <- 1:ncol(myMat)
for(i in 1:nrow(myMat)){
  grp <- find.clusters(v_genlight, n.pca = 40, choose.n.clust = FALSE,  max.n.clust = maxK)
  myMat[i,] <- grp$Kstat
}

my_df <- melt(myMat)
colnames(my_df)[1:3] <- c("Group", "K", "BIC")
my_df$K <- as.factor(my_df$K)
head(my_df)

p1 <- ggplot(my_df, aes(x = K, y = BIC))
p1 <- p1 + geom_boxplot()
p1 <- p1 + theme_bw()
p1 <- p1 + xlab("Number of groups (K)")
p1

#This sets the range of K to run
my_k <- 2:10 #change this depending on what you want the top right panel to be

grp_l <- vector(mode = "list", length = length(my_k))
dapc_l <- vector(mode = "list", length = length(my_k))

for(i in 1:length(dapc_l)){
  grp_l[[i]] <- find.clusters(v_genlight, n.pca = 40, n.clust = my_k[i])
  dapc_l[[i]] <- dapc(v_genlight, pop = grp_l[[i]]$grp, n.pca = 40, n.da = my_k[i])
}


my_df <- as.data.frame(dapc_l[[ length(dapc_l) ]]$ind.coord)
my_df$Group <- dapc_l[[ length(dapc_l) ]]$grp
head(my_df)

#Use either RcolourBrewer or Pals to get colours 
#my_pal <- RColorBrewer::brewer.pal(n = 10, name = "Set3")
my_pal <- cols25(n=10)

p2 <- ggplot(my_df, aes(x = LD1, y = LD2, color = Group, fill = Group))
p2 <- p2 + geom_point(size = 4, shape = 21)
p2 <- p2 + theme_bw()
p2 <- p2 + scale_color_manual(values=c(my_pal))
p2 <- p2 + scale_fill_manual(values=c(paste(my_pal, "66", sep = "")))
p2


tmp <- as.data.frame(dapc_l[[1]]$posterior)
tmp$K <- my_k[1]
tmp$Population <- rownames(tmp)
tmp <- melt(tmp, id = c("Population", "K"))
names(tmp)[3:4] <- c("Group", "Posterior")
tmp$Region <- pop.data$Population #change if needed
my_df <- tmp

for(i in 2:length(dapc_l)){
  tmp <- as.data.frame(dapc_l[[i]]$posterior)
  tmp$K <- my_k[i]
  tmp$Population <- rownames(tmp)
  tmp <- melt(tmp, id = c("Population", "K"))
  names(tmp)[3:4] <- c("Group", "Posterior")
  tmp$Region <- pop.data$Population #change if needed
  
  my_df <- rbind(my_df, tmp)
}

grp.labs <- paste("K =", my_k)
names(grp.labs) <- my_k

p3 <- ggplot(my_df, aes(x = Population, y = Posterior, fill = Group))
p3 <- p3 + geom_bar(stat = "identity")
p3 <- p3 + facet_grid(K ~ Region, scales = "free_x", space = "free", 
                      labeller = labeller(K = grp.labs))
p3 <- p3 + theme_bw()
p3 <- p3 + ylab("Posterior membership probability")
p3 <- p3 + theme(legend.position='none')
p3 <- p3 + scale_fill_manual(values=c(my_pal))
p3 <- p3 + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p3 <- p3 + theme(axis.text.y = element_text(size = 4, colour = "black"))
p3

ggarrange(ggarrange(p1,
                    p2,
                    ncol = 2, labels = c("A", "B")),
          p3,
          nrow = 2,
          labels = c("", "C"),
          heights = c(1, 2)
)


