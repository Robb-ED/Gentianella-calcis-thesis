#####Outlier Detection######
#Following this tutorial https://bcm-uga.github.io/pcadapt/articles/pcadapt.html
#Set working directory, read in data
setwd("~/R/R_stuff/Gentianella/Gentianella_project/Whole_Data/50_percent_250miss_3cov_NOREPS/")
library(vcfR)
library(adegenet)
library(dartR)

#colours
calcis_colours <- c("#EAB64D","#AFC494","#0BB24E","#5c6f9d","deepskyblue","#DA8051","#857E61"
                    ,"#90C8AD","orchid4","yellowgreen","#C97E8C","#CBD4C2")

##Read the files
ecology <- read.vcfR("wholedata_50_250_maf0.05_depth3_nomiss_repsremoved.vcf")

#Population data
pop.data <- read.table("popmap_full_altsorting_50percentrem.txt", sep = "\t", header = TRUE)

v_genlight <- vcfR2genlight(ecology)
pop(v_genlight) <- pop.data$Subspecies


##Using PCAdapt####
library(dartR)
library(pcadapt)

#Set the path to the VCF file
path <- "~/R/R_stuff/Gentianella/Gentianella_project/Whole_Data/50percent_250miss_3cov/wholedata_50percent_less250miss_maf0.05_nomiss_3cov.vcf"

#Read in the file
v_adapt <- read.pcadapt(path, type = ("vcf"))

##Perform the pcadapt - specify a number of axes to retain initially
pc_v <- pcadapt(input=v_adapt,K=20)

#Check the axes using a screeplot. We are looking for where the points deviate from the straight
#Line (i.e., Cattell's rule)
plot(pc_v, option ="screeplot")
#Based on this it looks like K=10 is most appropriate

#Repeat with K=10
pc_v <- pcadapt(input=v_adapt,K=10)
plot(pc_v,option="manhattan")
plot(pc_v, option = "qqplot")

#Plot the p values, those with really small ones are outliers
hist(pc_v$pvalues, xlab = "p-values", main = "Frequency of p-values", breaks = 50, col = "orchid3")

#Add the P values to their own dataframe, rename the column
pval <- as.data.frame(pc_v$pvalues)
colnames(pval)[1] <- "p_value"
#Add the SNP number to the next column (assumes the pvalues kept the order of the SNPS)
pval$SNP <- v_genlight$loc.names
#Perform bonferroni method which is a conservative way of finding SNP outliers
pval$bf <- p.adjust(pval$p_value, method="bonferroni")
#add 1s to identify which SNPs are outliers and which aren't (0)
pval$outlier <- ifelse(pval$bf < 0.1, '1',
                       ifelse(pval$p_value > 0.1, '0', '0'))
#Check how many there are
length(which(pval$outlier==1))

#Separate the Two dataframes
SNP_ouliers <- pval[which(pval$outlier==1),]
SNP_normals <- pval[which(pval$outlier==0),]

#Create a new GL using dartR that contains only the outlier loci
gl_outliers <- gl.keep.loc(v_genlight, SNP_ouliers$SNP)

##Now run it back through the PCA code 
o_PCA <- glPca(gl_outliers, center = TRUE, scale = FALSE, nf = 5, loadings = TRUE,
               alleleAsUnit = FALSE, useC = TRUE, parallel = FALSE,
               n.cores = NULL, returnDotProd=FALSE, matDotProd=NULL)


#Variation explained by axes
eig <- round((100*o_PCA$eig/sum(o_PCA$eig)), digits = 1)

#Extract PCA scores and merge with the population data
pca_scores <- as.data.frame(o_PCA$scores)
pca_scores["Lab_NO"] <- rownames(pca_scores)
row.names(pca_scores) <- NULL
pca_scores <- merge(pca_scores, pop.data, by ="Lab_NO")

#PLOT - the sprintf function allows us to use a variable in a string 
outlierPCA <- ggplot(pca_scores, aes(x=PC1, y = PC2)) + 
  geom_point(aes(color=Population, shape = Taxon), size =2) +
  stat_ellipse(aes(group = Population, color = Population)) +
  xlab(sprintf("PC1 %s %%", eig[1])) +
  ylab(sprintf("PC2 %s %%", eig[2])) +
  ggtitle("PCA of outliers") +
  scale_colour_manual(values=calcis_colours)

outlierPCA

##Repeat this but with neutral loci
gl_normals <- gl.keep.loc(v_genlight, SNP_normals$SNP)
##Now run it back through the PCA code 
n_PCA <- glPca(gl_normals, center = TRUE, scale = FALSE, nf = 5, loadings = TRUE,
               alleleAsUnit = FALSE, useC = TRUE, parallel = FALSE,
               n.cores = NULL, returnDotProd=FALSE, matDotProd=NULL)
 

#Variation explained by axes
eig <- round((100*n_PCA$eig/sum(n_PCA$eig)), digits = 1)

#Extract PCA scores and merge with the population data
pca_scores <- as.data.frame(n_PCA$scores)
pca_scores["Lab_NO"] <- rownames(pca_scores)
row.names(pca_scores) <- NULL
pca_scores <- merge(pca_scores, pop.data, by ="Lab_NO")

#PLOT - the sprintf function allows us to use a variable in a string 
normalPCA <- ggplot(pca_scores, aes(x=PC1, y = PC2)) + 
  geom_point(aes(color=Population, shape=Taxon), size =2) +
  stat_ellipse(aes(group = Population, color = Population)) +
  xlab(sprintf("PC1 %s %%", eig[1])) +
  ylab(sprintf("PC2 %s %%", eig[2])) +
  ggtitle("PCA of neutral SNPS") +
  scale_colour_manual(values = calcis_colours)

normalPCA

#Write these outliers to their own VCF
gl_normals$other$loc.metrics <- as.data.frame(gl_normals$other$loc.metrics)
gl_normals$loc.all <- rep("G/C",nLoc(gl_normals))


gl2vcf(
  gl_normals,
  plink_path = getwd(),
  outfile = "normal_vcf",
  outpath = getwd(),
  snp_pos = "0",
  snp_chr = "0",
  chr_format = "character"
)

#Now do the outliers
gl_outliers$other$loc.metrics <- as.data.frame(gl_outliers$other$loc.metrics)
gl_outliers$loc.all <- rep("G/C",nLoc(gl_outliers))

gl2vcf(
  gl_outliers,
  plink_path = getwd(),
  outfile = "outlier_vcf",
  outpath = getwd(),
  snp_pos = "0",
  snp_chr = "0",
  chr_format = "character"
)


####Fst outlier#####
#Using these tutorials https://baylab.github.io/MarineGenomics/week-7-fst-and-outlier-analysis.html
#http://rstudio-pubs-static.s3.amazonaws.com/305384_9aee1c1046394fb9bd8e449453d72847.html

#Note that for the thesis I didn't end up using this, but I've included the code I created to show
#how I could have done it

#What we will need for this
library(dartR)
library(OutFLANK)
library(qvalue)


##Use the gl.outflank wrapper from dartR to change the genlight into a outflank object
outflank <- gl.outflank(v_genlight, Hmin = 0.1, plot = TRUE, k=3, RightTrimFraction = 0.2, 
                        LeftTrimFraction = 0.1)
head(outflank)


outflank$OutlierFlag %>% summary

#plot to see if there are any low het loci with high fst
plot(outflank$outflank$results$He, outflank$outflank$results$FSTNoCorr)

#Check to see if the corrected fst from uncorrected fst deviate the same amount
plot(outflank$outflank$results$FST, outflank$outflank$results$FSTNoCorr, col ="#5c6f9d")
abline(0,1)

#summary
summary(outflank$outflank$results$FST)

#Histogram - check the Fst
hist(outflank$outflank$results$FSTNoCorr, breaks = seq(0,1.0, by=0.005),
     col = "#5c6f9d")


OutFLANKResultsPlotter(outflank$outflank,
                       NoCorr = TRUE, Hmin = 0.0, binwidth = 0.001, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)
##Different type of histogram if you would like it
hist(X,freq=TRUE,
     col="grey")# prob=TRUE for probabilities not counts
lines(density(X,na.rm= TRUE), col="blue", lwd=2) # add a density estimate with defaults

##Right tail
hist(outflank$outflank$results$pvaluesRightTail)

##Manhattan Plots(uncorrected estimate)
my_out <- outflank$outflank$results$OutlierFlag==TRUE
plot(outflank$outflank$results$He, outflank$outflank$results$FST, pch=19, col=rgb(0,0,0,0.1))
points(outflank$outflank$results$He[my_out], outflank$outflank$results$FST[my_out], col="blue")

##Corrected estimate
plot(outflank$outflank$results$LocusName[outflank$outflank$results$He>0.15], 
     outflank$outflank$results$FST[outflank$outflank$results$He>0.15],
     xlab="Position", ylab="FST", col=rgb(0,0,0,0.2))
points(outflank$outflank$results$LocusName[my_out], outflank$outflank$results$FST[my_out], 
       col="deepskyblue", pch=22)

#Check how many loci have been flagged as outliers 
sum(outflank$outflank$results$OutlierFlag==TRUE)
#Write the outlier loci to a new dataframe
top_candidates <- outflank$outflank$results$OutlierFlag==TRUE
topcan <- outflank$outflank$results[top_candidates,]
topcan[order(topcan$LocusName),]

###Compare and contrast the two outlier methods####
#See what loci the two methods differ by
length(setdiff(SNP_ouliers$SNP, topcan$LocusName))
length(setdiff(topcan$LocusName, SNP_ouliers$SNP))
length(intersect(SNP_ouliers$SNP, topcan$LocusName))

#What we need to do now is combine the two datasets to see cumulatively if the loci reveal anything

#Create a list of all outlier loci
newloci <- setdiff(topcan$LocusName, SNP_ouliers$SNP)
newloci <- append(newloci, (setdiff(SNP_ouliers$SNP, topcan$LocusName)))
#Check to make sure that this is the sum of the loci found by pcadapt and those that are unique
length(newloci)


