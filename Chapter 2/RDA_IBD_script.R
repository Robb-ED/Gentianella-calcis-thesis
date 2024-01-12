#########RDA (Redundancy Analysis Script)############
setwd("~/R/R_stuff/Gentianella/Gentianella_project/Whole_Data/50_percent_250miss_3cov_NOREPS/")

#Packages
library(vcfR)
library(adespatial)
library(codep)
library(adegraphics)
library(vegan)
library(psych)
library(adegenet)

#Set the seed
set.seed(2)

#Load in VCF, convert to dataframe
calcis <- read.vcfR("wholedata_50_250_maf0.05_depth3_nomiss_repsremoved.vcf")
genlight <- vcfR2genlight(calcis)

#Pop data 
pop.data <- read.table("popmap_full_altsorting_50percentrem.txt", sep = "\t", header = TRUE)

pop(genlight) <- pop.data$Population

#Convert genlight object to matrix and genind
x.mat <- as.matrix(genlight) # x is a genlight object
x.mat[x.mat == 0] <- "1/1" # homozygote reference
x.mat[x.mat == 1] <- "1/2" # heterozygote
x.mat[x.mat == 2] <- "2/2" # homozygote alternate
gentian_genetics <- df2genind(x.mat, sep = "/", pop = pop(genlight), ploidy = 2)

#Environmental Variables
env_data <- read.table("RDA_info.txt", sep ="\t", header = TRUE)

#Get them in the same order so that they correspond to samples
Popandenv <- merge(env_data, pop.data, by = "Lab_NO")

#isolate the distance metrics
gentian_distances <- data.frame(matrix(ncol = 2, nrow = 92))
colnames(gentian_distances) <- c("long","lat")
gentian_distances$long <- env_data$Long
gentian_distances$lat <- env_data$Lat


#investigate their distribution
plot(gentian_distances)

#Spatial distances among sites 
dist.spatial <- dist(gentian_distances) 

#do the dbmem - might need some kind of model selection, retain ones that are important (most parsimonious)
db.mem <- dbmem(dist.spatial, silent = FALSE)
summary(db.mem)

#Check to see if there needs to be detrending
Dmod <- rda(gentian_genetics ~ env_data$Lat + env_data$Long)
ordistep(Dmod, permutations = how(nperm = 99))
#Looks like there is a significant linear gradient 

dtgenmatrix <- resid(Dmod) #detrended genetic matrix

## Examine the MEMs - use the same model selection process on the detrended data
MEMmod<-rda(dtgenmatrix~., data=db.mem)
ordistep(MEMmod, permutations=how(nperm=99))
#Looks like all are significant 

#Predictor isolation
pred <- Popandenv[,c(8,9,10,11,15)]

#Test to see viff
testrda <- rda(gentian_genetics~., data=pred)
vif.cca(testrda)

ordistep(testrda)


##add to dist data
distances.mem <- gentian_distances
distances.mem <- cbind(gentian_distances, db.mem)


#Perform RDA accounting for geography
dist.rda <- rda(gentian_genetics~.+Condition(as.matrix(distances.mem)), data=pred)
summary(dist.rda)


#Summary stats etc
anova.cca(dist.rda, by = "margin", permutations=how(nperm=9999))
RsquareAdj(dist.rda)
summary(eigenvals(dist.rda, model = "constrained"))
screeplot(dist.rda)
vif.cca(dist.rda) #VIF scores will all be really high because distances have been included


#Make a pretty plot
pop <- Popandenv$Population
calcis_colours <- c("#EAB64D","#8e8f8c","#0BB24E","#5c6f9d","deepskyblue","#DA8051","#857E61"
                    ,"#90C8AD","orchid4","yellowgreen","#C97E8C","#CBD4C2")

poplist <- unique(Popandenv$Population)

# axes 1 & 2
plot(dist.rda, type="n", scaling=3)
points(dist.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(dist.rda, display="sites", pch=21, cex=1.3, scaling=3, bg = calcis_colours[as.factor(Popandenv$Population)]) # the Gentians
text(dist.rda, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("topleft", legend=poplist, bty="o", col="gray32", pch=21, bg="lightgrey",
pt.cex = 1, cex=0.75, title = "Populations", title.cex = 1, title.adj = 0.1,
pt.bg=calcis_colours[as.factor(unique(Popandenv$Population))]) 


# axes 3 & 4
plot(dist.rda, type="n", scaling=3, choices=c(3,4))
points(dist.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(dist.rda, display="sites", pch=21, cex=1.3, scaling=3, bg = calcis_colours[as.factor(Popandenv$Population)]) # the Gentians
text(dist.rda, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("topleft", legend=poplist, bty="o", col="gray32", pch=21, bg="lightgrey",
       pt.cex = 1, cex=0.75, title = "Populations", title.cex = 1, title.adj = 0.1,
       pt.bg=calcis_colours[as.factor(unique(Popandenv$Population))])


## see what variation is shared
geography <- distances.mem
stone_chemistry <- env_data[,8:9]
enviro_variables <- env_data[,c(10,11,15)]

#parition variance
gentian_var <- varpart(gentian_genetics, geography, stone_chemistry, enviro_variables)
gentian_var

#Venn
plot(gentian_var, digits = 2, Xnames = c('Neutral Genetic Structure', 'Limestone', "Climate"), 
     bg=c("orchid","limegreen","deepskyblue"))


####Do IBD (Isolation by distance)#######
library(dartR)

#Load in VCF of full dataset, convert to genlight
coords <- data.frame(matrix(ncol = 2, nrow = 92))
colnames(coords) <- c("lat","lon")
coords$lat <- env_data$Lat
coords$lon <- env_data$Long

ibd_result <- gl.ibd(genlight,
                     distance = 'euclidean',
                     coordinates = coords,
                     permutations = 99999)

ibd_result
