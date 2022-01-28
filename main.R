####################################################
### Title: main.R                                ###
### Description: Analysis for Pantel et al. 2022 ###
### Date: 28-01-2022                             ###
####################################################

### License
# zoop-niche: Data and analysis of niche occupancy and niche shifts in rural and urban environments for freshwater zooplankton.

# Copyright (C) 2022  Jelena Holly Pantel
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

## Package dependencies
library(car)
library(missMDA)
library(vegan)
library(ecospat)
library(geometry)
library(factoextra)
library(fitdistrplus)
library(dismo)
library(rJava)
library(geosphere)
library(rgl)
library(reshape2)
library(maps)
library(mapdata)
library(maptools)  #for shapefiles
library(scales)  #for transparency
library(rgdal)
library(MASS)

## Set working directory
setwd("/Users/jhpantel/Documents/GitHub/zoop-niche/data_files")

# Section 1. Import and transform data --------------------------------------------------

## SPACE
# spa: An 81x2 data frame of the lat-long coordinates of each site. 81 sites x 2 spatial coordinates (LAT, LONG).
spa <- read.csv("geographical.csv")
rownames(spa)<-spa[,1]
spa<-spa[,-1]
spa <- spa[-(82:84),]

## COMMUNITY
# spp: A 319x5 dataframe that lists each occurence instance in the dataset, for all 28 species and for 80 sites (site PL4_GRE has no species present). Each row in the data frame gives X (a unique identifier for the species), sites (the unique identifer for the site), LAT and LONG (the latitude and longitude coordinates for the site), and species (the full name of the species).
spp<- read.csv("spp-10-7-2015.csv", header=T)
spp <- spp[spp$sites != "GENT_B" & spp$sites != "MECH_KT2" & spp$sites != "TP_BLAP1_RIV",]
spp$sites <- factor(spp$sites)

# ZP2: A 81x28 matrix with the presence (1) or absence (0) of 28 Cladoceran zooplankton species, at each of 81 sites. ZP3 is a 81x16 matrix with the presence-absence of the 16 species that appear in 5 or more sites.
ZP.pa <- read.csv("ZP-pa.csv",row.names=1)
ZP2<-t(ZP.pa)
ZP3<-ZP2[,which(colSums(ZP2)>=5)]#only species occurring in more than five sites

## ENVIRONMENT
# env: An 81x28 dataframe with observed values of environmental variables at each site
env <- read.csv("SPEEDY_R_continu2.csv", header=T,row.names=1)
env<-env[-(82:84),-(1:26)]

# env data transform
# env2: an 81x21 matrix of environmental variables
mat <- matrix(NA, nrow = 81, ncol =17 )
lgt.shade<-car::logit(env$shade)
mat[,1] <- as.vector(lgt.shade)
mat[,2] <- as.vector(log(env$sludge+1))
mat[,3] <- as.vector((env$sneller))
mat[,4] <- as.vector(log(env$cond+1))
mat[,5] <- as.vector((env$pH))
mat[,6] <- as.vector((env$Depth))
subset <- as.matrix(log(env[,c(14:23,25)]))
mat[,7:17]<-subset
colnames(mat)<-c("shade","sludge","sneller","cond","pH","max depth","Oxygen","Chl a","SPM","ALK","SULF","DOC","total N","total P","available N","available P","surface area")

macro2<-env[,c(1:3,6)]
nb <- estim_ncpPCA(mat,ncp.max=10)
env2<-imputePCA(data.matrix(mat),ncp=2)
env2<-env2$completeObs
macro2<-imputePCA(data.matrix(macro2),ncp=2)
macro2<-macro2$completeObs
macro3<-car::logit(macro2)
env2<-cbind(env2,macro3)
#env2<-decostand(env2[,],"standardize")

# Urbanization information
# groups: An 81x7 dataframe placing all sites in 3 categories (low, med, high) at each of 7 spatial scales (50, 100, 200, 400, 800, 1600, 3200 m)
GBG<-read.csv("GBG.csv",header=T) # GBG values at different perimeters
rownames(GBG)<-GBG[,1]
GBG<-GBG[,-1]
GBG2<-log(GBG+1)
GBG3<-car::logit(GBG)

BVA<-read.csv("BVA.csv") #BIOLOGICALLY VALUABLE AREA
rownames(BVA)<-BVA[,1]
BVA<-BVA[,-1]
#NEW CATEGORIES BASED ON PERCENTAGE BU AREA AND BVA
#3 categories of urbanization and one extra based on BVA
low50<-(which(GBG [,1]<5 & BVA[,7]>20))
lowmed50<-(which(GBG[,1]<5 & BVA[,7]<20))
med50<-(which(GBG[,1]>5 & GBG[,1]<10))
high50<-(which(GBG[,1]>10))
low100<-(which(GBG[,2]<5 & BVA[,7]>20))
lowmed100<-(which(GBG[,2]<5 & BVA[,7]<20))
med100<-(which(GBG[,2]>5 & GBG[,2]<10))
high100<-(which(GBG[,2]>10))
low200<-(which(GBG[,3]<5 & BVA[,7]>20))
lowmed200<-(which(GBG[,3]<5 & BVA[,7]<20))
med200<-(which(GBG[,3]>5 & GBG[,3]<10))
high200<-(which(GBG[,3]>10))
low400<-(which(GBG[,4]<5 & BVA[,7]>20))
lowmed400<-(which(GBG[,4]<5 & BVA[,7]<20))
med400<-(which(GBG[,4]>5 & GBG[,4]<10))
high400<-(which(GBG[,4]>10))
low800<-(which(GBG[,5]<5 & BVA[,7]>20))
lowmed800<-(which(GBG[,5]<5 & BVA[,7]<20))
med800<-(which(GBG[,5]>5 & GBG[,5]<10))
high800<-(which(GBG[,5]>10))
low1600<-(which(GBG[,6]<5 & BVA[,7]>20))
lowmed1600<-(which(GBG[,6]<5 & BVA[,7]<20))
med1600<-(which(GBG[,6]>5 & GBG[,6]<10))
high1600<-(which(GBG[,6]>10))
low3200<-(which(GBG[,7]<5 & BVA[,7]>20))
lowmed3200<-(which(GBG[,7]<5 & BVA[,7]<20))
med3200<-(which(GBG[,7]>5 & GBG[,7]<10))
high3200<-(which(GBG[,7]>10))

groups<-matrix(nrow=81,ncol=7)
groups[low50,1]<-"low"
groups[c(lowmed50),1]<-"crops"
groups[c(med50,lowmed50),1]<-"med"
groups[high50,1]<-"high"
groups[low100,2]<-"low"
groups[c(med100,lowmed100),2]<-"med"
groups[high100,2]<-"high"
groups[low200,3]<-"low"
groups[c(med200,lowmed200),3]<-"med"
groups[high200,3]<-"high"
groups[low400,4]<-"low"
groups[c(med400,lowmed400),4]<-"med"
groups[high400,4]<-"high"
groups[low800,5]<-"low"
groups[c(med800,lowmed800),5]<-"med"
groups[high800,5]<-"high"
groups[low1600,6]<-"low"
groups[c(med1600,lowmed1600),6]<-"med"
groups[high1600,6]<-"high"
groups[low3200,7]<-"low"
groups[c(med3200,lowmed3200),7]<-"med"
groups[high3200,7]<-"high"

colnames(groups)<-c("cat50","cat100","cat200","cat400","cat800","cat1600","cat3200")
name <- gsub("-","_",rownames(GBG))
rownames(groups)<-name

cat50= factor(gsub('(.).', '\\1', groups[,1]))
# cat100= factor(gsub('(.).', '\\1', groups[,2]))
# cat200= factor(gsub('(.).', '\\1', groups[,3]))
# cat400= factor(gsub('(.).', '\\1', groups[,4]))
# cat800= factor(gsub('(.).', '\\1', groups[,5]))
# cat1600= factor(gsub('(.).', '\\1', groups[,6]))
cat3200= factor(gsub('(.).', '\\1', groups[,7]))

## SPECIES OCCURENCE + ENV
# clim: An 81x23 dataframe that gives the LAT, LONG, and value of 21 environmental variables (cscaled to mean 0, sd 1) for all sites
clim<-cbind(spa,env2)

# occ.sp: A 319x28 dataframe with X (a unique identifier for the species), sites (the unique identifier for the site), LAT and LONG, 21 environmental variables, the urbanization level (low, med, high) at 50 and 3200m, and the name of the species. This is structured to have one line for each occurence of species in the 81-site dataset.
occ.sp_test <- na.exclude(ecospat.sample.envar(dfsp=spp,colspxy=3:4,colspkept=1:4,dfvar=clim,colvarxy=1:2,colvar="all",resolution=100))
occ.sp_test$cat50 <- NA
occ.sp_test$cat3200 <- NA
for (k in 1:dim(groups)[1]){
  id2<-which(rownames(groups)[k]==as.character(occ.sp_test[,2]))
  occ.sp_test$cat50[id2] <- as.character(groups[k,1])
  occ.sp_test$cat3200[id2] <- as.character(groups[k,7])
}
occ.sp_test$cat3200 <- as.factor(occ.sp_test$cat3200)
occ.sp_test$cat50 <- as.factor(occ.sp_test$cat50)

ID<-as.numeric(rownames(occ.sp_test))
occ.sp<-cbind(occ.sp_test,spp[rownames(spp) == ID,5]) #add species names

# sp.list: list of species in the dataset
sp.list<-levels(as.factor((occ.sp[,28])))
sp.nbocc<-c()
for (i in 1:length(sp.list)){sp.nbocc<-c(sp.nbocc,length(which(occ.sp[,28] == sp.list[i])))}

#calculate the nb of occurences per species
sp.list <- sp.list[sp.nbocc>4] # remove species with less than 5 occurences
nb.sp <- length(sp.list) #nb of species

# Calculate the number of occurrences in rural and urban habitats
ru_tab <- array(NA,dim=c(16,4),dimnames=list(sp.list,c("rural50","rural3200","urban50","urban3200")))

for(i in 1:length(sp.list)){
  ru_tab[i,1] <- length(which(occ.sp[,28] == sp.list[i] & occ.sp[,26]=="low"))
  ru_tab[i,2] <- length(which(occ.sp[,28] == sp.list[i] & occ.sp[,27]=="low"))
  ru_tab[i,3] <- length(which(occ.sp[,28] == sp.list[i] & occ.sp[,26]=="high"))
  ru_tab[i,4] <- length(which(occ.sp[,28] == sp.list[i] & occ.sp[,27]=="high"))
}

# Section 2. Data Preparation ------------------------------------------------

# Prep 1. Alpha, beta, gamma diversity --------------------------------
# Multipart function gives average alpha, beta and gamma
## Multiplicative diversity partitioning
div <- multipart(ZP2[-10,],nsimul=1000,scales=1,relative=F,global=T,method="r0")

# Double-checking we have these done properly
# See Jost 2007 for explanations. We use order 1, which is Shannon entropy.
## alpha
al <- mean(exp(diversity(ZP2[-10,])))
## gamma
ga <- exp(diversity(colSums(ZP2[-10,])))
## beta
ga/al

# Prep 2. Niche volume ------------------------------------------------
# PCA-ENVIRONMENT
data<-rbind(occ.sp[,5:25],clim[,3:23])
w<-c(rep(0,nrow(occ.sp)),rep(1,nrow(clim)))
## Specify here to use 2 or 4 PC axes for niche volume
pca.cal <-dudi.pca(data, row.w = w, center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)

# convex hull volume
vol<-matrix(nrow=16,ncol=16)
for (i in 1:16){
  for (j in 1:16){
    row.sp1<-which(occ.sp[,28] == sp.list[i]) # rows in data corresponding to sp1
    row.sp2<-which(occ.sp[,28] == sp.list[j]) # rows in data corresponding to sp2
    name.sp1<-sp.list[i]
    name.sp2<-sp.list[j]
    # predict the scores on the axes
    scores.clim<- pca.cal$li[(nrow(occ.sp)+1):nrow(data),] #scores for global climate
    scores.sp1<- pca.cal$li[row.sp1,] #scores for sp1
    scores.sp2<- pca.cal$li[row.sp2,] #scores for sp2
    scores.sp1_2<-rbind(scores.sp1,scores.sp2)
    dim(scores.sp1_2)
    x<-convhulln(scores.sp1_2, "FA")
    vol[i,j]<-x$vol
  }
}

rownames(vol)<-sp.list
colnames(vol)<-sp.list
diag(vol) #diagonal of vol matrix represents the volume of the species alone (species itself)

vol2 <- vol[upper.tri(vol)]

x <- convhulln(scores.clim, "FA")
x$vol

# rural and urban niche volume
r50_vol <- convhulln(scores.clim[groups[,1] == "low",],"FA")
r3200_vol <- convhulln(scores.clim[groups[,7] == "low",],"FA")

u50_vol <- convhulln(scores.clim[groups[,1] == "high",],"FA")
u3200_vol <- convhulln(scores.clim[groups[,7] == "high",],"FA")

# Prep 3. Body size clustering ---------------------------------
BodySize <- read.csv(file = "zoop_info.csv", header = T,row.names=1)
size_dist <- dist(BodySize$Body_size_mm)
grp <- hclust(size_dist)
plot(grp)
BodySize$Body_size_category[c(6:9,15:16)] <- "Large"
BodySize$Body_size_category[c(1:5,10:14)] <- "Small"

# Prep 4. Niche overlap ---------------------------------
# PCA-ENVIRONMENT
data<-rbind(occ.sp[,5:25],clim[,3:23])
w<-c(rep(0,nrow(occ.sp)),rep(1,nrow(clim)))
pca.cal <-dudi.pca(data, row.w = w, center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)

# PCA overview (Table S1)
res.var <- get_pca_var(pca.cal)
eig.val <- get_eigenvalue(pca.cal)
eig.val
res.var$coord
res.var$contrib

## Species combinations - niche overlap, tests of niche equivalency
sp.combn<-combn(c(1:16),2)

## The niche overlap and equivalency tests take a few hours to run. The p-values given in the manuscript will change (slightly, given randomization) each time you run the analysis. I have saved the values of variables that take a long time to run in an R environment called large_vars.RData. I load that here. The code needed to produce I, D, and test variables is given as I_D_test.R.

## To run niche overlap and equivalency tests, uncomment line below
# source(I_D_test.R)
## This produces the variables test, I, and D. The large_vars.RData file also contains additional large variables: Geomean, Geomean_names, and geo_mean. I indicate below the code that was used to produce these.
## To load large variables test, I, D, Geomean, Geomean_names, geo_mean
load("large_vars.RData")

# Prep 5. Co-occurrence patterns  ---------------------------------

# C scores for all species pairs with more than 5 occurences
# To generate new C scores, use the command below. In the manuscript, we use the values that have been saved as Cscores.txt
# Cscore_res <- ecospat.Cscore(ZP3,nperm=10000,outpath=getwd())
# Used in manuscript
Cscore <- read.table("Cscores.txt",header=T,sep="\t")

# Descriptions of distributions, Supplement S1. Figure S1.
descdist(Cscore$obs.C.score)
Cgam <- fitdist(Cscore$obs.C.score,distr="gamma",method="mme")

descdist(Cscore$SES_Cscore)
SESCnorm <- fitdist(Cscore$SES_Cscore,distr="norm",method="mme")

# Number of observed co-occurences for each species pair
co_tab <- array(NA,dim=c(ncol(sp.combn),3),dimnames=list(NULL,c("sp1","sp2","co")))

for(i in 1:(ncol(sp.combn))){
  co_tab[i,1] <- sum(ZP3[,sp.combn[1,i]])
  co_tab[i,2] <- sum(ZP3[,sp.combn[1,i]])
  co_tab[i,3] <- sum(ZP3[,sp.combn[1,i]] == 1 & ZP3[,sp.combn[2,i]] == 1)
}



# Prep 6. Drivers of co-occurence patterns ---------------------------------

## CscoreSES: the Cscore standardized effect size for all pairs of species (for species that appear in 5 or more sites)
CscoreSES <- Cscore$SES_Cscore

# Dscore: a 120 length vector with the observed Dscore for the 120 species pairs
Dscore <- D

# Geomean: a 120 length vector with the mean probability of occurrence of each species pair
## The variables Geomean, Geomean_names, and geo_mean were created in an external file called Geomean.R (which uses the .csv file spat_interpol.csv, as well as two functions in the R files split_test_train.R and spat_interpolate.R). This code takes a very long time to run, so I included these variables when loading large_vars.RData. Thus:

## Dependencies for Geomean.R
## adehabitatMA, adehabitatHR, raster, SDMTools
## To run code to produce Geomean:
# source("Geomean.R")

# Habmean: a 120 length vector with the mean of a species' site-specific habitat suitability value, for each species pair.
Habmean_sp <- array(NA,dim=dim(ZP3)[2])
env3<-decostand(env2[,],"standardize")
clim2<-cbind(spa,env3)
env_use <- clim2[,3:23]

for(i in 1:dim(ZP3)[2]){
  use <- clim[,c(2,1)]
  use$occ <- ZP3[,i]
  me <- dismo::maxent(env_use,use[,3])
  r <- predict(me, env_use)
  Habmean_sp[i] <- mean(r)
}

Habmean_names <- array(NA,dim=c(ncol(sp.combn),2),dimnames=list(NULL,c("sp1","sp2")))
Habmean <- array(NA, dim=ncol(sp.combn))
for(i in 1:(ncol(sp.combn))){
  Habmean_names[i,1] <- sp.list[sp.combn[1,i]]
  Habmean_names[i,2] <- sp.list[sp.combn[2,i]]
  Habmean[i] <- mean(c(Habmean_sp[which(rownames(geo_mean) == sp.list[sp.combn[1,i]])],Habmean_sp[which(rownames(geo_mean) == sp.list[sp.combn[2,i]])]))
}

## Volsum: a 120 length vector with the convex hull volume of all sites where species pairs occurred, in the space of the first 2 PC axes.
Volsum_names <- array(NA,dim=c(ncol(sp.combn),2),dimnames=list(NULL,c("sp1","sp2")))
Volsum <- array(NA, dim=ncol(sp.combn))
for(i in 1:(ncol(sp.combn))){
  Volsum_names[i,1] <- sp.list[sp.combn[1,i]]
  Volsum_names[i,2] <- sp.list[sp.combn[2,i]]
  Volsum[i] <- vol[which(rownames(vol) == sp.list[sp.combn[1,i]]),which(colnames(vol) == sp.list[sp.combn[2,i]])]
}

# Section 3. Statistical tests ------------------------------------------------

# Test 1. a. Niche equivalency test, b. Niche similarity test --------------------------------

## These values are indicated in Table_S2
# Test 1a. SUBSET OF SIGNIFICANT DIFFERENCES IN NICHE EQUIVALENCE:
# Note these values come from the tests that were run in Prep 4. Niche overlap
# using D metric
eq_pval.D <- matrix(ncol=4,nrow=120)
eq_pval.D <- cbind(sp.list[sp.combn[1,]],sp.list[sp.combn[2,]],test[1,1,2,],test[1,2,2,])
id_greater <- which(eq_pval.D[,3]<0.05)
id_lower <- which(eq_pval.D[,4]<0.05)
cbind(eq_pval.D[id_greater,c(1:2,3)],D[id_greater])
cbind(eq_pval.D[id_lower,c(1:2,4)],D[id_lower])

# Test 1b. SUBSET OF SIGNIFICANT DIFFERENCES IN NICHE SIMILARITY:
# Note these values come from the tests that were run in Prep 4. Niche overlap
# using D metric
sim_pval.D<-matrix(ncol=4,nrow=120)
sim_pval.D<-cbind(sp.list[sp.combn[1,]],sp.list[sp.combn[2,]],test[2,1,2,],test[2,2,2,])
id_greater <- which(sim_pval.D[,3]<0.05)
id_lower <- which(sim_pval.D[,4]<0.05)
cbind(sim_pval.D[id_greater,c(1:2,3)],D[id_greater])
cbind(sim_pval.D[id_lower,c(1:2,4)],D[id_lower])

# Test 2. Discriminant Function Analysis and MANOVA --------------------------------

#To be sure the rows are aligned by species
cbind(sp.list,rownames(BodySize))
BodySize$Body_size_category <- factor(BodySize$Body_size_category)

scores.sp<-matrix(nrow=16,ncol=2)
for(i in 1:length(sp.list)) {
  row.sp1<-which(occ.sp[,28] == sp.list[i]) # rows in data corresponding to sp1
  name.sp1<-sp.list[i]
  # predict the scores on the axes
  scores.sp[i,1]<- mean(pca.cal$li[row.sp1,1]) #mean PC1 score for all sites where species occurs
  scores.sp[i,2]<- mean(pca.cal$li[row.sp1,2]) #mean PC2 score for all sites where species occurs
}

dfa <- lda(BodySize$Body_size_category ~ scores.sp[,1] + scores.sp[,2]) #discriminant analysis on mean pc scores of sites where each species occurs

dfa_results <- predict(dfa) #predictions to which subset: large vs small

man_test <- manova(cbind(scores.sp[,1], scores.sp[,2]) ~ BodySize$Body_size_category)
summary(man_test, "Wilks")

# Test 3. Generalist vs. specialist niche overlap --------------------------------
# Specialists - 4 species with lowest volume (< 9.5). Generalists - 4 species with highest volume (> 30.5)
special <- diag(vol)[diag(vol) < 9.5]
general <- diag(vol)[diag(vol) > 30.5]

# Generate all pairs of specialist species
ss_pairs <- c(which(sp.combn[1,] == 7 & sp.combn[2,] == 10),which(sp.combn[1,] == 7 & sp.combn[2,] == 12),which(sp.combn[1,] == 7 & sp.combn[2,] == 13),which(sp.combn[1,] == 10 & sp.combn[2,] == 12),which(sp.combn[1,] == 10 & sp.combn[2,] == 13),which(sp.combn[1,] == 12 & sp.combn[2,] == 13))
# Check the pairs are correct
Volsum_names[ss_pairs,]

# Generate all pairs of generalist species
gg_pairs <- c(which(sp.combn[1,] == 5 & sp.combn[2,] == 6),which(sp.combn[1,] == 5 & sp.combn[2,] == 9),which(sp.combn[1,] == 5 & sp.combn[2,] == 15),which(sp.combn[1,] == 6 & sp.combn[2,] == 9),which(sp.combn[1,] == 6 & sp.combn[2,] == 15),which(sp.combn[1,] == 9 & sp.combn[2,] == 15))
# Check the pairs are correct
Volsum_names[gg_pairs,]

# Generate all pairs of specialist-generalist species
sg <- expand.grid(c(7,10,12,13),c(5,6,9,15))
sg <- t(apply(sg,1,sort))

sg_pairs <- numeric()
for(i in 1:nrow(sg)){
  sg_pairs <- c(sg_pairs,which(sp.combn[1,] == sg[i,1] & sp.combn[2,] == sg[i,2]))
}
# Check the pairs are correct
Volsum_names[sg_pairs,]

# Assemble D values for ANOVA
d_anova <- array(NA,c(28,4),dimnames=list(NULL,c("sp1","sp2","D","pair")))
d_anova <- as.data.frame(d_anova)

d_anova[1:6,1] <- Volsum_names[ss_pairs,1]
d_anova[1:6,2] <- Volsum_names[ss_pairs,2]
d_anova[1:6,3] <- D[ss_pairs]
d_anova[1:6,4] <- "ss"

d_anova[7:12,1] <- Volsum_names[gg_pairs,1]
d_anova[7:12,2] <- Volsum_names[gg_pairs,2]
d_anova[7:12,3] <- D[gg_pairs]
d_anova[7:12,4] <- "gg"

d_anova[13:28,1] <- Volsum_names[sg_pairs,1]
d_anova[13:28,2] <- Volsum_names[sg_pairs,2]
d_anova[13:28,3] <- D[sg_pairs]
d_anova[13:28,4] <- "sg"

d_anova$pair <- as.factor(d_anova$pair)

# Run the ANOVA
d_aov_res <- aov(d_anova$D ~ d_anova$pair)
# Post-hoc t-tests
d_t_res <- pairwise.t.test(d_anova$D, d_anova$pair, p.adjust.method = "none", pool.sd=F)

# Test 4. Urban vs. rural Niche similarity test --------------------------------

## The variables produced here are rural_urban_50 and rural_urban_3200, which are records of each species, their D score for rural vs. urban habitats, the p-values for tests of equivalency and similarity. This also takes awhile to run, and the results included in the manuscript are loaded in large_vars.RData.
## To run the code that produced these results:
# source("rural_urban_test.R")

# Test 5. Urban vs. rural nestedness --------------------------------
# Nestedness temnperature calculation
res <- nestedtemp(ZP3)
nest <- oecosimu(ZP3,nestedtemp,method="r0",nsimul=999)

# Test 6. Spatial isolation test --------------------------------
## This test originates from Pantel et al. 2011, please cite if using this spatial randomization test.
## Pantel, J.H., T.E. Juenger, and M.A. Leibold. 2011. Environmental gradients structure Daphnia pulex × pulicaria clonal distribution. Journal of Evolutionary Biology 24(4): 723-732.

## Calculate spatial distance among sites
space_euc <- vegdist(spa, "euclidean")
space_dist <- as.matrix(space_euc)
spa2<-spa[,c(2,1)]
space_hav<-distm (spa2, fun = distHaversine)
space_dist_hav<-as.matrix(space_hav)

## standard error function
std.error <- function(x){
  st_err <- sd(x) / sqrt(length(x))
  return(st_err)
}

## Calculate observed average distance between all pairs of occupied sites
# We will use Haversine distance
distance <- space_dist_hav
spec <- colnames(ZP3)

#Pre-allocate a record of the average distance between all pairs of sites each species occurs at
obs_dist_avg <- rep(NA,length(spec))

for(i in 1:length(spec)){
  # Calculate the number of sites a species occurs in
  occ <- sum(ZP3[,i])
  
  # Get the row of the sites occupied by this species
  site <- which(ZP3[,i] ==1)
  
  # Place the distance matrix of all sites where this species occurs in a new matrix
  avg_dist <- distance[site,site]
  
  # Record the average distance among all pairs of sites where this species is present
  obs_dist_avg[i] <- mean(avg_dist[lower.tri(avg_dist)])
}

## Generate distribution of randomized average distance between all pairs of occupied sites, preserving number of sites where species occurs

# Pre-allocate a record of the randomized average distance between all pairs of sites each species occurs at
# Number of randomizations
N <- 1000
null_dist_avg <- array(NA,dim=c(length(spec),N),dimnames=list(spec,NULL))

for(i in 1:length(spec)){
  #Calculate the number of sites a species occurs in
  occ <- sum(ZP3[,i])
  for(j in 1:N){
    #Randomly choose which sites this species is found at
    site_rand <- sample(1:81,occ)
    
    #Place the distance matrix of all sites where this species occurs in a new matrix
    avg_dist <-distance[site_rand,site_rand]
    
    #Record the average distance among all pairs of sites where this species is present
    null_dist_avg[i,j] <- mean(avg_dist[lower.tri(avg_dist)])
  }
}

## Compare the observed value to the null distribution
# Pre-allocate a record of the average distance between all pairs of sites each species occurs at
pval <- array(NA,c(length(spec),2),dimnames=list(spec,c("greater","less")))

for(i in 1:length(spec)){
  sub <- null_dist_avg[i,]
  val <- length(sub[sub >= obs_dist_avg[i]])
  pval[i,1] <- val/N
  
  val <- length(sub[sub <= obs_dist_avg[i]])
  pval[i,2] <- val/N
}

## Three of the species show evidence of dispersal limitation, but they are not necessarily the rarest species
lim <- which(pval<= .05)  # Identifies which species are dispersal limited
spec[c(2,3,11)] # Prints name of dispersal limited species
c(sum(ZP3[,lim[1]]),sum(ZP3[,lim[2]]),sum(ZP3[,lim[3]])) # Prints number of occurences for dispersal limited species

table <- cbind(spec,pval[,1],colSums(ZP3))
colnames(table) <- c("species","p value","n")

# Test 7. Drivers of co-occurence linear model --------------------------------
# full model
mod1 <- lm(CscoreSES ~ Dscore + Geomean + Habmean + Volsum)

# null model
mod1a <- lm(CscoreSES ~ 1)

# one predictor models
mod1b <- lm(CscoreSES ~ Dscore)
mod1c <- lm(CscoreSES ~ Geomean)
mod1d <- lm(CscoreSES ~ Habmean)
mod1e <- lm(CscoreSES ~ Volsum)

# add predictors separately to Dscore model
mod1f <- lm(CscoreSES ~ Dscore + Geomean)
mod1g <- lm(CscoreSES ~ Dscore + Habmean)
mod1h <- lm(CscoreSES ~ Dscore + Volsum)

# add predictors separately to Dscore + Volsum model
mod1i <- lm(CscoreSES ~ Dscore + Volsum + Geomean)
mod1j <- lm(CscoreSES ~ Dscore + Volsum + Habmean)

AIC(mod1)
AIC(mod1a)
AIC(mod1b)
AIC(mod1c)
AIC(mod1d)
AIC(mod1e)
AIC(mod1f)
AIC(mod1g)
AIC(mod1h)
AIC(mod1i)
AIC(mod1j)

# likelihood ratio test for regression model comparison
# Single predictors compared to null model
pchisq(2 * (logLik(mod1b) - logLik(mod1a)), df = 1, lower.tail=FALSE)
pchisq(2 * (logLik(mod1c) - logLik(mod1a)), df = 1, lower.tail=FALSE)
pchisq(2 * (logLik(mod1d) - logLik(mod1a)), df = 1, lower.tail=FALSE)
pchisq(2 * (logLik(mod1e) - logLik(mod1a)), df = 1, lower.tail=FALSE)

# Single predictor (Dscore) compared to two predictors
pchisq(2 * (logLik(mod1f) - logLik(mod1b)), df = 1, lower.tail=FALSE)
pchisq(2 * (logLik(mod1g) - logLik(mod1b)), df = 1, lower.tail=FALSE)
pchisq(2 * (logLik(mod1h) - logLik(mod1b)), df = 1, lower.tail=FALSE)

# Two predictors (Dscore, Volsum) compared to three predictors
pchisq(2 * (logLik(mod1i) - logLik(mod1h)), df = 1, lower.tail=FALSE)
pchisq(2 * (logLik(mod1j) - logLik(mod1h)), df = 1, lower.tail=FALSE)

# Comparison of best supported model to null model
pchisq(2 * (logLik(mod1h) - logLik(mod1a)), df = 2, lower.tail=FALSE)

# Comparison of full model to null model
pchisq(2 * (logLik(mod1) - logLik(mod1a)), df = 4, lower.tail=FALSE)

## Post-hoc comparison of Dscore for species with significantly different (higher and lower) Cscore
sigC_lower <- which(Cscore$pval_less <= 0.05)
sigC_greater <- which(Cscore$pval_greater <= 0.05)

# t-test comparing Dscores of species pairs with lower and higher C scores than expected by chance (less or more coccurence than expected by chance)
t.test(Dscore[sigC_lower],Dscore[sigC_greater])
par(mfrow=c(1,2))
plot(Dscore,CscoreSES,pch=19,col="black")
abline(lm(CscoreSES~Dscore))
points(Dscore[sigC_lower],CscoreSES[sigC_lower],pch=19,col="blue")
points(Dscore[sigC_greater],CscoreSES[sigC_greater],pch=19,col="red")
boxplot(Dscore[sigC_lower],Dscore[sigC_greater],col=c("blue","red"),names=c("Cscore lower","Cscore higher"),ylab="Dscore")
par(mfrow=c(1,1))


# Figure 1. 2D niche use plot --------------------------------
z_uncor<-array(NA,dim=c(16,100,100),dimnames=list(sp.list,rep(paste("y",1:100,sep="")),rep(paste("x",1:100,sep=""))))

z_cor<-array(NA,dim=c(16,100,100),dimnames=list(sp.list,rep(paste("y",1:100,sep="")),rep(paste("x",1:100,sep=""))))

Z<-array(NA,dim=c(16,100,100),dimnames=list(sp.list,rep(paste("y",1:100,sep="")),rep(paste("x",1:100,sep=""))))

for (i in 1:16){
  row.sp1 <- which(occ.sp[,28] == sp.list[i]) # rows in data corresponding to sp1
  name.sp1 <- sp.list[i]
  # predict the scores on the axes
  scores.clim <- pca.cal$li[(nrow(occ.sp)+1):nrow(data),] #scores for global climate
  scores.sp1 <- pca.cal$li[row.sp1,] #scores for sp1
  
  z1 <- ecospat.grid.clim.dyn(scores.clim[,1:2], scores.clim[,1:2], scores.sp1[,1:2],R=100)
  #A grid z of RxR pixels (or a vector of R pixels) with z.uncor being the density of occurrence 
  #of the species, and z.cor the occupancy of the environment by the species 
  #(density of occurrences divided by the desinty of environment in the study area).
  z_uncor[i,,] <- t(as.matrix(z1$z.uncor))[, nrow(as.matrix(z1$z.uncor)):1]
  z_cor[i,,] <- t(as.matrix(z1$z.cor))[, nrow(as.matrix(z1$z.cor)):1]
  Z[i,,] <- t(as.matrix(z1$Z))[, nrow(as.matrix(z1$Z)):1]
}

# Contours only
z <- z_uncor[1,,]
z[z==0] <- NA

svg(paste(getwd(),"/raw_output/Fig_1.svg",sep=""),height=(984/300),width=(1476/300))
par(mar=c(3,3,2,1))
name.axis1 <- "PC-env axis 1"
name.axis2 <- "PC-env axis 2"
image(x = z1$x, y = z1$y, z = z, col = "white", zlim = c(1e-05, 1), xlab = name.axis1, ylab = name.axis2, useRaster=T)

# Graphic settings
cont_col <- rep("red",16)
cont_col[BodySize$Body_size_category == "Small"] <- "blue"

# PC-env lines
alt_name <- c("shade","sludge","snell","cond","pH","depth","O2","Chla","SPM","ALK","","DOC","tN","tP","aN","aP","SA","vSUB","","","vOV")

line <- array(NA,dim=c(2*(length(rownames(pca.cal$co))),3),dimnames=list(c(paste(alt_name,"0",sep=""),paste(alt_name,"1",sep="")),c("x","y","z")))
line[1:length(rownames(pca.cal$co)),1:3] <- 0
line[(length((rownames(pca.cal$co))) + 1):length(line[,1]),1] <- pca.cal$co[,1]/1.5
line[(length((rownames(pca.cal$co))) + 1):length(line[,1]),2] <- pca.cal$co[,2]/1.5
line[(length((rownames(pca.cal$co))) + 1):length(line[,1]),3] <- 0
ord <- c(seq(1,42,by=2),seq(2,42,by=2))
line <- line[order(ord),]

segments(0,0,line[,1],line[,2],col=rgb(0,0,0,100,maxColorValue = 255),lwd=1)

text <- cbind(pca.cal$co[,1:2]/1.5,rep(0,length(pca.cal$co[,1])))
text(text[,1],text[,2],alt_name,col="black",cex=.8,lwd=2,adj=c(1,1))

# Species contours at 97.5%
quant <- .975
for(i in 1:16){
  sub <- z_uncor[i,,]
  z <- z_uncor[i,,]
  z[z==0] <- NA
  contour(x = z1$x, y = z1$y, z =z, levels = quantile(sub[sub > 0], c(0, quant)), drawlabels = FALSE, lty = 1, col = cont_col[i], add=T, lwd=1.5)
}
dev.off()

## This is a 3-dimensional version of the niche plot (not included in the manuscript). I recommend saving your progress up to this point, as from time to time my R Studio session was aborted. Lines 697-728.
#### 3D, by body size category
#Example plot: 3D
#Create record of desired color values for each point
col_possible <- c("green","red","blue","darkmagenta","darkorange","gray","lightcoral","cyan","mediumorchid1","deepskyblue","goldenrod","navajowhite","mediumaquamarine","plum1","darkcyan","darkgreen","dimgray","firebrick4","gainsboro","darkseagreen4")
# Small
ramp_1 <- colorRamp(c(col_possible[1],"white"))
col_1 <- ramp_1(seq(0, 1, length = 100))
col_2 <- rgb(col_1, max = 255)

z <- z_uncor[1,,]
z[z==0] <- NA

persp3d(z1$x, z1$y, z, theta=50, phi=25, expand=0.75, col=rev(col_2),alpha=0.5,ticktype="detailed", xlab="", ylab="time", zlab="",axes=TRUE,NAcol=rgb(255,255,255,0,max=255),xlim=c(-5,5.4),ylim=c(-5,5.4))

small_rec <- which(BodySize$Body_size_category == "Small")
l_small <- length(small_rec)

for(i in 2:l_small){
  ramp_1 <- colorRamp(c(col_possible[small_rec[i]],"white"))
  col_1 <- ramp_1(seq(0, 1, length = 100))
  col_2 <- rgb(col_1, max = 255)
  
  z <- z_uncor[small_rec[i],,]
  z[z == 0] <- NA
  
  persp3d(z1$x, z1$y, z, theta=50, phi=25, expand=0.75, col=rev(col_2),alpha=0.5,ticktype="detailed", xlab="", ylab="time", zlab="",axes=F,NAcol=rgb(255,255,255,0,max=255),xlim=c(-5,5.4),ylim=c(-5,5.4),add=T)
}

for(i in 1:16){
# for easier visualization, use
#for(i in 1:8){  
  sub <- z_uncor[i,,]
  contour(x = z1$x, y = z1$y, z =z, add = TRUE, levels = quantile(sub[sub > 0], c(0, quant)), drawlabels = FALSE, lty = c(1,2), col = "black")
}

alt_name3d <- c("shade","sludge","snell","cond","pH","depth","O2","Chla","SPM","ALK","","DOC","tN","tP","aN","aP","SA","vSUB","","","vOV")
line <- array(NA,dim=c(2*(length(rownames(pca.cal$co))),3),dimnames=list(c(paste(alt_name,"0",sep=""),paste(alt_name,"1",sep="")),c("x","y","z")))
line[1:length(rownames(pca.cal$co)),1:3] <- 0
line[(length((rownames(pca.cal$co))) + 1):length(line[,1]),1] <- pca.cal$co[,1]/1.5
line[(length((rownames(pca.cal$co))) + 1):length(line[,1]),2] <- pca.cal$co[,2]/1.5
line[(length((rownames(pca.cal$co))) + 1):length(line[,1]),3] <- 0
ord <- c(seq(1,42,by=2),seq(2,42,by=2))
line <- line[order(ord),]
segments3d(line,col="black",add=T,lwd=3)
text <- cbind(pca.cal$co[,1:2]/1.5,rep(0,length(pca.cal$co[,1])))
text3d(text[,1],text[,2],text[,3],alt_name3d,col="black",cex=1,lwd=3,add=T)


# Figure 2. LDA density plot of body size ---------------------------------
svg(paste(getwd(),"/raw_output/Fig_2.svg",sep=""),height=(984/300),width=(1476/300))

par(mar=c(4, 4, 0, 2) + 0.1)
x1 <- dfa_results$x[BodySize$Body_size_category=="Large"]
x2 <- dfa_results$x[BodySize$Body_size_category=="Small"]
x <-  dfa_results$x

hist(x,freq=F,breaks=20,ylim=c(0,.8),xlim=c(-2,4),axes=F,xlab="LDA 1",ylab="Density",border="white",main="",col="white")
axis(side=2, at=c(0,.2,.4,.6,.8))
axis(side=1, at=c(-2,-1,0,1,2,3), pos=-.08)

d1 <- density(x1)
d2 <- density(x2)
lines(d1)
lines(d2)
polygon(d1,col=rgb(255,0,0,75,max=255))
polygon(d2,col=rgb(0,0,255,75,max=255))

segments(x1,rep(0,length(x1)),x1,rep(-.05,length(x1)),col="red")
segments(x2,rep(0,length(x2)),x2,rep(-.05,length(x2)),col="blue")

#large
text(x1,rep(-.005,length(x1)),labels=sp.list[c(6:9,15:16)],srt=45,adj=c(-.1,0),cex=.6,col="red")
#small
text(x2,rep(-.005,length(x2)),labels=sp.list[c(1:5,10:14)],srt=45,adj=c(-.1,0),cex=.6,col="blue")

dev.off()

## In case you are interested to see how the environmental variables scale on the LD1 axis. This code is derived from Pantel et al. 2022, Ecological Monographs. Please consider citing if this ends up being used in a publication (if not, I hope you are citing this code repository at least!).
## Scaling of env vars with LD1
x1 <- -10*dfa$scaling[1]
y1 <- -10*dfa$scaling[2]
x2 <- 10*dfa$scaling[1]
y2 <- 10*dfa$scaling[2]
# Here you have the slope and intercept to add a line onto the PCA biplot that indicates where the DFA line lies
slope <- (y2 - y1)/(x2-x1) 
intercept <- y1 - slope*x1
env_points <- pca.cal$co[,1:2]
# Here you would have the x-y coordinates to add the location of a given environmental variable to the DFA line.
shade <- c(pca.cal$co[1,1]*dfa$scaling[1],pca.cal$co[1,2]*dfa$scaling[2])


# Figure 3. Rural and urban niche shift plots, 50m and 3200m scale --------
for (i in 1:length(sp.list)){
  # 50m scale
  scores.clim<- pca.cal$li[(nrow(occ.sp)+1):nrow(data),] #scores for global climate
  scores.clim1<- pca.cal$li[((nrow(occ.sp))+which(cat50=="lw")),] #scores for rural climate
  scores.clim2<- pca.cal$li[((nrow(occ.sp))+which(cat50=="hg")),] #scores for urban climate
  
  row.sp1<-which(occ.sp[,28] == sp.list[i] & occ.sp[,26]=="low") # rows in data corresponding to rural
  row.sp2<-which(occ.sp[,28] == sp.list[i] & occ.sp[,26]=="high") # rows in data corresponding to urban
  scores.sp1<- pca.cal$li[row.sp1,] #scores for rural
  scores.sp2<- pca.cal$li[row.sp2,] #scores for urban
  
  if(length(row.sp1)>=5 & length(row.sp2)>=5){
    #calculate z-scores for high and low urbanization separately
    z1<- ecospat.grid.clim.dyn(scores.clim, scores.clim1, scores.sp1,R=100)
    z2<- ecospat.grid.clim.dyn(scores.clim, scores.clim2, scores.sp2,R=100)
    
    #plot niche shift from low to high urbanization
    svg(paste(getwd(),"/raw_output/Fig_3/Fig_3_",paste(substr(strsplit(sp.list[i]," ")[[1]][1],1,3),substr(strsplit(sp.list[i]," ")[[1]][2],1,2),sep="."),"_50.svg",sep=""),height=8,width=12)
    
    ecospat.plot.niche.dyn (z1, z2,quant=0.95, interest=1,title="",name.axis1="PC1",name.axis2="PC2",colinter="lightgray")
    
    ecospat.shift.centroids(scores.sp1, scores.sp2, scores.clim1, scores.clim2)
    dev.off()
  }
  
  # 3200m scale
  scores.clim<- pca.cal$li[(nrow(occ.sp)+1):nrow(data),] #scores for global climate
  scores.clim1<- pca.cal$li[((nrow(occ.sp))+which(cat3200=="lw")),] #scores for global climate
  scores.clim2<- pca.cal$li[((nrow(occ.sp))+which(cat3200=="hg")),]
  
  row.sp1<-which(occ.sp[,28] == sp.list[i] & occ.sp[,27]=="low") # rows in data corresponding to rural
  row.sp2<-which(occ.sp[,28] == sp.list[i] & occ.sp[,27]=="high") # rows in data corresponding to urban
  scores.sp1<- pca.cal$li[row.sp1,] #scores for rural
  scores.sp2<- pca.cal$li[row.sp2,] #scores for urban
  
  if(length(row.sp1)>=5 & length(row.sp2)>=5){
    #calculate z-scores for high and low urbanization separately
    z1<- ecospat.grid.clim.dyn(scores.clim, scores.clim1, scores.sp1,R=100)
    z2<- ecospat.grid.clim.dyn(scores.clim, scores.clim2, scores.sp2,R=100)
    
    #plot niche shift from low to high urbanization
    svg(paste(getwd(),"/raw_output/Fig_3/Fig_3_",paste(substr(strsplit(sp.list[i]," ")[[1]][1],1,3),substr(strsplit(sp.list[i]," ")[[1]][2],1,2),sep="."),"_3200.svg",sep=""),height=8,width=12)
    
    ecospat.plot.niche.dyn (z1, z2,quant=0.95, interest=1,title="",name.axis1="PC1",name.axis2="PC2",colinter="lightgray")
    
    ecospat.shift.centroids(scores.sp1, scores.sp2, scores.clim1, scores.clim2)
    dev.off()
  }
}

# Figure 4. Scatterplot of CscoreSES and Dscore ---------------------------
svg(paste(getwd(),"/raw_output/Fig_4.svg",sep=""),height=(945/300),width=(945/300),pointsize=8)
par(mar=c(4,4,2,2))
plot(Dscore,CscoreSES,pch=19,col="black")
abline(lm(CscoreSES~Dscore))
points(Dscore[sigC_lower],CscoreSES[sigC_lower],pch=19,col="blue")
points(Dscore[sigC_greater],CscoreSES[sigC_greater],pch=19,col="red")
dev.off()

# Figure S1. Map of Belgium w/ site richness --------------------------------------
## Site map colored by species richness
p <- readOGR(paste(getwd(),"/gadm36_BEL_shp/gadm36_BEL_2.shp",sep=""))   #layer of data for species range
pal <- colorRampPalette(c("white", "blue"))

svg(filename=paste(getwd(),"/raw_output/Fig_S1.svg",sep=""),height=8,width=8)
options(warn=-1)
plot(p,col="grey90", fill=TRUE, warnings=F)
options(warn=0)
points(spa$LONG, spa$LAT,pch = 21,bg = pal(9)[rowSums(ZP2)],col="black")
legend("topright", legend = levels(as.factor(rowSums(ZP2))),pt.cex = 2,pt.bg=pal(9),pch = 21,bg="transparent",bty="n")
scalebar(10,type="bar")
dev.off()

# Figure S2. Distribution of Scores ---------------------------------------
svg(filename=paste(getwd(),"/raw_output/Fig_S2.svg",sep=""),height=2.05,width=3.15,pointsize=5)
par(mfrow=c(2,2),mar = c(4, 4, 0.1, 0.1),lwd=.5)
denscomp(Cgam,fitlwd=.5)
denscomp(SESCnorm,fitlwd=.5)
cdfcomp(Cgam,fitlwd=.5)
cdfcomp(SESCnorm,fitlwd=.5)
dev.off()

# Figure S3. Map of SDM habitat suitability ------------------------------
auc <- rep(NA,dim(ZP3)[2])
svg(paste(getwd(),"/raw_output/Fig_S3.svg",sep=""),height=8,width=12)
par(mfrow=c(4,4))
for(i in 1:dim(ZP3)[2]){
  use <- clim[,c(2,1)]
  use$occ <- ZP3[,i]
  me <- dismo::maxent(env_use,use[,3])
  r <- predict(me, env_use)
  auc[i] <- me@results[5]
  
  # Plot
  # Draw a spatial plot and fill in habitat suitability values.
  colfunc <- colorRampPalette(c("white", "black"))(1000)
  #palette <- colorRampPalette(c("blue","red"))(maxColorValue)
  plot(use[,1:2],col="white",bty="n",main=paste(sp.list[i]," (AUC = ",auc[i],")"))
  points(use[use[,3] == 0,1], use[use[,3] == 0,2], pch=19,col = colfunc[cut(r, 1000)][use[,3] == 0],cex=2)
  points(use[use[,3] == 1,1], use[use[,3] == 1,2], pch=15,col = colfunc[cut(r, 1000)][use[,3] == 1],cex=2)
}
dev.off()

# legend
svg(paste(getwd(),"/raw_output/Fig_S3_legend.svg",sep=""),height=2,width=2,pointsize=6)
colfunc <- colorRampPalette(c("white", "black"))(1000)
par(mfrow=c(1,1))
legend_image <- as.raster(matrix(rev(colfunc), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'habitat suitability')
text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1,l=5))
rasterImage(legend_image, 0, 0, 1,1)
dev.off()

# STOPPED HERE for image sizing #
# Figure S4. Species site occurences --------------------------------
comm <- ZP2
rownames(comm) <- 1:81
colnames(comm) <- 1:28
lis <- melt(comm)

#Plot p-a matrix
svg(paste(getwd(),"/raw_output/Fig_S4.svg",sep=""),height=8,width=12)
par(mar=c(5, 12, 4, 2) + 0.1)
plot(lis[,1],lis[,2],xlab="site",ylab="",yaxt="n")
axis(side=2,at=1:28,labels=colnames(ZP2),las=1)
points(lis[lis$value==1,1], lis[lis$value==1,2], pch=19, col="black")
dev.off()

# Figure S5. Boxplot of D scores by species pair type ---------------------
meanD <- rep(NA,length(sp.list))
for(i in 1:length(sp.list)){
  meanD[i] <- mean(D[which(apply(sp.combn, 2, function(r) any(r %in% i)))])
}

svg(filename=paste(getwd(),"/raw_output/Fig_S5.svg",sep=""),height=4.725,width=3.15,pointsize=5)
par(mfrow=c(2,1))
boxplot(d_anova$D ~ d_anova$pair,ylab="D score",xaxt="n",cex.axis=1.5)
plot(sp.nbocc[sp.nbocc>4],meanD,ylab = "",xlab="",bty="n",ylim=c(.2,.5),cex=2,cex.axis=1.5)
points(sp.nbocc[sp.nbocc>4][c(5,6,9,15)],meanD[c(5,6,9,15)],pch=19,col="black",cex=2)
points(sp.nbocc[sp.nbocc>4][c(7,10,12,13)],meanD[c(7,10,12,13)],pch=19,col="orange",cex=2)
dev.off()


# Figure S6. 2D niche use plot, generalists and specialists --------------------------------
#Plot 1st species, specialist
#Declare contrasting-ish colors
col_possible <- c("green","red","blue","darkmagenta","darkorange","gray","lightcoral","cyan","mediumorchid1","deepskyblue","goldenrod","navajowhite","mediumaquamarine","plum1","darkcyan","darkgreen","dimgray","firebrick4","gainsboro","darkseagreen4")
col_possible <- col_possible[1:16]
name.axis1 <- "Axis 1"
name.axis2 <- "Axis 2"

quant <- .975
ramp_1 <- colorRamp(c(col_possible[10],"white"))
col_1 <- rgb(ramp_1(seq(0, 1, length = 100)), max = 255)

svg(paste(getwd(),"/raw_output/Fig_S6.svg",sep=""),height=12,width=8)
par(mfrow=c(2,1))
image(x = z1$x, y = z1$y, z = z_uncor[10,,], col = rev(col_1), zlim = c(1e-05, 1), xlab = name.axis1, ylab = name.axis2)

#Plot all other specialist species in a loop
for (i in c(7,12,13)){
  ramp <- colorRamp(c(col_possible[i],"white"))
  col_n <- rgb(ramp(seq(0, 1, length = 100)), max = 255)
  image(x = z1$x, y = z1$y, z = z_uncor[i,,], col = rev(col_n), zlim = c(1e-05, 1), xlab = name.axis1, ylab = name.axis2, add=T)
}

for(i in c(7,10,12,13)){
  contour(x = z1$x, y = z1$y, z = z_uncor[i,,], add = TRUE, levels = quantile(z_uncor[z_uncor[i,,]>0], c(0, quant)), drawlabels = FALSE, lty = c(1,2), col = "black")
}

# Generalists
ramp_1 <- colorRamp(c(col_possible[5],"white"))
col_1 <- rgb(ramp_1(seq(0, 1, length = 100)), max = 255)
image(x = z1$x, y = z1$y, z = z_uncor[5,,], col = rev(col_1), zlim = c(1e-05, 1), xlab = name.axis1, ylab = name.axis2)

#Plot all other species in a loop
for (i in c(6,9,15)){
  ramp <- colorRamp(c(col_possible[i],"white"))
  col_n <- rgb(ramp(seq(0, 1, length = 100)), max = 255)
  image(x = z1$x, y = z1$y, z = z_uncor[i,,], col = rev(col_n), zlim = c(1e-05, 1), xlab = name.axis1, ylab = name.axis2, add=T)
}

for(i in c(5,6,9,15)){
  contour(x = z1$x, y = z1$y, z = z_uncor[i,,], add = TRUE, levels = quantile(z_uncor[z_uncor[i,,]>0], c(0, quant)), drawlabels = FALSE, lty = c(1,2), col = "black")
}
dev.off()

# Figure S7. PCenv biplot with species scores and LD1 axis ----------------
## Scale env vars with LD1
#plot scaling onto PC biplot
alt_name <- c("shade","sludge","sneller","cond","pH","depth","O2","Chla","SPM","ALK","SULF","DOC","tN","tP","aN","aP","SA","vSUB","vFLOAT","vEMERSE","vOV")

svg(paste(getwd(),"/raw_output/Fig_S7.svg",sep=""),height=8,width=12)
plot(scores.sp[,1], scores.sp[,2], col = "white",xlab = "PC1", ylab = "PC2",xlim=c(-5,5.4),ylim=c(-5,5.4))
alt.sp.name1 <- c("Dlc","Dm","Do","Dp","Se","Sv")
alt.sp.name2 <- c("Ag","Ar","Bl","Ce","Cs","Gt","Lq","Pa","Pt","Sm")
#points(scores.sp[c(6:9,15:16),1], scores.sp[c(6:9,15:16),2],pch=19,col="red",cex=2) # large
#points(scores.sp[c(1:5,10:14),1], scores.sp[c(1:5,10:14),2],pch=19,col="blue",cex=2) # small
arrows(0,0,pca.cal$co[,1]/1.5,pca.cal$co[,2]/1.5,col="darkgray")
text(pca.cal$co/1.5, alt_name,col="black",pos=3)
segments((-10*dfa$scaling[1]), (-10*dfa$scaling[2]), (10*dfa$scaling[1]), (10*dfa$scaling[2]), lwd=2,col="black")
text(scores.sp[c(6:9,15:16),1], scores.sp[c(6:9,15:16),2],labels=alt.sp.name1,col="red")
text(scores.sp[c(1:5,10:14),1], scores.sp[c(1:5,10:14),2],labels=alt.sp.name2,col="blue")
dev.off()




# Figure S8. PCA for rural vs. urban niche axes ----------------
for_dfapc_50 <- clim[c(which(cat50=="lw"), which(cat50=="hg")),3:23]
for_dfapc_3200 <- clim[c(which(cat3200=="lw"), which(cat3200=="hg")),3:23]
pc_50 <- prcomp(for_dfapc_50,center=T,scale=T)
pc_3200 <- prcomp(for_dfapc_3200,center=T,scale=T)

urban_rec <- cat50[c(which(cat50=="lw"), which(cat50=="hg"))]
urban_rec_3200 <- cat3200[c(which(cat3200=="lw"), which(cat3200=="hg"))]
# Plot of groupings
svg(paste(getwd(),"/raw_output/Fig_S8.svg",sep=""),height=8,width=8)
par(mfrow=c(2,1))
# 50 m
plot(pc_50$x[,1:2],type="n",bty="none",main="50m",font.main=1,family.main="sans")
points(pc_50$x[which(urban_rec=="lw"),1:2],,pch=19,col=rgb(34,139,34,max=255))
points(pc_50$x[which(urban_rec=="hg"),1:2],,pch=19,col=rgb(255,0,0,max=255))
x1 <- chull(pc_50$x[which(urban_rec=="lw"),1:2])
x1 <- c(x1,x1[1])
lines(pc_50$x[x1,1:2])
polygon(pc_50$x[x1,1:2],col=rgb(34,139,34,75,max=255))

x2 <- chull(pc_50$x[which(urban_rec=="hg"),1:2])
x2 <- c(x2,x2[1])
lines(pc_50$x[which(urban_rec=="hg")[x2],1:2])
polygon(pc_50$x[which(urban_rec=="hg")[x2],1:2],col=rgb(255,0,0,75,max=255))

# 3200 m
plot(pc_3200$x[,1:2],type="n",bty="none",main="3200m",font.main=1,family.main="sans")
points(pc_3200$x[which(urban_rec_3200=="lw"),1:2],,pch=19,col=rgb(34,139,34,max=255))
points(pc_3200$x[which(urban_rec_3200=="hg"),1:2],,pch=19,col=rgb(255,0,0,max=255))
x1 <- chull(pc_3200$x[which(urban_rec_3200=="lw"),1:2])
x1 <- c(x1,x1[1])
lines(pc_3200$x[x1,1:2])
polygon(pc_3200$x[x1,1:2],col=rgb(34,139,34,75,max=255))

x2 <- chull(pc_3200$x[which(urban_rec_3200=="hg"),1:2])
x2 <- c(x2,x2[1])
lines(pc_3200$x[which(urban_rec_3200=="hg")[x2],1:2])
polygon(pc_3200$x[which(urban_rec_3200=="hg")[x2],1:2],col=rgb(255,0,0,75,max=255))

legend(2,4,legend=c("rural","urban"),pch=19,col=c("forestgreen","red"),bty="n",cex=1.5,y.intersp=.5)

dev.off()

# Figure S9. Nestedness temperature plot ----------------------------------
# Points depict scale of %BUA (Built Up Area)
#Import BUA data
BUA <- read.csv("BUA.csv")
gray_pal <- colorRampPalette(c("white","black"))
BUA$Col_50 <- gray_pal(50)[as.numeric(cut(BUA$X50,breaks = 50))]
BUA$Col_3200 <- gray_pal(50)[as.numeric(cut(BUA$X3200,breaks = 50))]

#re-arrange BUA values to match the ordering in res
x <- as.character(BUA$PONDCODE)
y <- rownames(res$u)
order_50 <- BUA$Col_50[order(match(x,y))]
order_3200 <- BUA$Col_3200[order(match(x,y))]

#Color text of species name by body size category
lab_cols <- as.data.frame(cbind(names(res$c),rep(NA,length(names(res$c)))))
colnames(lab_cols) <- c("spec","size")
lab_cols$size <- factor(lab_cols$size,levels=c("small","large"))
lab_cols$size <- c("small","large","large","large","small","small","large","small","small","large","small","small","large","small","small","small")
lab_cols$color <- rep(NA,dim(lab_cols)[1])
lab_cols$color[lab_cols$size=="small"] <- "blue"
lab_cols$color[lab_cols$size=="large"] <- "red"

svg(paste(getwd(),"/raw_output/Fig_S9.svg",sep=""),height=12,width=12)
par(mar=c(5,4,6,2)+0.1)
x <- res
z <- x$u
z <- t(z[nrow(z):1,])
col = rev(heat.colors(100))
image(z,axes=F,col=col,ylim=c(-.01,1.3),xlim=c(-.06,1.2))
lines(x$smooth$x, 1 - x$smooth$y)

z <- res$comm
z <- t(z[nrow(z):1,])
points(rep(-.05,length((seq(1, 0, len = ncol(z))))),seq(1, 0, len = ncol(z)),pch=19,col=order_50,cex=.9)
points(rep(-.03,length((seq(1, 0, len = ncol(z))))),seq(1, 0, len = ncol(z)),pch=19,col=order_3200,cex=.9)
text(seq(1, 0, len = nrow(z)), rep(1,length((seq(1, 0, len = nrow(z))))),labels=colnames(res$comm),col=lab_cols$color,srt=45,adj=c(-.1,.5))
text(c(-.05,-.03),c(1,1),labels=c("BUA 50","BUA 3200"),col="black",srt=45,adj=c(-.1,.1))
dev.off()

# Figure S10. Niche plots for species pairs with limiting similarity -------
### Chydorus sphaericus and Daphnia pulex
row.sp1 <- which(occ.sp[,28] == sp.list[5]) # rows in data corresponding to sp1
row.sp2 <- which(occ.sp[,28] == sp.list[9]) # rows in data corresponding to sp1
# predict the scores on the axes
scores.clim<- pca.cal$li[(nrow(occ.sp)+1):nrow(data),] #scores for global climate
scores.sp1<- pca.cal$li[row.sp1,] #scores for sp1
scores.sp2<- pca.cal$li[row.sp2,] #scores for sp2
# calculation of occurence density and test of niche equivalency and similarity
z1<- ecospat.grid.clim.dyn(scores.clim, scores.clim, scores.sp1,R=100)
z2<- ecospat.grid.clim.dyn(scores.clim, scores.clim, scores.sp2,R=100)

#plot niche differences
svg(paste(getwd(),"/raw_output/Fig_S10a.svg",sep=""),height=8,width=12)
ecospat.plot.niche.dyn (z1, z2,quant=1, title=paste("niche of ", sp.list[5], " and ", sp.list[9],sep=""),name.axis1="PC1",name.axis2="PC2",colinter="lightgray")
points(scores.sp1,col="green3",pch=19)
points(scores.sp2,col="red3",pch=19)

segments(0,0,line[,1],line[,2],col=rgb(0,0,0,100,maxColorValue = 255),lwd=1)

text(text[,1],text[,2],alt_name,col="black",cex=1,lwd=3,adj=c(1,1))
dev.off()

### Chydorus sphaericus and Daphnia longispina complex
row.sp1 <- which(occ.sp[,28] == sp.list[5]) # rows in data corresponding to sp1
row.sp2 <- which(occ.sp[,28] == sp.list[6]) # rows in data corresponding to sp1
# predict the scores on the axes
scores.clim<- pca.cal$li[(nrow(occ.sp)+1):nrow(data),] #scores for global climate
scores.sp1<- pca.cal$li[row.sp1,] #scores for sp1
scores.sp2<- pca.cal$li[row.sp2,] #scores for sp2
# calculation of occurence density and test of niche equivalency and similarity
z1<- ecospat.grid.clim.dyn(scores.clim, scores.clim, scores.sp1,R=100)
z2<- ecospat.grid.clim.dyn(scores.clim, scores.clim, scores.sp2,R=100)

#plot niche differences
svg(paste(getwd(),"/raw_output/Fig_S10b.svg",sep=""),height=8,width=12)
ecospat.plot.niche.dyn (z1, z2,quant=1, title=paste("niche of ", sp.list[5], " and ", sp.list[6],sep=""),name.axis1="PC1",name.axis2="PC2",colinter="lightgray")
points(scores.sp1,col="green3",pch=19)
points(scores.sp2,col="red3",pch=19)
segments(0,0,line[,1],line[,2],col=rgb(0,0,0,100,maxColorValue = 255),lwd=1)
text(text[,1],text[,2],alt_name,col="black",cex=1,lwd=3,adj=c(1,1))
dev.off()

### Alona guttata and Daphnia pulex
row.sp1 <- which(occ.sp[,28] == sp.list[1]) # rows in data corresponding to sp1
row.sp2 <- which(occ.sp[,28] == sp.list[9]) # rows in data corresponding to sp1
# predict the scores on the axes
scores.clim<- pca.cal$li[(nrow(occ.sp)+1):nrow(data),] #scores for global climate
scores.sp1<- pca.cal$li[row.sp1,] #scores for sp1
scores.sp2<- pca.cal$li[row.sp2,] #scores for sp2
# calculation of occurence density and test of niche equivalency and similarity
z1<- ecospat.grid.clim.dyn(scores.clim, scores.clim, scores.sp1,R=100)
z2<- ecospat.grid.clim.dyn(scores.clim, scores.clim, scores.sp2,R=100)

#plot niche differences
svg(paste(getwd(),"/raw_output/Fig_S10c.svg",sep=""),height=8,width=12)
ecospat.plot.niche.dyn (z1, z2,quant=1, title=paste("niche of ", sp.list[1], " and ", sp.list[9],sep=""),name.axis1="PC1",name.axis2="PC2",colinter="lightgray")
points(scores.sp1,col="green3",pch=19)
points(scores.sp2,col="red3",pch=19)
segments(0,0,line[,1],line[,2],col=rgb(0,0,0,100,maxColorValue = 255),lwd=1)
text(text[,1],text[,2],alt_name,col="black",cex=1,lwd=3,adj=c(1,1))
dev.off()

########## Extra material for revisions ##########
### 1. Plots of rural and urban environment at 50m and 3200m scale
scores.clim<- pca.cal$li[(nrow(occ.sp)+1):nrow(data),] #scores for global climate
scores.clim1<- pca.cal$li[((nrow(occ.sp))+which(cat50=="lw")),] #scores for rural climate
scores.clim2<- pca.cal$li[((nrow(occ.sp))+which(cat50=="hg")),] #scores for urban climate
z1<- ecospat.grid.clim.dyn(scores.clim, scores.clim, scores.clim1,R=100)
z2<- ecospat.grid.clim.dyn(scores.clim, scores.clim, scores.clim2,R=100)
ecospat.plot.niche.dyn (z1, z2,quant=1, title="Rural (green) and urban (red) niche, 3200m",name.axis1="PC1",name.axis2="PC2",colinter="lightgray")
segments(0,0,line[,1],line[,2],col=rgb(0,0,0,100,maxColorValue = 255),lwd=1)
text(text[,1],text[,2],alt_name,col="black",cex=1,lwd=3,adj=c(1,1))

scores.clim1<- pca.cal$li[((nrow(occ.sp))+which(cat3200=="lw")),] #scores for rural climate
scores.clim2<- pca.cal$li[((nrow(occ.sp))+which(cat3200=="hg")),] #scores for urban climate
## Boxplot of submerged vegetation
boxplot(env$submerse[cat50 == "lw"],env$submerse[cat50 == "hg"],main="vSUB at 50m (left is rural, right is urban)")
boxplot(env$submerse[cat3200 == "lw"],env$submerse[cat3200 == "hg"],main="vSUB at 3200m (left is rural, right is urban)")

### 2. Evaluation of z values and their relationship with overall occupancy
sum_z_list <- rep(NA,length(sp.list))
for(i in 1:length(sp.list)){
  sp.list[i]
  row.sp1<-which(occ.sp[,28] == sp.list[i])
  scores.clim<- pca.cal$li[(nrow(occ.sp)+1):nrow(data),] #scores for global climate
  scores.sp1<- pca.cal$li[row.sp1,] #scores for sp1
  z1<- ecospat.grid.clim.dyn(scores.clim, scores.clim, scores.sp1,R=100)
  sum_z_list[i] <- sum(z1$z.cor[,,1])
}
# We can see that total z *is* associated strongly with the number of occupied sites.
plot(sum_z_list,sp.nbocc[sp.nbocc>4])

# Next we check whether D is associated with the difference in the number of occupied sites between pairs of species
d12 <- rep(NA,dim(sp.combn)[2])
for(i in 1:(ncol(sp.combn))) {
  nm1 <- sp.nbocc[sp.nbocc>4][sp.combn[1,i]]
  nm2 <- sp.nbocc[sp.nbocc>4][sp.combn[2,i]]
  d12[i] <- abs(nm1-nm2)
}
  
plot(d12,D)
# Color points by whether they are gen-gen, gen-spec, spec-spec species pairs
points(d12[ss_pairs],D[ss_pairs],pch=19,col="black")
points(d12[gg_pairs],D[gg_pairs],pch=19,col="orange")
points(d12[sg_pairs],D[sg_pairs],pch=19,col="skyblue")

# Now I check whether a species average D is associated with its # sites occupied
meanD <- rep(NA,length(sp.list))
for(i in 1:length(sp.list)){
  meanD[i] <- mean(D[which(apply(sp.combn, 2, function(r) any(r %in% i)))])
}

plot(sp.nbocc[sp.nbocc>4],meanD,ylab = "mean D",xlab="number of site occurences",bty="n",ylim=c(.2,.5))
points(sp.nbocc[sp.nbocc>4][c(5,6,8,9)],meanD[c(5,6,8,9)],pch=19,col="black")
points(sp.nbocc[sp.nbocc>4][c(7,10:12)],meanD[c(7,10:12)],pch=19,col="orange")
# This is now incorporated into Figure S2.

### 3. Stats analysis whether sites differed in environmental variables.
mod50 <- betadisper(dist(scale(for_dfapc_50)),urban_rec)
anova(mod50)
mod3200 <- betadisper(dist(scale(for_dfapc_3200)),urban_rec_3200)
anova(mod3200)