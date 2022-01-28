#FOR LOCAL URBANIZATION / 50m radius
#calculate scores of entire study, of only low and only high for 50m urbanization
scores.clim<- pca.cal$li[(nrow(occ.sp)+1):nrow(data),] #scores for global climate
scores.clim1<- pca.cal$li[((nrow(occ.sp))+which(cat50=="lw")),] #scores for rural climate
scores.clim2<- pca.cal$li[((nrow(occ.sp))+which(cat50=="hg")),] #scores for urban climate

# Record for results of all tests
rural_urban_50 <- data.frame(sp=character(),D_score=numeric(),equiv_pval_low=numeric(),equiv_pval_high=numeric(),sim_pval_low=numeric(),sim_pval_high=numeric(),expansion=numeric(),stability=numeric(),unfilling=numeric())

for (i in 1:length(sp.list)){
  row.sp1<-which(occ.sp[,28] == sp.list[i] & occ.sp[,26]=="low") # rows in data corresponding to rural
  row.sp2<-which(occ.sp[,28] == sp.list[i] & occ.sp[,26]=="high") # rows in data corresponding to urban
  scores.sp1<- pca.cal$li[row.sp1,] #scores for rural
  scores.sp2<- pca.cal$li[row.sp2,] #scores for urban
  
  if(length(row.sp1)>=5 & length(row.sp2)>=5){
    df <- array(NA,c(1,9))
    df[1] <- as.character(sp.list[i])
    
    #calculate z-scores for high and low urbanization separately
    z1<- ecospat.grid.clim.dyn(scores.clim, scores.clim1, scores.sp1,R=100)
    z2<- ecospat.grid.clim.dyn(scores.clim, scores.clim2, scores.sp2,R=100)
    
    #calculate overlap (D score)
    df[2] <- ecospat.niche.overlap(z1,z2,cor=T)$D
    
    x1 <- ecospat.niche.equivalency.test(z1,z2,rep=100,alternative = "lower")
    df[3] <- x1$p.D
    
    x2<-ecospat.niche.equivalency.test(z1,z2,rep=100,alternative = "greater")
    df[4] <- x2$p.D
    
    x3 <- ecospat.niche.similarity.test(z1,z2,rep=100,rand.type=1, alternative = "lower")
    df[5] <- x3$p.D
    x4 <- ecospat.niche.similarity.test(z1,z2,rep=100,rand.type=1, alternative = "greater")
    df[6] <- x4$p.D
    
    #plot niche shift from low to high urbanization
    niche.dyn <- ecospat.niche.dyn.index (z1, z2, intersection = 0.1)
    df[7:9] <- niche.dyn$dynamic.index.w
    ecospat.plot.niche.dyn (z1, z2,quant=0.95, interest=1,title=paste("niche of ", sp.list[i], "for low and high urbanization"),name.axis1="PC1",name.axis2="PC2")
    
    ecospat.shift.centroids(scores.sp1, scores.sp2, scores.clim1, scores.clim2)
    rural_urban_50 <- rbind(rural_urban_50,df)
  }
}
colnames(rural_urban_50) <- c("sp","D_score","equiv_pval_low","equiv_pval_high","sim_pval_low","sim_pval_high","expansion","stability","unfilling")

#FOR REGIONAL URBANIZATION
#calculate scores of entire study, of only low and only high for 50m urbanization
scores.clim<- pca.cal$li[(nrow(occ.sp)+1):nrow(data),] #scores for global climate
scores.clim1<- pca.cal$li[((nrow(occ.sp))+which(cat3200=="lw")),] #scores for global climate
scores.clim2<- pca.cal$li[((nrow(occ.sp))+which(cat3200=="hg")),]

# Record for results of all tests
rural_urban_3200 <- data.frame(sp=character(),D_score=numeric(),equiv_pval_low=numeric(),equiv_pval_high=numeric(),sim_pval_low=numeric(),sim_pval_high=numeric(),expansion=numeric(),stability=numeric(),unfilling=numeric())

for (i in 1:length(sp.list)){
  row.sp1<-which(occ.sp[,28] == sp.list[i] & occ.sp[,27]=="low") # rows in data corresponding to rural
  row.sp2<-which(occ.sp[,28] == sp.list[i] & occ.sp[,27]=="high") # rows in data corresponding to urban
  scores.sp1<- pca.cal$li[row.sp1,] #scores for rural
  scores.sp2<- pca.cal$li[row.sp2,] #scores for urban
  
  if(length(row.sp1)>=5 & length(row.sp2)>=5){
    df <- array(NA,c(1,9))
    df[1] <- as.character(sp.list[i])
    
    #calculate z-scores for high and low urbanization separately
    z1<- ecospat.grid.clim.dyn(scores.clim, scores.clim1, scores.sp1,R=100)
    z2<- ecospat.grid.clim.dyn(scores.clim, scores.clim2, scores.sp2,R=100)
    
    #calculate overlap (D score)
    df[2] <- ecospat.niche.overlap(z1,z2,cor=T)$D
    
    x1 <- ecospat.niche.equivalency.test(z1,z2,rep=100,alternative = "lower")
    df[3] <- x1$p.D
    
    x2<-ecospat.niche.equivalency.test(z1,z2,rep=100,alternative = "greater")
    df[4] <- x2$p.D
    
    x3 <- ecospat.niche.similarity.test(z1,z2,rep=100,rand.type=1, alternative = "lower")
    df[5] <- x3$p.D
    x4 <- ecospat.niche.similarity.test(z1,z2,rep=100,rand.type=1, alternative = "greater")
    df[6] <- x4$p.D
    
    #plot niche shift from low to high urbanization
    niche.dyn <- ecospat.niche.dyn.index (z1, z2, intersection = 0.1)
    df[7:9] <- niche.dyn$dynamic.index.w
    ecospat.plot.niche.dyn (z1, z2,quant=0.1, interest=1,title=paste("niche of ", sp.list[i], "for low and high urbanization"),name.axis1="PC1",name.axis2="PC2")
    
    ecospat.shift.centroids(scores.sp1, scores.sp2, scores.clim1, scores.clim2)
    rural_urban_3200 <- rbind(rural_urban_3200,df)
  }
}

colnames(rural_urban_3200) <- c("sp","D_score","equiv_pval_low","equiv_pval_high","sim_pval_low","sim_pval_high","expansion","stability","unfilling")

## Note the intersect value is 0.1, which corresponds to a 90% inclusion of the two types of habitats.