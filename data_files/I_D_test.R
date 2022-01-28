# I and D are observed values of I and D scores for all species pairs
I <- rep(NA,ncol(sp.combn))
D <- rep(NA,ncol(sp.combn))
# test, a 2x2x2x120 array that holds the p values of niche equivalency and similarity tests, for greater than and lower than, for the I and D measures of niche overlap.
test <- array(NA,dim=c(2,2,2,120),dimnames=list(c("equivalency","similarity"),c("greater","lower"),c("pI","pD"),NULL))

for(i in 1:(ncol(sp.combn))) {
  row.sp1<-which(occ.sp[,28] == sp.list[sp.combn[1,i]]) # rows in data corresponding to sp1
  row.sp2<-which(occ.sp[,28] == sp.list[sp.combn[2,i]]) # rows in data corresponding to sp2
  name.sp1<-sp.list[sp.combn[1,i]]
  name.sp2<-sp.list[sp.combn[2,i]]
  
  # predict the scores on the axes
  scores.clim<- pca.cal$li[(nrow(occ.sp)+1):nrow(data),] #scores for global climate
  scores.sp1<- pca.cal$li[row.sp1,] #scores for sp1
  scores.sp2<- pca.cal$li[row.sp2,] #scores for sp2
  
  # calculation of occurence density and test of niche equivalency and similarity
  z1<- ecospat.grid.clim.dyn(scores.clim, scores.clim, scores.sp1,R=100)
  z2<- ecospat.grid.clim.dyn(scores.clim, scores.clim, scores.sp2,R=100)
  
  ## calculate the overlap metrics D and I (see Warren et al 2008) based on two species occurrence density grids z1 and z2 created by grid.clim
  ## cor=T correct occurrence densities of each species by the prevalence of the environments in their range
  ecospat.niche.overlap(z1,z2,cor=T)
  
  ## runs niche equivalency test(see Warren et al 2008) based on two species occurrence density grids
  ## compares the observed niche overlap between z1 and z2 to overlaps between random niches z1.sim and z2.sim.
  ## z1.sim and z2.sim are built from random reallocations of occurences of z1 and z2
  ## rep is the number of iterations
  x1<-ecospat.niche.equivalency.test(z1,z2,rep=100,alternative = "greater")
  x2<-ecospat.niche.equivalency.test(z1,z2,rep=100,alternative = "lower")
  
  I[i] <- x1$obs$I
  D[i] <- x1$obs$D
  test[1,1,1,i] <- x1$p.I   # equivalency, greater, p-value of I
  test[1,1,2,i] <- x1$p.D   # equivalency, greater, p-value of D
  test[1,2,1,i] <- x2$p.I   # equivalency, lower, p-value of I
  test[1,2,2,i] <- x2$p.D   # equivalency, lower, p-value of D
  
  ## runs niche similarity test(see Warren et al 2008) based on two species occurrence density grids
  ## compares the observed niche overlap between z1 and z2 to overlaps between z1 and random niches (z2.sim) in available in the range of z2 (z2$Z) 
  ## z2.sim have the same patterns as z2 but their center are randomly translatated in the availabe z2$Z space and weighted by z2$Z densities
  ## rep is the number of iterations
  x3 <- ecospat.niche.similarity.test(z1,z2,rep=100,rand.type=1, alternative = "greater")
  x4 <- ecospat.niche.similarity.test(z1,z2,rep=100,rand.type=1, alternative = "lower")
  
  test[2,1,1,i] <- x3$p.I   # similarity, greater, p-value of I
  test[2,1,2,i] <- x3$p.D   # similarity, greater, p-value of D
  test[2,2,1,i] <- x4$p.I   # similarity, lower, p-value of I
  test[2,2,2,i] <- x4$p.D   # similarity, lower, p-value of D
}
