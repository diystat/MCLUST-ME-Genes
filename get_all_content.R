### Run the following function to obtain all graphs and tables in the paper
get.paper.content = function(){
  ## First source all relative functions
  source("all_core_functions.R")
  source("all_other_functions.R")
  
  #------ Figure 1 ------#
  # Import clustering results
  setwd("/Users/wzhang/Research Project/Simulations for Paper/Example 1/p=0.5/Results")
  load("res.RData")
  tmp = out$sim.result
  
  simrun7 = list(res.mcme=tmp[[22]],res.mevvv=tmp[[23]],z.ini=tmp[$z.ini][24]],
    rand.samples=tmp[[25]],errmat=tmp[[26]],index=tmp[[27]],k=tmp[[28]])

  par(mfrow=c(1,2))
  plot.boundary(simrun7,7)
  par(mfrow=c(1,1))
  
  
  #------ Figure 2 ------#
  mcme7 = tmp[[22]]
  mclust7 = tmp[[23]]
  unc.mcme = mcme7$uncertainty
  unc.mclust = numeric(200)
  for(i in 1:200){
    unc.mclust[i] = 1-max(mclust7$z[i,])
  }
  point.size = 0.3+unc.mcme*4
  point.size.mclust = 0.3+unc.mclust*4

  plot.boundary.new(simrun7,point.size,point.size.mclust)
  
  
  #------ Figure 3 ------#
  # Import simulation results
  setwd("/Users/wzhang/Research Project/Simulations for Paper/Example 1/p=0.1/Results")
  load("res.RData")
  rr = out$rand.raw
  r1.mcme = rr[,1]
  r1.mevvv = rr[,4]
  fr1.mcme = rr[,7]
  fr1.mevvv = rr[,8]

  setwd("/Users/wzhang/Research Project/Simulations for Paper/Example 1/p=0.3/Results")
  load("res.RData")
  rr = out$rand.raw
  r3.mcme = rr[,1]
  r3.mevvv = rr[,4]
  fr3.mcme = rr[,7]
  fr3.mevvv = rr[,8]

  setwd("/Users/wzhang/Research Project/Simulations for Paper/Example 1/p=0.5/Results")
  load("res.RData")
  rr = out$rand.raw
  r5.mcme = rr[,1]
  r5.mevvv = rr[,4]
  fr5.mcme = rr[,7]
  fr5.mevvv = rr[,8]

  setwd("/Users/wzhang/Research Project/Simulations for Paper/Example 1/p=0.7/Results")
  load("res.RData")
  rr = out$rand.raw
  r7.mcme = rr[,1]
  r7.mevvv = rr[,4]
  fr7.mcme = rr[,7]
  fr7.mevvv = rr[,8]

  setwd("/Users/wzhang/Research Project/Simulations for Paper/Example 1/p=0.9/Results")
  load("res.RData")
  rr = out$rand.raw
  r9.mcme = rr[,1]
  r9.mevvv = rr[,4]
  fr9.mcme = rr[,7]
  fr9.mevvv = rr[,8]

  rand.mcme = c(r1.mcme,r3.mcme,r5.mcme,r7.mcme,r9.mcme)
  rand.mevvv = c(r1.mevvv,r3.mevvv,r5.mevvv,r7.mevvv,r9.mevvv)
  frand.mcme = c(fr1.mcme,fr3.mcme,fr5.mcme,fr7.mcme,fr9.mcme)
  frand.mevvv = c(fr1.mevvv,fr3.mevvv,fr5.mevvv,fr7.mevvv,fr9.mevvv)

  group = c(rep(0.1,100),rep(0.3,100),rep(0.5,100),rep(0.7,100),rep(0.9,100))

  # Individual Rand/Fuzzy Rand indices
  par(mfrow=c(2,2))
  boxplot(rand.mcme~group,main="Rand, MCLUST-ME",xlab="proportion of obs with errors")
  boxplot(rand.mevvv~group,main="Rand, MCLUST",xlab="proportion of obs with errors")
  boxplot(frand.mcme~group,main="Fuzzy Rand, MCLUST-ME",xlab="proportion of obs with errors")
  boxplot(frand.mevvv~group,main="Fuzzy Rand, MCLUST",xlab="proportion of obs with errors")
  par(mfrow=c(1,1))
  
  
  
  #------ Figure 4 ------#
  # Pairwise difference for all p
  rand.diff = rand.mcme - rand.mevvv
  frand.diff = frand.mcme - frand.mevvv
  par(mfrow=c(1,2))
  boxplot(rand.diff~group,main="Rand index",xlab="proportion of obs with errors")
  abline(0,0,lty="dashed")
  boxplot(frand.diff~group,main="Rand index",xlab="proportion of obs with errors")
  abline(0,0,lty="dashed")
  par(mfrow=c(1,1))


  
  #------ Figure 5 ------#
  # Pairwise comparison for p=0.5
  par(mfrow=c(1,2))
  plot(r5.mevvv,r5.mcme-r5.mevvv)
  abline(0,0,lty="dashed")

  plot(fr5.mevvv,fr5.mcme-fr5.mevvv)
  abline(0,0,lty="dashed")
  par(mfrow=c(1,1))
  
  
  
  #------ Figure 6 ------#
  
  
  
  
  #------ Figure 7 ------#
  # Read simulation results
  setwd("/Users/wzhang/Research Project/Simulations for Paper/Example 2/p=0.1/Results")
  load("res.RData")
  rr = out$rand.raw
  r1.mcme = rr[,1]
  r1.mevvv = rr[,4]
  fr1.mcme = rr[,7]
  fr1.mevvv = rr[,8]

  setwd("/Users/wzhang/Research Project/Simulations for Paper/Example 2/p=0.3/Results")
  load("res.RData")
  rr = out$rand.raw
  r3.mcme = rr[,1]
  r3.mevvv = rr[,4]
  fr3.mcme = rr[,7]
  fr3.mevvv = rr[,8]

  setwd("/Users/wzhang/Research Project/Simulations for Paper/Example 2/p=0.5/Results")
  load("res.RData")
  rr = out$rand.raw
  r5.mcme = rr[,1]
  r5.mevvv = rr[,4]
  fr5.mcme = rr[,7]
  fr5.mevvv = rr[,8]

  setwd("/Users/wzhang/Research Project/Simulations for Paper/Example 2/p=0.7/Results")
  load("res.RData")
  rr = out$rand.raw
  r7.mcme = rr[,1]
  r7.mevvv = rr[,4]
  fr7.mcme = rr[,7]
  fr7.mevvv = rr[,8]

  setwd("/Users/wzhang/Research Project/Simulations for Paper/Example 2/p=0.9/Results")
  load("res.RData")
  rr = out$rand.raw
  r9.mcme = rr[,1]
  r9.mevvv = rr[,4]
  fr9.mcme = rr[,7]
  fr9.mevvv = rr[,8]

  rand.mcme = c(r1.mcme,r3.mcme,r5.mcme,r7.mcme,r9.mcme)
  rand.mevvv = c(r1.mevvv,r3.mevvv,r5.mevvv,r7.mevvv,r9.mevvv)
  frand.mcme = c(fr1.mcme,fr3.mcme,fr5.mcme,fr7.mcme,fr9.mcme)
  frand.mevvv = c(fr1.mevvv,fr3.mevvv,fr5.mevvv,fr7.mevvv,fr9.mevvv)
  group = c(rep(0.1,100),rep(0.3,100),rep(0.5,100),rep(0.7,100),rep(0.9,100))

  # Individual Rand/Fuzzy Rand indices
  par(mfrow=c(2,2))
  boxplot(rand.mcme~group,main="Rand, MCLUST-ME",xlab="proportion of obs with errors")
  boxplot(rand.mevvv~group,main="Rand, MCLUST",xlab="proportion of obs with errors")
  boxplot(frand.mcme~group,main="Fuzzy Rand, MCLUST-ME",xlab="proportion of obs with errors")
  boxplot(frand.mevvv~group,main="Fuzzy Rand, MCLUST",xlab="proportion of obs with errors")
  par(mfrow=c(1,1))


  
  #------ Figure 8 ------#
  # Pairwise difference for all p
  rand.diff = rand.mcme - rand.mevvv
  frand.diff = frand.mcme - frand.mevvv

  par(mfrow=c(2,2))
  boxplot(rand.diff~group,main="Rand index",xlab="proportion of obs with errors")
  abline(0,0,lty="dashed")
  boxplot(frand.diff~group,main="Fuzzy Rand index",xlab="proportion of obs with errors")
  abline(0,0,lty="dashed")

  # Pairwise comparison for p=0.5
  plot(r5.mevvv,r5.mcme-r5.mevvv,xlab="MCLUST",ylab="MCLUST-ME--MCLUST",main="Rand index, p = 0.5")
  abline(0,0,lty="dashed")

  plot(fr5.mevvv,fr5.mcme-fr5.mevvv,xlab="MCLUST",ylab="MCLUST-ME--MCLUST",main="Fuzzy Rand index, p = 0.5")
  abline(0,0,lty="dashed")
  par(mfrow=c(1,1))

  
  
  #------ Figure 9 ------#
  # Plot for paper:
  layout(matrix(1:6, 3, 2, byrow = TRUE))
  setwd("/Users/wzhang/Research Project/Simulations for Paper/BIC simulation/Results")
  load("bic_well_separated.RData")
  plot.wellsep(101,bicres$res5)

  load("bic2_79.RData")
  plot.2close(79,bic2_79)
  
  load("bic3_51.RData")
  plot.3close(51,bic3_51)
  par(mfrow=c(1,1))
  
  
  
  #------ Figure 10 ------#
  ### Raw data scatterplot matrix:
  pairs(obs,labels=c("10min","1h","3h","6h","12h"),pch=16,cex=0.6)
  
  
}  
  
  