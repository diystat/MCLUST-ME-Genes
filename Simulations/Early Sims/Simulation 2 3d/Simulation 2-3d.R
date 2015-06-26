

          #---------------------Cluster Simulation 2-3d --------------------------#

          ### 3-dimensional clusters
          ### One cluster has measurement error, the other does not.


library(mclust)
library(gdata)
library(MASS)
library(phyclust)
library(lattice)

col.blue = rgb(30,144,255,max=255)


sim2.3d = function(seed,n,lambda,imgname){

  n = n
  set.seed(seed)
  
## Simulate two clusters far away from each other.
  mu1 = c(-3,0,0)
  mu2 = c(3,0,0)
  sig = matrix(c(4,0,0,0,4,0,0,0,4),nrow=3)
  
  I = matrix(c(1,0,0,0,1,0,0,0,1),nrow=3)
  E = (lambda^2)*I

  # Sample size:
  n1 = n2 = n
  s = n1 + n2

  # Generate samples:
  dat1 = mvrnorm(n1, mu1, (sig+E))
  dat2 = mvrnorm(n2, mu2, sig)
  
  # Combine two groups:
  dat = rbind(dat1, dat2)
  
 
  
  
  #-------------------------------------------------------------#




  # Construct initial classification. Use true membership:
  z.ini = matrix(c(rep(c(1,0),n1),rep(c(0,1),n2)), byrow=T, nrow=s)

  # Construct the true error matrices:
  errmat = array(0,dim=c(3,3,s))
  
  for(i in 1:n1){
    errmat[,,i] = E
  }


  
  #-------------------------------------------------------------#



  # Use meVVV
  res.me = meVVV(dat,z.ini)
  z1 = res.me$z

  t = nrow(z1)
  G = ncol(z1)
  
  for(i in 1:t){
    rowmax = max(z1[i,])
    for(k in 1:G){
      z1[i,k] = ifelse(z1[i,k]==rowmax,1,0)
    }  
  }
    
 
  
  
    #-------------------------------------------------------------#

  
  
  
  ## Use MCME:
    
  # Perform our clustering method:
  my.result = mcmeVVV(dat, z.ini, errmat)
  z2 = my.result$z
  
  t = nrow(z2)
  G = ncol(z2)
  
  for(i in 1:t){
    rowmax = max(z2[i,])
    for(k in 1:G){
      z2[i,k] = ifelse(z2[i,k]==rowmax,1,0)
    }  
  }
    
 
  
  
  #-----------------------------------------------------------------------#


  
  ## Plot scatterplot matrix for both results:
  labels = c(z.ini[,1],z2[,1],z1[,1])
  temp = cbind(rbind(dat,dat,dat),labels)
  dat.comb = data.frame(x1 = temp[,1],x2 = temp[,2],x3=temp[,3],label=temp[,4])
  group = c(rep("Generated Clusters",s),rep("mcme",s),rep("meVVV",s))

  temp = trellis.par.get()
  supsym <- trellis.par.get("superpose.symbol")
    supsym$col <- c(col.blue,"red")
    supsym$pch <- c(19,0)
  trellis.par.set("superpose.symbol",supsym)
    	  
  a = splom(~dat.comb[,1:3]|group, data=dat.comb, groups=dat.comb[,4])
  
  png(imgname,1200,400)
  plot(a)
  dev.off()
  
  trellis.par.set(temp)


  #-----------------------------------------------------------------------#


  
  ## Compute Rand index
  z.ini = z.ini + 1
  z1 = z1 + 1
  z2 = z2 + 1
  rand.mevv = RRand(z.ini,z1)
  rand.mcme = RRand(z.ini,z2)
  
  

  #-----------------------------------------------------------------------#



  ## Misclassification rate
  tab1 = table(z.ini,z1)
  tab2 = table(z.ini,z2)
  mis.mevv = 1-sum(diag(tab1))/sum(tab1)
  mis.mcme = 1-sum(diag(tab2))/sum(tab2)
  z.ini = z.ini - 1
  z1 = z1 - 1
  z2 = z2 - 1
  
  out = list(seed=seed,n=n,rand.mevv=rand.mevv$Rand,rand.mcme=rand.mcme$Rand,
    mis.mevv=mis.mevv,mis.mcme=mis.mcme)
  return(out)
  
  
}


  #-----------------------------------------------------------------------#



## Vector of random seeds:
seedvec = c(91,92,93,94,95)



## lambda = 1:
imgname50 = sapply(seedvec,paste,"-50.png",sep="")
res50 = list()
for(i in 1:5){
  res50[[i]] = sim2.3d(seedvec[i],n=50,lambda=1,imgname50[i])
}




imgname100 = sapply(seedvec,paste,"-100.png",sep="")
res100 = list()
for(i in 1:5){
  res100[[i]] = sim2.3d(seedvec[i],n=100,lambda=1,imgname100[i])
}



imgname200 = sapply(seedvec,paste,"-200.png",sep="")
res200 = list()
for(i in 1:5){
  res200[[i]] = sim2.3d(seedvec[i],n=200,lambda=1,imgname200[i])
}



imgname500 = sapply(seedvec,paste,"-500.png",sep="")
res500 = list()
for(i in 1:5){
  res500[[i]] = sim2.3d(seedvec[i],n=500,lambda=1,imgname500[i])
}



imgname800 = sapply(seedvec,paste,"-800.png",sep="")
res800 = list()
for(i in 1:5){
  res800[[i]] = sim2.3d(seedvec[i],n=800,lambda=1,imgname800[i])
}




out1 = list(res50=res50,res100=res100,res200=res200,
  res500=res500,res800=res800)
save(out1,file="sim2-3d-1")








## lambda = 0.1:
imgname50 = sapply(seedvec,paste,"-50.1.png",sep="")
res50.1 = list()
for(i in 1:5){
  res50.1[[i]] = sim2.3d(seedvec[i],n=50,lambda=0.1,imgname50[i])
}




imgname100 = sapply(seedvec,paste,"-100.1.png",sep="")
res100.1 = list()
for(i in 1:5){
  res100.1[[i]] = sim2.3d(seedvec[i],n=100,lambda=0.1,imgname100[i])
}



imgname200 = sapply(seedvec,paste,"-200.1.png",sep="")
res200.1 = list()
for(i in 1:5){
  res200.1[[i]] = sim2.3d(seedvec[i],n=200,lambda=0.1,imgname200[i])
}



imgname500 = sapply(seedvec,paste,"-500.1.png",sep="")
res500.1 = list()
for(i in 1:5){
  res500.1[[i]] = sim2.3d(seedvec[i],n=500,lambda=0.1,imgname500[i])
}



imgname800 = sapply(seedvec,paste,"-800.1.png",sep="")
res800.1 = list()
for(i in 1:5){
  res800.1[[i]] = sim2.3d(seedvec[i],n=800,lambda=0.1,imgname800[i])
}



out2 = list(res50=res50.1,res100=res100.1,res200=res200.1,
  res500=res500.1,res800=res800.1)
save(out2,file="sim2-3d-2")





