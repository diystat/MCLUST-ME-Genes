

## Simulate two clusters that are clearly separated but close to each other.
  mu1 = c(-3,0)
  mu2 = c(3,0)
  # 3.2 is chosen as square root of chi(0.01) with df=2, so that 99% of data are expected
  # to fall into the region
  sig = matrix(c(1,0,0,1),nrow=2)

  n = 60

  dat1 = mvrnorm(n, mu1, sig)
  dat2 = mvrnorm(n, mu2, sig)

  plot(dat1, xlim = c(-7,7), ylim = c(-5,5), col="blue")
  points(dat2, col="red")
  # Observe that the two clusters are separated but close.
  
  
  # Now we truncate each group and keep the observations with x-coord <=-3 and >=3:
  dat1.tr = dat1[which(dat1[,1]<=-3),]
  dat2.tr = dat2[which(dat2[,1]>=3),]
  
  plot(dat1.tr, xlim=c(-7,7), ylim=c(-5,5), col="blue")
  points(dat2.tr, col="red")
  abline(v=3, lty="dashed")
  abline(v=-3, lty="dashed")

  
  # Define an error matrix:
  err = matrix(c(4,0,0,4),nrow=2)
  
  sig.new = sig + err
  dat1.new = mvrnorm(n, mu1, sig.new)
  dat2.new = mvrnorm(n, mu2, sig.new)
  
  dat1.new.tr = dat1.new[which(dat1.new[,1]>-3),]
  dat2.new.tr = dat2.new[which(dat2.new[,1]<3),]
  
  group1 = rbind(dat1.tr, dat1.new.tr)
  group2 = rbind(dat2.tr, dat2.new.tr)
  
  plot(group1, xlim=c(-7,7), ylim=c(-7,7), pch=0)
  points(group2, pch=2)
  abline(v=3, lty="dashed")
  abline(v=-3, lty="dashed")
  # The combined dataset has points in the middle are commingled.
  
  
  #-----------------------------------------------------------------------------#
  
  
  ## Automatic clustering by MCLUST
  # Combine two groups:
  dat = rbind(group1, group2)
  
  # Run Mclust on the whole data set:
  library(mclust)
  result = Mclust(dat)
  result$modelName
  result$G
  
  result$classification
  datSum = summary(result)
  
  mclust2Dplot(dat, classification=datSum$classification,
               parameters=datSum$parameters)
  abline(v=-3,lty="dashed")
  abline(v=3, lty="dashed")
  # Mclust recognizes three groups
  
  
  
  #-------------------------------------------------------------#
  
  
  ## First, try "meclust"
  s = nrow(dat)
  
  # Construct initial classification:
  z.ini = matrix(c(rep(c(1,0),nrow(dat1.tr)),rep(c(0,1),nrow(dat2.tr)),
    rep(c(1,0),nrow(dat1.new.tr)),rep(c(0,1),nrow(dat2.new.tr))), byrow=T, nrow=s)
  
  res0 = meVVV(dat, z.ini)
  z0 = res0$z
  
  t = nrow(z0)
  G = ncol(z0)
  
  for(i in 1:t){
    rowmax = max(z0[i,])
    for(k in 1:G){
      z0[i,k] = ifelse(z0[i,k]==rowmax,1,0)
    }  
  }
  
  dat.new = cbind(dat,z0)
  
  newgroup1 = dat.new[which(dat.new[,3]==1),][,1:2]
  newgroup2 = dat.new[which(dat.new[,3]==0),][,1:2]
  
  plot(newgroup1,col="blue",xlim=c(-7,7),ylim=c(-7,7))
  points(newgroup2,col="red")
  abline(v=-3,lty="dashed")
  abline(v=3,lty="dashed")
  
  
  ## Then use our method, assuming all errors are zero:
  zeroerr = array(0,dim=c(2,2,s))
  myres = ME.VVV.err(dat, z.ini, zeroerr, "identical")
  names(myres)
  gr = myres$z
  
  t = nrow(gr)
  G = ncol(gr)
  
  for(i in 1:t){
    rowmax = max(gr[i,])
    for(k in 1:G){
      gr[i,k] = ifelse(gr[i,k]==rowmax,1,0)
    }  
  }
  
  dat.new = cbind(dat,gr)
  
  newgroup1 = dat.new[which(dat.new[,3]==1),][,1:2]
  newgroup2 = dat.new[which(dat.new[,3]==0),][,1:2]
  
  plot(newgroup1,col="blue",xlim=c(-7,7),ylim=c(-7,7))
  points(newgroup2,col="red")
  abline(v=-3,lty="dashed")
  abline(v=3,lty="dashed")
  # Same results, and neither is a good one.
  
  
  #--------------------------------------------------------------------------#
  
  
  ## Use MEVVV:
  # Construct the true error matrices:
  temp = array(0,dim=c(2,3,s))
  
  for(i in 1:s){
    temp[,1,i] = dat[i,]
  }
  
  for(i in 1:s){
    if(temp[1,1,i]>-3 && temp[1,1,i]<3){
      temp[1,2,i] = temp[2,3,i] = 4
    }
  }
  
  errmat = array(0,dim=c(2,2,s))
  for(i in 1:s){
    errmat[,,i] = temp[,2:3,i]
  }
  
  
  # Perform our clustering method:
  my.result = ME.VVV.err(dat, z.ini, errmat, "none")
  names(my.result)
  gr = my.result$z
  
  t = nrow(gr)
  G = ncol(gr)
  
  for(i in 1:t){
    rowmax = max(gr[i,])
    for(k in 1:G){
      gr[i,k] = ifelse(gr[i,k]==rowmax,1,0)
    }  
  }
  
  
  ## Visualizing the new clusters:
  # Shapes indicate true grouping; colors indicate grouping by our method.
  n1 = nrow(group1)
  group1.new = cbind(group1,gr[1:n1,])
  group2.new = cbind(group2,gr[(n1+1):s,])
  
  plot(group1, xlim=c(-7,7), ylim=c(-8,8), pch=0, col=ifelse(group1.new[,3]==1,"blue","red"),xlab="x1", ylab="x2")
  points(group2, pch=2, col=ifelse(group2.new[,3]==1,"blue","red"))
  abline(v=3, lty="dashed")
  abline(v=-3, lty="dashed")
  
  