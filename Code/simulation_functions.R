
 
#######################################################################
###------------------- Functions for Simulation --------------------###
#######################################################################



    ###########################################################
    ################ Two-Cluster Simulation ###################
    ###########################################################
    
### Simulation 1 wrapper    
sim1 = function(p){
  N = 200
  tau = 0.5
  mu1 = c(0,0)
  mu2 = c(8,0)
  sig1 = matrix(c(64,0,0,64),nrow=2)
  sig2 = matrix(c(16,0,0,16),nrow=2)
  k = 36
  nseed = 100

  # Check if singularity occurs for any seed
  check.seed(N,tau,mu1,mu2,sig1,sig2,k,p)

  # Run the simulation
  out = sim.driver(N,tau,mu1,mu2,sig1,sig2,k,p,nseed)
  return(out)
}        
    
    

### Simulation 2 wrapper    
sim2 = function(p){
  N = 200
  tau = 0.5
  mu1 = c(0,0)
  mu2 = c(8,0)
  sig1 = matrix(c(64,0,0,64),nrow=2)
  sig2 = matrix(c(16,0,0,16),nrow=2)
  k = 9
  nseed = 100

  # Check if singularity occurs for any seed
  check.seed(N,tau,mu1,mu2,sig1,sig2,k,p)

  # Run the simulation
  out = sim.driver(N,tau,mu1,mu2,sig1,sig2,k,p,nseed)
  return(out)
}        
    



## Examine clustering uncertainties
plot.uncr = function(seed,out){
  k = which(out$seed==seed)
  
  # Extract clustering result for seed="seed"
  i = (k-1)*7+1
  res.mcme = out$sim.result[[i]]
  res.mevvv = out$sim.result[[(i+1)]]
  z.ini = out$sim.result[[(i+2)]]
  rand.samples = out$sim.result[[(i+3)]]
  errmat = out$sim.result[[(i+4)]]
  index = out$sim.result[[(i+5)]]
  k = out$sim.result[[(i+6)]]
  simrun = list(res.mcme=res.mcme,res.mevvv=res.mevvv,
    errmat=errmat,z.ini=z.ini,rand.samples=rand.samples,index=index,k=k)

  # Obtain classification uncertainty vector
  unc.mcme = res.mcme$uncertainty
  unc.mclust = numeric(200)
  for(i in 1:200){
    unc.mclust[i] = 1-max(res.mevvv$z[i,])
  }

  point.size = 0.3+unc.mcme*4
  point.size.mclust = 0.3+unc.mclust*4

  plot.boundary.new(simrun,point.size,point.size.mclust)
}





### implementation of fuzzy rand index described in Campello(2006)
fuzzyrand = function(R,Q,tnorm="min"){
  
  # "tnorm" denotes the triangular norm in fuzzy set theory
  # User can choose to use min() or product as tnorm.
  # default is min().
  
  # Note: dimension of R and Q is GxN, so transpose the membership matrix
  # before using as input
  
  if(tnorm=="min"){
    tnorm = function(a,b){
      return(min(a,b))
    }
  } else if(tnorm=="product"){
    tnorm = function(a,b){
      return(a*b)
    }
  }
  
  k = dim(R)[1]
  N = dim(R)[2]
  v = dim(Q)[1]
  
  ## construct fuzzy set of data pairs:
  V = X = Y = Z = matrix(0,N,N)
  for(j1 in 1:N){
    for(j2 in 1:N){
      
      temp1 = numeric(k)
      temp2 = numeric(v)
      
      for(i in 1:k){
        temp1[i] = tnorm(R[i,j1],R[i,j2])
      }
      
      for(j in 1:k){
        temp2[j] = tnorm(Q[j,j1],Q[j,j2])
      }
      
      V[j1,j2] = max(temp1)
      Y[j1,j2] = max(temp2)
      
      temp3 = matrix(0,k,k)
      temp4 = matrix(0,v,v)
      for(i in 1:k){
        for(j in 1:k){
          temp3[i,j] = tnorm(R[i,j1],R[j,j2])
        }
      }
      
      for(i in 1:v){
        for(j in 1:v){
          temp4[i,j] = tnorm(Q[i,j1],Q[j,j2])
        }
      }
      
      diag(temp3) = -1
      diag(temp4) = -1
      
      X[j1,j2] = max(temp3)
      Z[j1,j2] = max(temp4)
    }
  }
  
  ## calculate a,b,c,d:
  require(gdata)
  A = B = C = D = matrix(0,N,N)
  for(j1 in 1:N){
    for(j2 in 1:N){
      A[j1,j2] = tnorm(V[j1,j2],Y[j1,j2])
      B[j1,j2] = tnorm(V[j1,j2],Z[j1,j2])
      C[j1,j2] = tnorm(X[j1,j2],Y[j1,j2])
      D[j1,j2] = tnorm(X[j1,j2],Z[j1,j2])
    }
  }
  a = sum(upperTriangle(A))
  b = sum(upperTriangle(B))
  c = sum(upperTriangle(C))
  d = sum(upperTriangle(D))
  
  ## calculate fuzzy rand index:
  out = (a+d)/(a+b+c+d)
  return(out)
}



### implementation of weighted fuzzy rand index described in Campello(2006)
wtdfuzzyrand = function(alpha,R,Q,tnorm="min",err){
  
  # "tnorm" denotes the triangular norm in fuzzy set theory
  # User can choose to use min() or product as tnorm.
  # default is min().
  
  if(tnorm=="min"){
    tnorm = function(a,b){
      return(min(a,b))
    }
  } else if(tnorm=="product"){
    tnorm = function(a,b){
      return(a*b)
    }
  }
  
  k = dim(R)[1]
  N = dim(R)[2]
  v = dim(Q)[1]
  
  tv = numeric(N)
  for(i in 1:N){
    tv[i] = det(err[,,i])
  }
  
  for(i in 1:N){
    Q[,i] = Q[,i]/(1+alpha*tv[i])
    #R[,i] = R[,i]/(1+alpha*tv[i])
  }
 
  ## construct fuzzy set of data pairs:
  V = X = Y = Z = matrix(0,N,N)
  for(j1 in 1:N){
    for(j2 in 1:N){
      
      temp1 = numeric(k)
      temp2 = numeric(v)
      
      for(i in 1:k){
        temp1[i] = tnorm(R[i,j1],R[i,j2])
      }
      
      for(j in 1:k){
        temp2[j] = tnorm(Q[j,j1],R[j,j2])
      }
      
      V[j1,j2] = max(temp1)
      Y[j1,j2] = max(temp2)
      
      temp3 = matrix(0,k,k)
      temp4 = matrix(0,v,v)
      for(i in 1:k){
        for(j in 1:k){
          temp3[i,j] = tnorm(R[i,j1],R[j,j2])
        }
      }
      
      for(i in 1:v){
        for(j in 1:v){
          temp4[i,j] = tnorm(Q[i,j1],Q[j,j2])
        }
      }
      
      diag(temp3) = -1
      diag(temp4) = -1
      
      X[j1,j2] = max(temp3)
      Z[j1,j2] = max(temp4)
    }
  }
  
  ## calculate a,b,c,d:
  require(gdata)
  A = B = C = D = matrix(0,N,N)
  for(j1 in 1:N){
    for(j2 in 1:N){
      A[j1,j2] = tnorm(V[j1,j2],Y[j1,j2])
      B[j1,j2] = tnorm(V[j1,j2],Z[j1,j2])
      C[j1,j2] = tnorm(X[j1,j2],Y[j1,j2])
      D[j1,j2] = tnorm(X[j1,j2],Z[j1,j2])
    }
  }
  
  a = sum(upperTriangle(A))
  b = sum(upperTriangle(B))
  c = sum(upperTriangle(C))
  d = sum(upperTriangle(D))
  
  ## calculate fuzzy rand index:
  out = (a+d)/(a+b+c+d)
  return(out)
}


## Function to convert membership matrix into vector
## of classification uncertainties:
unc = function(z){
  ## Input argument: z---membership matrix
  n = nrow(z)
  uncertainty = numeric() # records classification uncertainty of each obs.
  for(i in 1:n){
    rowmax = max(z[i,])
    uncertainty[i] = 1-rowmax
  }
  return(uncertainty)
}


# Calculate rand indices
get.rand = function(simrun){
  
  z.ini = simrun$z.ini
  z.mcme = simrun$res.mcme$z
  z.mevvv = simrun$res.mevvv$z
  index = simrun$index
  
  N = nrow(z.mcme)
  G = ncol(z.mcme)
  
  # Discretize membership matrix:
  z1 = z2 = matrix(0,N,G)
  for(i in 1:N){
    rm.mcme = max(z.mcme[i,])
    rm.mevvv = max(z.mevvv[i,])
    for(k in 1:G){
      z1[i,k] = ifelse(z.mcme[i,k]==rm.mcme,1,0)
      z2[i,k] = ifelse(z.mevvv[i,k]==rm.mevvv,1,0)
    }  
  }
  
  z.ini1 = z.ini[,1] + 1
  z1 = z1[,1] + 1
  z2 = z2[,1] + 1
  
  ## Calculate Rand index:
  mcme = round(RRand(z.ini1,z1)$Rand,4)
  mcme.me = round(RRand(z.ini1[which(index==1)],z1[which(index==1)])$Rand,4)
  mcme.nn = round(RRand(z.ini1[which(index==0)],z1[which(index==0)])$Rand,4)
  mevvv = round(RRand(z.ini1,z2)$Rand,4)
  mevvv.me = round(RRand(z.ini1[which(index==1)],z2[which(index==1)])$Rand,4)
  mevvv.nn = round(RRand(z.ini1[which(index==0)],z2[which(index==0)])$Rand,4)
  mcme.f = round(fuzzyrand(t(z.ini),t(z.mcme)),4)
  mevvv.f = round(fuzzyrand(t(z.ini),t(z.mevvv)),4)
  
  out = c(mcme,mcme.me,mcme.nn,mevvv,mevvv.me,mevvv.nn,mcme.f,mevvv.f)
  return(out)
}

# Driver function for simulation 1 and 2
sim.driver = function(N,tau,mu1,mu2,sig1,sig2,k,p,nseed){

  simres = list()
  randres = matrix(0,nrow=nseed,ncol=8)
  
  v = seq(1,2*nseed,by=1)
  seedvec = sample(v,size=nseed)
  
  for(i in 1:nseed){
    
    print(paste("seed=",seedvec[i],sep=""))
    
    # obtain and store clustering results    
    data.par = sim.par(N,tau,mu1,mu2,sig1,sig2,k,p)
    set.seed(seedvec[i])
    simdata = sim.data(data.par)
    simrun = sim.run(simdata)
    simres = c(simres,simrun)
    
    # obtain and store rand indices
    rindex = get.rand(simrun)
    randres[i,] = rindex
    
    # create and store boundary plots
    file = paste("bdry",seedvec[i],".png",sep="")
    png(filename=file,1500,600)
    plot.boundary(simrun,seedvec[i])
    dev.off()
  }
  
  # calculate average rand index
  rand.avg = colMeans(randres)
  
  out = list(sim.result=simres, rand.raw=randres, rand.avg=rand.avg, seed=seedvec)
  save(out,file="res.RData")
}


## Specifying model parameters used to generate data
sim.par = function(N,tau,mu1,mu2,sig1,sig2,k,p){
  out = list(N=N,tau=tau,mu1=mu1,mu2=mu2,sig1=sig1,sig2=sig2,k=k,p=p)
  return(out)
}

## Generate data using given parameters
sim.data = function(simpar){  
  N = simpar$N; tau = simpar$tau
  mu1 = simpar$mu1; mu2 = simpar$mu2
  sig1 = simpar$sig1; sig2 = simpar$sig2
  k = simpar$k; p = simpar$p  
  E = k*matrix(c(1,0,0,1),nrow=2)
    
  ## Randomly assign constant measurement error to 100p% of observations:
  index = rbinom(N,1,p)  
  errmat = array(0,dim=c(2,2,N))
  for(i in 1:N){
    errmat[,,i] = index[i]*E
  }
 
  ## Generate sample from mixture distribution:
  U = runif(N)
  z.ini = matrix(0,N,2)
  rand.samples = matrix(0,N,2)
  for(i in 1:N){
    if(U[i]<tau){
      rand.samples[i,] = MASS::mvrnorm(1,mu1,(sig1+errmat[,,i]))
      z.ini[i,] = c(1,0)
    } else{
      rand.samples[i,] = MASS::mvrnorm(1,mu2,(sig2+errmat[,,i]))
      z.ini[i,] = c(0,1)
    }
  }

  out = list(data=rand.samples, z.ini=z.ini, err=errmat, index=index, k=k)
  return(out)
}

## Obtain results from MCME and meVVV
sim.run = function(simdata){
  rand.samples = simdata$data
  z.ini = simdata$z.ini
  errmat = simdata$err
  index = simdata$index
  k = simdata$k
  
  ## Run MCME:
  res.mcme = mcmeVVV(rand.samples, z.ini, errmat)
  
  ## Run mevvv:
  res.mevvv = mclust::meVVV(rand.samples,z.ini)
  
  out = list(res.mcme=res.mcme, res.mevvv=res.mevvv, z.ini=z.ini, rand.samples=rand.samples, errmat=errmat, index=index, k=k)
  return(out)
}


## Check if mcmeVVV has singularity issue for any of the random seeds
check.seed = function(N,tau,mu1,mu2,sig1,sig2,k,p){
  seedvec = seq(1,200,1)
  for(i in 1:200){
    data.par = sim.par(N,tau,mu1,mu2,sig1,sig2,k,p)
    print(seedvec[i])
    set.seed(seedvec[i])
    simdata = sim.data(data.par)
  
    rand.samples = simdata$data
    z.ini = simdata$z.ini
  
    res = meVVV(rand.samples,z.ini)
  }
}


## Obtain coordinates for two boundaries
get.boundary = function(y,simrun,k){
  E = k*matrix(c(1,0,0,1),nrow=2)
  param = simrun$res.mcme$param
  errmat = simrun$errmat
  tau1 = param$tauhat[1]
  tau2 = param$tauhat[2]
  sig1 = param$sigmahat[,,1]
  sig2 = param$sigmahat[,,2]
  mu1 = as.numeric(param$muhat[,1])
  mu2 = as.numeric(param$muhat[,2])
  
  part1 = tau1/sqrt(det(sig1+E))*exp((-1/2)*t(c(y[1],y[2])-mu1)%*%solve(sig1+E)%*%(c(y[1],y[2])-mu1))
  part2 = tau2/sqrt(det(sig2+E))*exp((-1/2)*t(c(y[1],y[2])-mu2)%*%solve(sig2+E)%*%(c(y[1],y[2])-mu2))
  f = part1 - part2
  return(f)
}


## Obtain boundary for meVVV:
get.boundary.mevvv = function(y,simrun){
  param = simrun$res.mevvv$param
  tau1 = param$pro[1]
  tau2 = param$pro[2]
  sig1 = param$var$sigma[,,1]
  sig2 = param$var$sigma[,,2]
  mu1 = as.numeric(param$mean[,1])
  mu2 = as.numeric(param$mean[,2])
  
  part1 = tau1/sqrt(det(sig1))*exp((-1/2)*t(c(y[1],y[2])-mu1)%*%solve(sig1)%*%(c(y[1],y[2])-mu1))
  part2 = tau2/sqrt(det(sig2))*exp((-1/2)*t(c(y[1],y[2])-mu2)%*%solve(sig2)%*%(c(y[1],y[2])-mu2))
  f = part1 - part2
  return(f)
}


## Plot boundaries for MCME and meVVV
plot.boundary = function(simrun,seed){
  ## obtain sim results:
  res.mcme = simrun$res.mcme
  z.ini = simrun$z.ini
  rand.samples = simrun$rand.samples
  
  ## specify graphing parameters:
  col.blue = rgb(30,144,255,max=255)  
  xlow = min(rand.samples[,1])-0.5
  xup = max(rand.samples[,1])+0.5
  ylow = min(rand.samples[,2])-0.2
  yup = max(rand.samples[,2])+0.2

  ## evaluate boundary function over grid:
  x = seq(xlow, xup, len=250)
  y = seq(ylow, yup, len=250)
  z <- outer(x, y,
  FUN=function(x, y, ...){
  apply(cbind(x,y), 1, get.boundary, ...)
  }, simrun=simrun, k=0)
  
  kk = simrun$k
  zz <- outer(x, y,
  FUN=function(x, y, ...){
  apply(cbind(x,y), 1, get.boundary, ...)
  }, simrun=simrun, k=kk)
  
  ze <- outer(x, y,
  FUN=function(x, y, ...){
  apply(cbind(x,y), 1, get.boundary.mevvv, ...)
  }, simrun=simrun)
  
  ## add clustering result to contour plot
  z.mcme = res.mcme$z
  t = nrow(z.mcme)
  G = ncol(z.mcme)

  # Discretize membership matrix:
  z1 = (z.mcme[,1]>z.mcme[,2])

  # distinguish error-free and erroneous obs:
  index = simrun$index

  rs.new = cbind(rand.samples,z1,index,z.ini[,1])  
  plot(rand.samples,main="MCLUST-ME",cex.main=1,pch=ifelse(rs.new[,3]==1,22,24),
    bg=ifelse(rs.new[,4]==1,"black","white"), xlab="",ylab="",xlim=c(xlow,xup),ylim=c(ylow,yup),cex=0.7)
  
  # boundary for error-free obs:
  contour(x,y,z,level=0,add=TRUE,drawlabels=FALSE,lty="dashed",lwd=2)
  # boundary for erroneous obs:
  contour(x,y,zz,level=0,add=TRUE,drawlabels=FALSE,lty="dotted",lwd=2)
  # boundary for mevvv:
  contour(x,y,ze,level=0,add=TRUE,drawlabels=FALSE,lwd=2)
  
  z.mevvv = simrun$res.mevvv$z
  z2 = (z.mevvv[,1]>z.mevvv[,2])
  rs.new1 = cbind(rand.samples,z2,index)  
  plot(rand.samples,pch=ifelse(rs.new1[,3]==1,22,24),main="MCLUST",cex.main=1,
    bg=ifelse(rs.new1[,4]==1,"black","white"),xlab="",ylab="",xlim=c(xlow,xup),ylim=c(ylow,yup),cex=0.7)
  contour(x,y,ze,level=0,add=TRUE,drawlabels=FALSE,lwd=2)
} 


## Plot boundaries for MCME and meVVV with pointsize reflecting clustering uncertainty
plot.boundary.new = function(simrun,pointsize,pointsize2){
  ## obtain sim results:
  res.mcme = simrun$res.mcme
  z.ini = simrun$z.ini
  rand.samples = simrun$rand.samples
  
  ## specify graphing parameters:
  col.blue = rgb(30,144,255,max=255)  
  xlow = min(rand.samples[,1])-0.5
  xup = max(rand.samples[,1])+0.5
  ylow = min(rand.samples[,2])-0.2
  yup = max(rand.samples[,2])+0.2

  ## evaluate boundary function over grid:
  x = seq(xlow, xup, len=250)
  y = seq(ylow, yup, len=250)
  z <- outer(x, y,
  FUN=function(x, y, ...){
  apply(cbind(x,y), 1, get.boundary, ...)
  }, simrun=simrun, k=0)
  
  kk = simrun$k
  zz <- outer(x, y,
  FUN=function(x, y, ...){
  apply(cbind(x,y), 1, get.boundary, ...)
  }, simrun=simrun, k=kk)
  
  ze <- outer(x, y,
  FUN=function(x, y, ...){
  apply(cbind(x,y), 1, get.boundary.mevvv, ...)
  }, simrun=simrun)
  

  ## add clustering result to contour plot
  z.mcme = res.mcme$z
  t = nrow(z.mcme)
  G = ncol(z.mcme)

  # Discretize membership matrix:
  z1 = (z.mcme[,1]>z.mcme[,2])

  # distinguish error-free and erroneous obs:
  index = simrun$index

  rs.new = cbind(rand.samples,z1,index,z.ini[,1])  
  plot(rand.samples,main="MCLUST-ME",cex.main=1,pch=ifelse(rs.new[,3]==1,22,24),
    bg=ifelse(rs.new[,4]==1,"black","white"), xlab="",ylab="",xlim=c(xlow,xup),ylim=c(ylow,yup),cex=pointsize)
  
  # boundary for error-free obs:
  contour(x,y,z,level=0,add=TRUE,drawlabels=FALSE,lty="dashed",lwd=2)
  # boundary for erroneous obs:
  contour(x,y,zz,level=0,add=TRUE,drawlabels=FALSE,lty="dotted",lwd=2)
  # boundary for mevvv:
  contour(x,y,ze,level=0,add=TRUE,drawlabels=FALSE,lwd=2)
  
  
  z.mevvv = simrun$res.mevvv$z
  z2 = (z.mevvv[,1]>z.mevvv[,2])
  rs.new1 = cbind(rand.samples,z2,index)  
  plot(rand.samples,pch=ifelse(rs.new1[,3]==1,22,24),main="MCLUST",cex.main=1,
    bg=ifelse(rs.new1[,4]==1,"black","white"),xlab="",ylab="",xlim=c(xlow,xup),ylim=c(ylow,yup),cex=pointsize2)
  contour(x,y,ze,level=0,add=TRUE,drawlabels=FALSE,lwd=2)
  #par(mfrow=c(1,1))
  
} 









      ########################################################
      ################### BIC Simulation #####################
      ########################################################

## Simulate 3-component data
## Specifying model parameters used to generate data
sim.par3 = function(N,tau1,tau2,mu1,mu2,mu3,sig1,sig2,sig3,k,p){
  out = list(N=N,tau1=tau1,tau2=tau2,mu1=mu1,mu2=mu2,mu3=mu3,
    sig1=sig1,sig2=sig2,sig3=sig3,k=k,p=p)
  return(out)
}

## Generate data using given parameters
sim.data3 = function(simpar){  
  N = simpar$N; tau1 = simpar$tau1; tau2 = simpar$tau2
  mu1 = simpar$mu1; mu2 = simpar$mu2; mu3 = simpar$mu3
  sig1 = simpar$sig1; sig2 = simpar$sig2; sig3 = simpar$sig3
  k = simpar$k; p = simpar$p  
  E = k*matrix(c(1,0,0,1),nrow=2)
    
  ## Randomly assign constant measurement error to 100p% of observations:
  index = rbinom(N,1,p)  
  errmat = array(0,dim=c(2,2,N))
  for(i in 1:N){
    errmat[,,i] = index[i]*E
  }
 
  ## Generate sample from mixture distribution:
  U = runif(N)
  z.ini = matrix(0,N,2)
  rand.samples = matrix(0,N,2)
  for(i in 1:N){
    if(U[i]<tau1){
      rand.samples[i,] = MASS::mvrnorm(1,mu1,(sig1+errmat[,,i]))
      z.ini[i,] = c(1,0)
    } else if(U[i]<(tau1+tau2)){
      rand.samples[i,] = MASS::mvrnorm(1,mu2,(sig2+errmat[,,i]))
      z.ini[i,] = c(0,1)
    } else{
      rand.samples[i,] = MASS::mvrnorm(1,mu3,(sig3+errmat[,,i]))
    }
  }

  out = list(data=rand.samples, z.ini=z.ini, err=errmat, index=index, k=k)
  return(out)
}


## Implement selection of number of clusters with BIC ##
nc.select = function(data, err, nc=1:5, d=1, itmax=Inf, lb=1e-3){
  ## "maxClnum" lets user specify the maximum number of clusters
  ## in consideration.
  
  G = Gmin = min(nc)
  Gmax = max(nc)
  bic.vec = numeric(length(nc))
  N = nrow(data)
  hcTree = mclust::hc(data)
  k = 1
  res = list()
  
  while(G <= Gmax){
    print(paste("Components = ", G, sep=""))
    
    # Generate initial membership matrix with hierarchical agglomeration
    cl = mclust::hclass(hcTree, G)
    z.ini = matrix(0,N,G)
    for(i in 1:N){
      for(j in 1:G){
        z.ini[i,j] = ifelse(cl[i]==j,1,0)
      }
    }

    # Run MCLUST-ME with chosen initial membership
    out = mcmeVVV(data, z.ini, err, d, itmax, lb)
    res = c(res,out)
    bic.vec[k] = out$BIC
    
    G = G + 1
    k = k + 1
  }
  
  ## Allow user to view results of the optimal model chosen
  n.opt = which(bic.vec==max(bic.vec))
  l = 14*n.opt - 13
  u = 14*n.opt
  res.opt = res[l:u]
  
  ## Output BIC values for each choice of number of clusters
  tbl = cbind(nc,bic.vec)

  out = list(G=nc[n.opt], bic=max(bic.vec), res.opt=res.opt, BIC=tbl,
    nc=nc, bicvec=bic.vec, res=res)
  return(out)
}  
 

## Plot BIC for each choice of number of clusters
plotbic = function(nc.result,title){
  nc = nc.result$nc
  bic.vec = nc.result$bicvec
  
  # get the range for the x and y axis
  xrange = range(nc)
  yrange = range(bic.vec)
  # set up the plot
  plot(xrange, yrange, type="n", xlab="number of components",
     ylab="BIC", xaxt="n", main=title)
  # add lines
  lines(nc, bic.vec, type="b", lwd=1.5, col="red")
  axis(1, at=nc)
}  
  


## Examine difference in model selection results by BIC, with 3 components
test.3group = function(seed,N,tau1,tau2,mu1,mu2,mu3,sig1,sig2,sig3,k,p){
  set.seed(seed)
  data.par = sim.par3(N,tau1,tau2,mu1,mu2,mu3,sig1,sig2,sig3,k,p)
  simdata = sim.data3(data.par)

  dat = simdata$data
  err = simdata$err

  out.mcme = nc.select(dat, err, 1:7)
  out.mclust = Mclust(dat, G=1:7, "VVV")
  
  par(mfrow=c(1,2))
  plotbic(out.mcme,"MCLUST-ME")
  plot(out.mclust,"BIC",legendArgs=list(plot=FALSE),main="MCLUST")
  par(mfrow=c(1,1))
  
  out = list(out.mcme = out.mcme$BIC, out.mclust = out.mclust$BIC)
}


  
# Three well-separated clusters:
bic_wellsep = function(seed){
  N = 300
  tau1 = 0.3
  tau2 = 0.4
  mu1 = c(-16,0)
  mu2 = c(16,0)
  mu3 = c(0,24)
  sig1 = matrix(c(64,0,0,64),nrow=2)
  sig2 = matrix(c(16,0,0,16),nrow=2)
  sig3 = matrix(c(36,0,0,36),nrow=2)
  k = 6
  p = 0.5

  set.seed(seed)
  data.par = sim.par3(N,tau1,tau2,mu1,mu2,mu3,sig1,sig2,sig3,k,p)
  simdata = sim.data3(data.par)

  dat = simdata$data
  err = simdata$err
  z = simdata$z.ini
  group = numeric(N)
  for(i in 1:N){
    if(z[i,1]==1 & z[i,2]==0){
      group[i] = 1
    }else if(z[i,1]==0 & z[i,2]==1){
      group[i] = 2
    }else{
      group[i] = 3
    }
  }
  
  gp1 = dat[group==1,]
  gp2 = dat[group==2,]
  gp3 = dat[group==3,]
  
  xmin = min(dat[,1])
  xmax = max(dat[,1])
  ymin = min(dat[,2])
  ymax = max(dat[,2])
  
  plot(gp1,xlim=c(xmin,xmax),ylim=c(ymin,ymax),pch=0,
    main="Three well-separated clusters",xlab="",ylab="")
  points(gp2,pch=1)
  points(gp3,pch=2)
  
  out = test.3group(seed,N,tau1,tau2,mu1,mu2,mu3,sig1,sig2,sig3,k,p)
  return(out)
}
  
  
# Two clusters close, far away from third one:
bic_2close = function(seed){  
  N = 300
  tau1 = 0.3
  tau2 = 0.4
  mu1 = c(-6,0)
  mu2 = c(6,0)
  mu3 = c(0,36)
  sig1 = matrix(c(64,0,0,64),nrow=2)
  sig2 = matrix(c(16,0,0,16),nrow=2)
  sig3 = matrix(c(36,0,0,36),nrow=2)
  k = 25
  p = 0.5

  set.seed(seed)
  data.par = sim.par3(N,tau1,tau2,mu1,mu2,mu3,sig1,sig2,sig3,k,p)
  simdata = sim.data3(data.par)

  dat = simdata$data
  err = simdata$err
  z = simdata$z.ini
  group = numeric(N)
  for(i in 1:N){
    if(z[i,1]==1 & z[i,2]==0){
      group[i] = 1
    }else if(z[i,1]==0 & z[i,2]==1){
      group[i] = 2
    }else{
      group[i] = 3
    }
  }
  
  gp1 = dat[group==1,]
  gp2 = dat[group==2,]
  gp3 = dat[group==3,]
  
  xmin = min(dat[,1])
  xmax = max(dat[,1])
  ymin = min(dat[,2])
  ymax = max(dat[,2])
  
  plot(gp1,xlim=c(xmin,xmax),ylim=c(ymin,ymax),pch=0,
    main="Two clusters close, one far away",xlab="",ylab="")
  points(gp2,pch=1)
  points(gp3,pch=2)
  
  out = test.3group(seed,N,tau1,tau2,mu1,mu2,mu3,sig1,sig2,sig3,k,p)
  return(out)
}  
  
  
# All three clusters close to each other:  
bic_3close = function(seed){  
  N = 300
  tau1 = 0.3
  tau2 = 0.4
  mu1 = c(-8,0)
  mu2 = c(9,0)
  mu3 = c(0,20)
  sig1 = matrix(c(64,0,0,64),nrow=2)
  sig2 = matrix(c(16,0,0,16),nrow=2)
  sig3 = matrix(c(36,0,0,36),nrow=2)
  k = 36
  p = 0.5

  set.seed(seed)
  data.par = sim.par3(N,tau1,tau2,mu1,mu2,mu3,sig1,sig2,sig3,k,p)
  simdata = sim.data3(data.par)

  dat = simdata$data
  err = simdata$err
  z = simdata$z.ini
  group = numeric(N)
  for(i in 1:N){
    if(z[i,1]==1 & z[i,2]==0){
      group[i] = 1
    }else if(z[i,1]==0 & z[i,2]==1){
      group[i] = 2
    }else{
      group[i] = 3
    }
  }
  
  gp1 = dat[group==1,]
  gp2 = dat[group==2,]
  gp3 = dat[group==3,]
  
  xmin = min(dat[,1])
  xmax = max(dat[,1])
  ymin = min(dat[,2])
  ymax = max(dat[,2])
  
  plot(gp1,xlim=c(xmin,xmax),ylim=c(ymin,ymax),pch=0,
    main="All clusters close",xlab="",ylab="")
  points(gp2,pch=1)
  points(gp3,pch=2)
  
  out = test.3group(seed,N,tau1,tau2,mu1,mu2,mu3,sig1,sig2,sig3,k,p)
  return(out)
}  
  


plot.wellsep = function(seed,bicws){
  N = 300
  tau1 = 0.3
  tau2 = 0.4
  mu1 = c(-16,0)
  mu2 = c(16,0)
  mu3 = c(0,24)
  sig1 = matrix(c(64,0,0,64),nrow=2)
  sig2 = matrix(c(16,0,0,16),nrow=2)
  sig3 = matrix(c(36,0,0,36),nrow=2)
  k = 6
  p = 0.5

  set.seed(seed)
  data.par = sim.par3(N,tau1,tau2,mu1,mu2,mu3,sig1,sig2,sig3,k,p)
  simdata = sim.data3(data.par)

  dat = simdata$data
  err = simdata$err
  z = simdata$z.ini
  group = numeric(N)
  for(i in 1:N){
    if(z[i,1]==1 & z[i,2]==0){
      group[i] = 1
    }else if(z[i,1]==0 & z[i,2]==1){
      group[i] = 2
    }else{
      group[i] = 3
    }
  }
  
  gp1 = dat[group==1,]
  gp2 = dat[group==2,]
  gp3 = dat[group==3,]
  
  xmin = min(dat[,1])
  xmax = max(dat[,1])
  ymin = min(dat[,2])
  ymax = max(dat[,2])
  
  bic.mcme = bicws$out.mcme
  bic.mclust = bicws$out.mclust
  allbic = c(bic.mcme[,2],bic.mclust[1:7])
  ymin.bic = min(allbic)
  ymax.bic = max(allbic)

  plot(gp1,xlim=c(xmin,xmax),ylim=c(ymin,ymax),pch=0,
      main="Case 1: Three well-separated clusters",xlab="",ylab="")
  points(gp2,pch=1)
  points(gp3,pch=2)

  plot(bic.mcme,type="b",lwd=2,ylim=c(ymin.bic,ymax.bic),
    xlab="number of components",ylab="",pch=0,main="BIC values")
  lines(1:7,bic.mclust[1:7],type="b",lty="dashed",pch=2,lwd=2)
  legend("bottomright",legend=c("MCLUST-ME","MCLUST"),
    pch=c(0,2),lty=c("solid","dashed"),cex=0.8)
}



plot.2close = function(seed,bic2){
  N = 300
  tau1 = 0.3
  tau2 = 0.4
  mu1 = c(-6,0)
  mu2 = c(6,0)
  mu3 = c(0,36)
  sig1 = matrix(c(64,0,0,64),nrow=2)
  sig2 = matrix(c(16,0,0,16),nrow=2)
  sig3 = matrix(c(36,0,0,36),nrow=2)
  k = 25
  p = 0.5

  set.seed(seed)
  data.par = sim.par3(N,tau1,tau2,mu1,mu2,mu3,sig1,sig2,sig3,k,p)
  simdata = sim.data3(data.par)

  dat = simdata$data
  err = simdata$err
  z = simdata$z.ini
  group = numeric(N)
  for(i in 1:N){
    if(z[i,1]==1 & z[i,2]==0){
      group[i] = 1
    }else if(z[i,1]==0 & z[i,2]==1){
      group[i] = 2
    }else{
      group[i] = 3
    }
  }
  
  gp1 = dat[group==1,]
  gp2 = dat[group==2,]
  gp3 = dat[group==3,]
  
  xmin = min(dat[,1])
  xmax = max(dat[,1])
  ymin = min(dat[,2])
  ymax = max(dat[,2])
  
  bic.mcme = bic2$out.mcme
  bic.mclust = bic2$out.mclust
  allbic = c(bic.mcme[,2],bic.mclust[1:7])
  ymin.bic = min(allbic)
  ymax.bic = max(allbic)

  plot(gp1,xlim=c(xmin,xmax),ylim=c(ymin,ymax),pch=0,
      main="Case 2: Two clusters close, one far away",xlab="",ylab="")
  points(gp2,pch=1)
  points(gp3,pch=2)

  plot(bic.mcme,type="b",lwd=2,ylim=c(ymin.bic,ymax.bic),
    xlab="number of components",ylab="",pch=0,main="BIC values")
  lines(1:7,bic.mclust[1:7],type="b",lty="dashed",pch=2,lwd=2)
  legend("bottomright",legend=c("MCLUST-ME","MCLUST"),
    pch=c(0,2),lty=c("solid","dashed"),cex=0.8)
}


plot.3close = function(seed,bic3){
  N = 300
  tau1 = 0.3
  tau2 = 0.4
  mu1 = c(-8,0)
  mu2 = c(9,0)
  mu3 = c(0,20)
  sig1 = matrix(c(64,0,0,64),nrow=2)
  sig2 = matrix(c(16,0,0,16),nrow=2)
  sig3 = matrix(c(36,0,0,36),nrow=2)
  k = 25
  p = 0.5

  set.seed(seed)
  data.par = sim.par3(N,tau1,tau2,mu1,mu2,mu3,sig1,sig2,sig3,k,p)
  simdata = sim.data3(data.par)

  dat = simdata$data
  err = simdata$err
  z = simdata$z.ini
  group = numeric(N)
  for(i in 1:N){
    if(z[i,1]==1 & z[i,2]==0){
      group[i] = 1
    }else if(z[i,1]==0 & z[i,2]==1){
      group[i] = 2
    }else{
      group[i] = 3
    }
  }
  
  gp1 = dat[group==1,]
  gp2 = dat[group==2,]
  gp3 = dat[group==3,]
  
  xmin = min(dat[,1])
  xmax = max(dat[,1])
  ymin = min(dat[,2])
  ymax = max(dat[,2])
  
  bic.mcme = bic3$out.mcme
  bic.mclust = bic3$out.mclust
  allbic = c(bic.mcme[,2],bic.mclust[1:7])
  ymin.bic = min(allbic)
  ymax.bic = max(allbic)

  plot(gp1,xlim=c(xmin,xmax),ylim=c(ymin,ymax),pch=0,
      main="Case 3: Three clusters close",xlab="",ylab="")
  points(gp2,pch=1)
  points(gp3,pch=2)

  plot(bic.mcme,type="b",lwd=2,ylim=c(ymin.bic,ymax.bic),
    xlab="number of components",ylab="",pch=0,main="BIC values")
  lines(1:7,bic.mclust[1:7],type="b",lty="dashed",pch=2,lwd=2)
  legend("bottomright",legend=c("MCLUST-ME","MCLUST"),
    pch=c(0,2),lty=c("solid","dashed"),cex=0.8)
}



  