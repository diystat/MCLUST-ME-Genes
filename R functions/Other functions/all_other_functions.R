


### implementation of fuzzy rand index described in Campello(2006)
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







uncstep = function(data,my.result,filename){
  n = nrow(data)
  palette = c("magenta",heat.colors(10)[2:8],heat.colors(10)[10],"floralwhite")
    
  xlow = min(data[,1])
  xup = max(data[,1])
  ylow = min(data[,2])
  yup = max(data[,2])
  
  xl = seq(0,1,0.1)[1:10]
  yb = rep(0,10)
  xr = seq(0,1,0.1)[2:11]
  yt = rep(0.005,10)
  
  it = my.result$iteration
  
  library("animation")
  oopt <- ani.options(interval = 0.2) # set time between frames
  saveGIF({
  for (i in 1:it) {
    
    z = my.result$member[(n*i+1):(n*i+n),]
    unc = unc(z)
    obs.col = character()
    for(i in 1:n){
      for(k in 2:11){
        if(unc[i]>=(1.1-0.1*k) & unc[i]<(1.2-0.1*k)) obs.col[i] = palette[(k-1)]
      }
    }
    
    par(fig=c(0,1,0,0.8))
    plot(c(xlow,xup),c(ylow,yup),type="n",xlab="",ylab="")
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black")
    points(data,col=obs.col,pch=16,cex=1.5)
  
    par(fig=c(0,1,0.55,1),new=T)
    plot(c(0,1),c(0,0.005),type="n",xlab="",ylab="",
      bty="n",axes=F)
    rect(xl,yb,xr,yt,col=rev(palette),border=NA)
    axis(3,at=seq(1,0,-0.1))

    
    ani.pause()
  }}, movie.name=paste(filename,".gif"), ani.width=1200, ani.height=600)
}  

  
  
  
  

## Plot function for classification uncertainty
uncmap = function(data,unc,filename){
  require(MASS)
  source("uncertainty.R")
  # data---data matrix of dim nxp
  # unc---vector of classification uncertainties
  #       can be obtained using self-defined function "uncertainty()"
  # filename---character string of filename for the plot
  
  n = nrow(data)
  palette = c("magenta",heat.colors(10)[2:8],heat.colors(10)[10],"floralwhite")
  obs.col = character()
  for(i in 1:n){
    for(k in 2:11){
      if(unc[i]>=(1.1-0.1*k) & unc[i]<(1.2-0.1*k)) obs.col[i] = palette[(k-1)]
    }
  }
  
  xlow = min(data[,1])
  xup = max(data[,1])
  ylow = min(data[,2])
  yup = max(data[,2])
  
 
  xl = seq(0,1,0.1)[1:10]
  yb = rep(0,10)
  xr = seq(0,1,0.1)[2:11]
  yt = rep(0.005,10)
  
  
  png(file=paste(filename,".png"),width=1000,height=500)
  par(fig=c(0,1,0,0.8))
  plot(c(xlow,xup),c(ylow,yup),type="n",xlab="",ylab="")
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black")
  points(data,col=obs.col,pch=16,cex=1.5)
  
  par(fig=c(0,1,0.55,1),new=T)
  plot(c(0,1),c(0,0.005),type="n",xlab="",ylab="",
    bty="n",axes=F)
  rect(xl,yb,xr,yt,col=rev(palette),border=NA)
  axis(3,at=seq(1,0,-0.1))
  dev.off()
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






  ## Plot stepwise clustering
  ## Current function only works on two target clusters.

plot.step = function(dat1, dat2, mu1, mu2, z.ini, errmat, itmax, filename){
  # Arguments: mevvv---clustering result object from stepwiseME()
  #            mcme---clustering result object from mcmeVVV()
  #            dat1,dat2---two groups of observations
  #            mu1,mu2---true means of each group
  #            z.ini---initial membership for all observations
  #            itmax---number of iterations
  #            filename---name of gif file produced
  
  col.blue = rgb(30,144,255,max=255)
  
  data = rbind(dat1,dat2)
  
  xlow = min(data[,1])-0.5
  xup = max(data[,1])+0.5
  ylow = min(data[,2])-0.2
  yup = max(data[,2])+0.2
  
  n1 = nrow(dat1)
  n2 = nrow(dat2)
  
  mevvv = stepwiseME(data,z.ini,itmax=itmax)
  mcme = mcmeVVV(data,z.ini,errmat,itmax=itmax)
  
  s = nrow(z.ini)
  G = ncol(z.ini)
  
  mb.mevvv = mevvv$member
  mb.mcme = mcme$member
  
  for(i in 1:nrow(mb.mevvv)){
    rowmax = max(mb.mevvv[i,])
    for(k in 1:G){
      mb.mevvv[i,k] = ifelse(mb.mevvv[i,k]==rowmax,1,0)
    }  
  }
  
  for(i in 1:nrow(mb.mcme)){
    rowmax = max(mb.mcme[i,])
    for(k in 1:G){
      mb.mcme[i,k] = ifelse(mb.mcme[i,k]==rowmax,1,0)
    }  
  }

  it = itmax
  mbary.mevvv = mbary.mcme = array(0,dim=c(s,G,(it+1)))
  for(i in 1:(it+1)){
    mbary.mevvv[,,i] = mb.mevvv[(i*s-s+1):(i*s),]
    mbary.mcme[,,i] = mb.mcme[(i*s-s+1):(i*s),]
  }


  # Estimated centers at each iteration
  cen.mevvv = mevvv$center
  cen.mcme = mcme$center
  p = ncol(data)
  cenary.mevvv = cenary.mcme = array(0,dim=c(p,G,(it+1)))
  for(i in 2:(it+1)){
    cenary.mevvv[,,i] = cen.mevvv[(i*p-p+1):(i*p),]
    cenary.mcme[,,i] = cen.mcme[(i*p-p+1):(i*p),]
  }
  
  cenary.mevvv[,,1] = cenary.mcme[,,1] = cbind(mu1,mu2)


  
  # Make an animation:
  library("animation")
  oopt <- ani.options(interval = 3) # set time between frames
  saveGIF({
  for (i in 1:(it+1)) {
        
    group1 = cbind(dat1,z.ini[1:n1,])
    group2 = cbind(dat2,z.ini[(n1+1):s,])
    
    group1.new = cbind(dat1,mbary.mevvv[1:n1,1,i],mbary.mcme[1:n1,1,i])
    group2.new = cbind(dat2,mbary.mevvv[(n1+1):s,1,i],mbary.mcme[(n1+1):s,1,i])
    
    #group1.mcme = cbind(dat1,mbary.mcme[1:n1,,i])
    #group2.mcme = cbind(dat2,mbary.mcme[(n1+1):s,,i])
    
    par(mfrow=c(1,3),cex=1.5)
    
    plot(dat1, xlim = c(xlow,xup), ylim = c(ylow,yup), pch=ifelse(group1[,3]==1,19,0),
      col=ifelse(group1[,3]==1,col.blue,"red") ,xlab="", ylab="",
      main="True Clusters")
    points(dat2, pch=ifelse(group2[,3]==1,19,0),
      col=ifelse(group2[,3]==1,col.blue,"red"))
    
    plot(dat1, xlim = c(xlow,xup), ylim = c(ylow,yup), pch=ifelse(group1.new[,3]==1,19,0),
      col=ifelse(group1.new[,3]==group1.new[,4],"gray",ifelse(group1.new[,3]==1,col.blue,"red")),
      xlab="", ylab="", main=paste("meVVV, iteration =",i-1))
    points(dat2, pch=ifelse(group2.new[,3]==1,19,0),
      col=ifelse(group2.new[,3]==group2.new[,4],"gray",ifelse(group2.new[,3]==1,col.blue,"red")))
    points(cenary.mevvv[1,,i],cenary.mevvv[2,,i],col=c("lawngreen","slateblue"),pch=4,cex=3,lwd=5)
    
    plot(dat1, xlim = c(xlow,xup), ylim = c(ylow,yup), pch=ifelse(group1.new[,4]==1,19,0),
      col=ifelse(group1.new[,3]==group1.new[,4],"gray",ifelse(group1.new[,4]==1,col.blue,"red")),
      xlab="", ylab="", main=paste("MCME, iteration =",i-1))
    points(dat2, pch=ifelse(group2.new[,4]==1,19,0),
      col=ifelse(group2.new[,3]==group2.new[,4],"gray",ifelse(group2.new[,4]==1,col.blue,"red")))
    points(cenary.mcme[1,,i],cenary.mcme[2,,i],col=c("lawngreen","slateblue"),pch=4,cex=3,lwd=5)
       
    par(mfrow=c(1,1),cex=1)
    
    ani.pause()
  }}, movie.name=paste(filename,".gif"), ani.width=1800, ani.height=600)
}
  



## Step-by-step iteration for MCLUST VVV algorithm:
stepwiseME = function(data, z, itmax=Inf){
  ## data---data matrix
  ## z---initial membership matrix
  
  n = nrow(data)
  p = ncol(data)
  G = ncol(z)
  
  member = z # matrix for storing membership estimates
  center = matrix(0,p,G) # matrix for storing mean estimates

  k = 2 # keeps track of number of iteration
  it = 1

  loglikelihood = NA
  parameters = list()
  zhat = matrix(0,n,G)
  
  errvec = numeric()

  llike = rep(0, 1000) # set convergence criterion
  #llike[1] = -10000
  #llike[2] = llike[1] + 1e-4
  FLMAX = 1.7976931348623157e308
  llike[1] = FLMAX/2
  llike[2] = FLMAX
 
  tol = 1e-5
  # while loop for iteration:
  repeat {
    
    #print(paste("iteration =",k-1)) # prints number of iterations
    
    thetahat = mstepVVV(data, z)$parameters # M-step
    
    temp = estepVVV(data, thetahat) # E-step
  
    zhat = temp$z # membership estimates
    
    muhat = temp$parameters$mean # mean estimates
    
    parameters = temp$parameters
    
    loglikelihood = temp$loglik # records log likelihood
          
    z = zhat # update membership matrix
    
    member = rbind(member,zhat) # store membership estimates
    
    center = rbind(center,muhat) # store mean estimates
    
    llike[k+1] = loglikelihood # update log likelihood of observed data
    
    errvec = c(errvec,llike[k]-llike[k-1])
    
    #print(loglikelihood) # prints current evaulated observed data log-likelihood
      
    k = k+1 # increment k
    
    it = it+1
    
    err = abs(llike[k-1]-llike[k])/(1+abs(llike[k]))
    
    if(err<tol) break;
  }
  
  error = llike[k]-llike[k-1]
  
  # edit output so it's basically consistent with meVVV() from MCLUST:
  out = list(modelname="VVV with est error", n=n, d=p, G=G, z=zhat, member=member, center=center,
    parameters=parameters, loglik=loglikelihood, iteration=k-2, error=error, likvec=llike[1:k][-(1:2)],
    errvec=errvec)
  
  return(out)
  
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




prop.plot = function(tau,seed=4){
  N = 200
  mu1 = c(0,0)
  mu2 = c(8,0)
  sig1 = matrix(c(64,0,0,64),nrow=2)
  sig2 = matrix(c(16,0,0,16),nrow=2)
  k = 36
  p = 0.3

  data.par = sim.par(N,tau,mu1,mu2,sig1,sig2,k,p)
  set.seed(seed)
  simdata = sim.data(data.par)
  rand.samples = simdata$data
  z.ini = simdata$z.ini
  index = simdata$index

  col.blue = rgb(30,144,255,max=255)  
  xlow = min(rand.samples[,1])-0.5
  xup = max(rand.samples[,1])+0.5
  ylow = min(rand.samples[,2])-0.2
  yup = max(rand.samples[,2])+0.2
  
  rs.new = cbind(rand.samples,index,z.ini[,1])  
  plot(rand.samples,pch=ifelse(rs.new[,3]==1,2,16),main=paste("tau=",tau,sep=""),
    col=ifelse(rs.new[,4]==1,col.blue,"red"),xlab="",ylab="",xlim=c(xlow,xup),ylim=c(ylow,yup))
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






meVVV.fc = function(data, z,itmax=Inf, mu1, mu2){
  # Argument "errstr" is user-specified error structure.
  # If errstr="identical", all errors are the same.
  # If errstr="cluster", errors are the same within each cluster.
  # If errstr="none", no constraints on error structure.
  # When errstr="none", "d" denotes the diagonal element, default = 1
  
  # Argument "lb" sets lower bound for diagonal elements of decomposed cov
  # matrices. Default set at 1e-3=0.001. Can be set to a larger number to
  # avoid singularity issues.
  
  source("E-step_VVV_error.R")
  source("M-step_VVV_error.R")
  source("log likelihood_VVV_error.R")
  source("objective function_VVV_error.R")
  source("W_k matrix.R")
  source("inipar.R")
  
  n = nrow(data)
  p = ncol(data)
  G = ncol(z)
  
  member = z # matrix for storing membership estimates
  center = matrix(0,p,G) # matrix for storing mean estimates

  k = 2 # keeps track of number of iteration

  piconst = n*p*log(2*pi)/2
  
  loglikelihood = NA
  parameters = list()
  zhat = matrix(0,n,G)

  llike = rep(0, 1000) # set convergence criterion
  FLMAX = 1.7976931348623157e308
  llike[1] = FLMAX/2
  llike[2] = FLMAX
 
  tol = 1e-5
  
  mu.fix = cbind(mu1,mu2)
  
  
  # while loop for iteration:
  repeat{
    
    print(paste("iteration =",k-1)) # prints number of iterations
      
    thetahat = mstepVVV(data, z) # M-step
    thetahat$parameters$mean = mu.fix
    
    temp = estepVVV(data, thetahat$parameters) # E-step
  
    zhat = temp$z # membership estimates
    
    parameters = temp$parameters # parameter estimates
    
    loglikelihood = temp$loglik # records log likelihood
          
    z = zhat # update membership matrix
          
    k = k+1 # increment k
    
    delta = abs(llike[k-1]-llike[k])/(1+abs(llike[k]))
    
    it = k-2
    
    if(delta<tol || it>=itmax) break;
  }
  
  error = llike[k]-llike[k-1]
  
  # edit output so it's basically consistent with meVVV() from MCLUST:
  out = list(modelname="VVV", n=n, d=p, G=G, z=zhat, parameters=parameters,
    loglik=loglikelihood, iteration=k-2, error=error)
  
  return(out)
}






mcmeVVV.fc = function(data, z, err, errstr="none", d=1, itmax=Inf, lb=1e-3, mu1, mu2){
  # Argument "errstr" is user-specified error structure.
  # If errstr="identical", all errors are the same.
  # If errstr="cluster", errors are the same within each cluster.
  # If errstr="none", no constraints on error structure.
  # When errstr="none", "d" denotes the diagonal element, default = 1
  
  # Argument "lb" sets lower bound for diagonal elements of decomposed cov
  # matrices. Default set at 1e-3=0.001. Can be set to a larger number to
  # avoid singularity issues.
  
  source("E-step_VVV_error.R")
  source("M-step_VVV_error.R")
  source("log likelihood_VVV_error.R")
  source("objective function_VVV_error.R")
  source("W_k matrix.R")
  source("inipar.R")
  
  n = nrow(data)
  p = ncol(data)
  G = ncol(z)
  
  member = z # matrix for storing membership estimates
  center = matrix(0,p,G) # matrix for storing mean estimates

  k = 2 # keeps track of number of iteration

  piconst = n*p*log(2*pi)/2
  
  loglikelihood = NA
  parameters = list()
  zhat = matrix(0,n,G)

  llike = rep(0, 1000) # set convergence criterion
  FLMAX = 1.7976931348623157e308
  llike[1] = FLMAX/2
  llike[2] = FLMAX
 
  tol = 1e-5
  
  mu.fix = cbind(mu1,mu2)
  
  
  # while loop for iteration:
  repeat{
    
    print(paste("iteration =",k-1)) # prints number of iterations
    
     # initial values for M-step
    if(errstr=="identical"){
      ini.par = ini.par.iderr(data, z, err)
    } else if(errstr=="cluster"){
      ini.par = ini.par.clust(data, z, err)
    } else if(errstr=="none"){
      ini.par = ini.par.no(data, z, d)
    }
      
    thetahat = MstepVVV.err(z, data, err, ini.par, lb) # M-step
    thetahat$muhat = mu.fix
    
    temp = EstepVVV.err(thetahat, data, err) # E-step
  
    zhat = temp[[1]] # membership estimates
    
    parameters = temp[[2]] # parameter estimates
    
    loglikelihood = temp[[3]] # records log likelihood
          
    z = zhat # update membership matrix
    
    member = rbind(member,zhat) # store membership estimates
    
    center = rbind(center,parameters$muhat) # store mean estimates
    
    llike[k+1] = loglikelihood # update log likelihood of observed data
    
    #print(paste("loglik =",loglikelihood))
          
    k = k+1 # increment k
    
    delta = abs(llike[k-1]-llike[k])/(1+abs(llike[k]))
    
    it = k-2
    
    if(delta<tol || it>=itmax) break;
  }
  
  error = llike[k]-llike[k-1]
  
  uncertainty = numeric() # records classification uncertainty of each obs.
  for(i in 1:n){
    rowmax = max(zhat[i,])
    uncertainty[i] = 1-rowmax
  }
  
  # edit output so it's basically consistent with meVVV() from MCLUST:
  out = list(modelname="VVV with est error", n=n, d=p, G=G, z=zhat, parameters=parameters, uncertainty=uncertainty,
    loglik=loglikelihood, iteration=k-2, error=error, member=member, center=center, likvec=llike[1:k][-(1:2)])
  
  return(out)
}




## Function that compares outputs of two functions, meVVV and stepwiseME
mecompare = function(data,z){
  temp1 = meVVV(data,z)
  temp2 = stepwiseME(data,z)
  itmax = min(as.numeric(attributes(temp1)$info[1]),as.numeric(temp2$it)) # number of iterations meVVV used
  
  pro = mean = sigma = result = list()
  
  for(i in 1:itmax){ 
    res.step = stepwiseME(data,z,itmax=i)
    pro.step = res.step$par$pro
    mean.step = res.step$par$mean
    sigma.step = res.step$par$var$sigma
    
    res.mevvv = meVVV(data,z,control=emControl(itmax=i))
    pro.mevvv = res.mevvv$par$pro
    mean.mevvv = res.mevvv$par$mean
    sigma.mevvv = res.mevvv$par$var$sigma
    
    pro = list(pro.step=pro.step,pro.mevvv=pro.mevvv)
    mean = list(mean.step=mean.step,mean.mevvv=mean.mevvv)
    sigma = list(sigma.step=sigma.step,sigma.mevvv=sigma.mevvv)
    
    combined = list(iteration=i,pro=pro,mean=mean,sigma=sigma)
    
    result = append(result,combined)
  }
    return(result)
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






sim.driver = function(N,tau,mu1,mu2,sig1,sig2,k,p,nseed){
  
  source('~/Research Project/Simulations/Sim/Functions/par_data_run.R')
  source('~/Research Project/R functions for MCME(VVV)/Other functions/fuzzyrand.R')
  source('~/Research Project/Simulations/Sim/Functions/rand.R')
  source('~/Research Project/Simulations/Sim/Functions/boundary.R')
  
  simres = list()
  randres = matrix(,nrow=nseed,ncol=8)
  
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

  #par(mfrow=c(1,2))
  #plot(rand.samples,pch=ifelse(rs.new[,4]==1,2,16),main=paste("True Grouping,seed=",seed,sep=""),
  #  col=ifelse(rs.new[,5]==1,col.blue,"red"),xlab="",ylab="",xlim=c(xlow,xup),ylim=c(ylow,yup),cex=1.6)
  
  #plot(rand.samples,pch=ifelse(rs.new[,4]==1,2,16),main=paste("MCLUST-ME,seed=",seed,sep=""),
  #  col=ifelse(rs.new[,3]==1,col.blue,"red"),xlab="",ylab="",xlim=c(xlow,xup),ylim=c(ylow,yup),cex=1.6)
  #plot(rand.samples,pch=ifelse(rs.new[,4]==1,2,16),main="MCLUST-ME",cex.main=2,
  #  col=ifelse(rs.new[,3]==1,col.blue,"red"),xlab="",ylab="",xlim=c(xlow,xup),ylim=c(ylow,yup),cex=1.7)
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
  #plot(rand.samples,pch=ifelse(rs.new1[,4]==1,2,16),main=paste("MCLUST,seed=",seed,sep=""),
  #  col=ifelse(rs.new1[,3]==1,col.blue,"red"),xlab="",ylab="",xlim=c(xlow,xup),ylim=c(ylow,yup),cex=1.6)
  plot(rand.samples,pch=ifelse(rs.new1[,3]==1,22,24),main="MCLUST",cex.main=1,
    bg=ifelse(rs.new1[,4]==1,"black","white"),xlab="",ylab="",xlim=c(xlow,xup),ylim=c(ylow,yup),cex=0.7)
  contour(x,y,ze,level=0,add=TRUE,drawlabels=FALSE,lwd=2)
  #par(mfrow=c(1,1))
  
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
 
  #par(mfrow=c(1,2))
  #plot(rand.samples,pch=ifelse(rs.new[,4]==1,2,16),main=paste("True Grouping,seed=",seed,sep=""),
  #  col=ifelse(rs.new[,5]==1,col.blue,"red"),xlab="",ylab="",xlim=c(xlow,xup),ylim=c(ylow,yup),cex=1.6)
  
  #plot(rand.samples,pch=ifelse(rs.new[,4]==1,2,16),main=paste("MCLUST-ME,seed=",seed,sep=""),
  #  col=ifelse(rs.new[,3]==1,col.blue,"red"),xlab="",ylab="",xlim=c(xlow,xup),ylim=c(ylow,yup),cex=1.6)
  #plot(rand.samples,pch=ifelse(rs.new[,4]==1,2,16),main="MCLUST-ME",cex.main=2,
  #  col=ifelse(rs.new[,3]==1,col.blue,"red"),xlab="",ylab="",xlim=c(xlow,xup),ylim=c(ylow,yup),cex=1.7)
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
  #plot(rand.samples,pch=ifelse(rs.new1[,4]==1,2,16),main=paste("MCLUST,seed=",seed,sep=""),
  #  col=ifelse(rs.new1[,3]==1,col.blue,"red"),xlab="",ylab="",xlim=c(xlow,xup),ylim=c(ylow,yup),cex=1.6)
  plot(rand.samples,pch=ifelse(rs.new1[,3]==1,22,24),main="MCLUST",cex.main=1,
    bg=ifelse(rs.new1[,4]==1,"black","white"),xlab="",ylab="",xlim=c(xlow,xup),ylim=c(ylow,yup),cex=pointsize2)
  contour(x,y,ze,level=0,add=TRUE,drawlabels=FALSE,lwd=2)
  #par(mfrow=c(1,1))
  
} 





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





## Examine difference in model selection results by BIC, with 3 components
test.3group = function(seed){

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

  out.mcme = nc.select(dat, err, 1:7)
  out.mclust = Mclust(dat, G=1:7, "VVV")
  
  par(mfrow=c(1,2))
  plotbic(out.mcme,"MCLUST-ME")
  plot(out.mclust,"BIC",legendArgs=list(plot=FALSE),main="MCLUST")
  par(mfrow=c(1,1))
  
  out = list(out.mcme = out.mcme$BIC, out.mclust = out.mclust$BIC)
  
}



  
  