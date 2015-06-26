
  

    #---- Compare likelihood of each point wrt two dist ----#
    #----           Using MCME or meVVV model           ----#



library(mvtnorm)
library(mclust)
library(gdata)
library(MASS)
library(phyclust)

col.blue = rgb(30,144,255,max=255)



mvndata = function(n,seed,lambda,k1,k2){
  
  mu1 = c(-3,0)
  mu2 = c(3,0)
  I = matrix(c(1,0,0,1),nrow=2)
  E = (lambda^2)*I
  sig1 = matrix(c(4,0,0,4),nrow=2)
  sig2 = sig1 + E  

  # Sample size:
  n1 = k1*n
  n2 = k2*n  

  # Generate samples:
  set.seed(seed); dat1 = mvrnorm(n1, mu1, sig1)
  set.seed(seed); dat2 = mvrnorm(n2, mu2, sig2)
  dat = rbind(dat1,dat2)
  
  out = list(dat=dat,E=E,n=n,k1=k1)
  return(out)
}




lik.compare = function(mvd,filename,method="MCME",it){

  dat = mvd$dat
  k1 = mvd$k1
  n = mvd$n
  
  # Sample size:
  s = nrow(dat)
  n1 = k1*n
  n2 = s - n1

  dat1 = dat[1:n1,]
  dat2 = dat[(n1+1):s,]
   
  # Construct initial classification. Use true membership:
  z.ini = matrix(c(rep(c(1,0),n1),rep(c(0,1),n2)), byrow=T, nrow=s)

  # Construct the true error matrices:
  errmat = array(0,dim=c(2,2,s))
  E = mvd$E
  for(i in (n1+1):s){
    errmat[,,i] = E
  }


  ## One-step iteration of M-step MCME:
  ini.par = ini.par.no(dat, z.ini, 1)
  z = z.ini
  theta = vector("list",it)
  for(i in 1:it){
    param = MstepVVV.err(z, dat, errmat, ini.par, lb=1e-3)
    z = EstepVVV.err(param, dat, errmat)[[1]]
    theta[[i]] = param
  }

  
  ## Component likelihood calculation:
  lik = function(x,mu,S,E){
    out = (2*pi)^(-1)*(det(S+E))^(-1/2)*exp((-1/2)*t(x-mu)%*%solve(S+E)%*%(x-mu))
    return(out)
  }
  
  
  ## Plotting function for scaled likelihoods:
  plot.lkhd = function(theta,it){
    
    xlow = min(dat[,1])-0.5
    xup = max(dat[,1])+0.5
    ylow = min(dat[,2])-0.2
    yup = max(dat[,2])+0.2
    
    ## Mean and covaraince estimates:
    muhat1 = theta$muhat[,1]
    muhat2 = theta$muhat[,2]
  
    S1 = theta$sigma[,,1]
    S2 = theta$sigma[,,2]
    
    p1 = theta$phat[1]
    p2 = theta$phat[2]
    #p1 = p2 = 0.5
  
    lik1 = lik2 = numeric(s)
    if(method=="MCME"){
      for(i in 1:s){
        lik1[i] = p1*lik(dat[i,],muhat1,S1,errmat[,,i])
        lik2[i] = p2*lik(dat[i,],muhat2,S2,errmat[,,i])
      }} else{
      for(i in 1:s){
        lik1[i] = p1*dmvnorm(dat[i,],mean=muhat1,sigma=S1)
        lik2[i] = p2*dmvnorm(dat[i,],mean=muhat2,sigma=S2)
      }}
  

    likratio = -log(lik1/lik2)  

    bdind = which(likratio>-0.2 & likratio<0.2)
    bdind1 = bdind[which(bdind<=n1)]
    bdind2 = bdind[which(bdind>n1)]-n1

    bdind.lt = which(likratio>=0.2)
    bdind.gt = which(likratio<=(-0.2))

    bdind.lt1 = bdind.lt[which(bdind.lt<=n1)]
    bdind.lt2 = bdind.lt[which(bdind.lt>n1)]-n1

    bdind.gt1 = bdind.gt[which(bdind.gt<=n1)]
    bdind.gt2 = bdind.gt[which(bdind.gt>n1)]-n1
  
      
    ## Plot boundary of scaled likelihood equality:
    plot(dat1,xlim = c(xlow,xup), ylim = c(ylow,yup),col="gray",pch=19
      ,xlab="",ylab="",main=paste("Same scaled likelihood (it=",it,")",sep=""))
    points(dat1[bdind1,],col="red",pch=19)
    points(dat2,col="gray",pch=0)
    points(dat2[bdind2,],col=col.blue,pch=0)

    ## Points whose scaled lkhd are larger in the right group:
    plot(dat1,xlim = c(xlow,xup), ylim = c(ylow,yup),col="gray",pch=19
      ,xlab="",ylab="",main=paste("Larger scld lkhd in right group (p1=",round(p1,3),")",sep=""))
    points(dat1[bdind.lt1,],col="red",pch=19)
    points(dat2,col="gray",pch=0)
    points(dat2[bdind.lt2,],col=col.blue,pch=0)

    ## Points whose scaled lkhd are larger in the left group:
    plot(dat1,xlim = c(xlow,xup), ylim = c(ylow,yup),col="gray",pch=19
      ,xlab="",ylab="",main="Larger scld lkhd in left group")
    points(dat1[bdind.gt1,],col="red",pch=19)
    points(dat2,col="gray",pch=0)
    points(dat2[bdind.gt2,],col=col.blue,pch=0)
  }
  
  
  h = 300*it
  png(paste(filename,".png",sep=""),900,height=h)
  par(mfrow=c(it,3))
  for(i in 1:it){
    plot.lkhd(theta[[i]],i)
  }
  par(mfrow=c(1,1))
  dev.off()
  
}


data = mvndata(300,91,3,1,1)
lik.compare(data,"new2",it=5)






