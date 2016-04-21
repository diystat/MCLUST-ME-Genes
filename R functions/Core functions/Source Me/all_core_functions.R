

#######################################################################
###------------------------ Core Functions -------------------------###
#######################################################################


### File description:
### This .R file contains all functions responsible for the EM algorithm
### of MCLUST-ME. The highest level wrapper is "mcmeVVV", which requires
### an initial membership matrix and an array of measurement covariances
### as input, and outputs membership probabilities and MLEs upon
### convergence.



#-----------------------Wrapper function-------------------------#

### manual EM algorithm for VVV with measurement error
mcmeVVV = function(data, z, err, d=1, itmax=Inf, lb=1e-3){
  # "d" denotes the diagonal element for initial value of input
  # decomposed matrix, default = 1
  
  # Argument "lb" sets lower bound for diagonal elements of decomposed cov
  # matrices. Default set at 1e-3=0.001. Can be set to a larger number to
  # avoid singularity issues.
  
  # Argument "itmax" sets an upper bound of the number of iterations
  # of the EM algorithm.
  
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
  
  # while loop for iteration:
  repeat{

    ini.par = ini.par.no(data, z, d) # initial values for M-step
    
    thetahat = MstepVVV.err(z, data, err, ini.par, lb) # M-step
    
    temp = EstepVVV.err(thetahat, data, err) # E-step
  
    zhat = temp[[1]] # membership estimates
    
    parameters = temp[[2]] # parameter estimates
    
    loglikelihood = temp[[3]] # records log likelihood
          
    z = zhat # update membership matrix
    
    member = rbind(member,zhat) # store membership estimates
    
    center = rbind(center,parameters$muhat) # store mean estimates
    
    llike[k+1] = loglikelihood # update log likelihood of observed data
              
    k = k+1 # increment k
    
    delta = abs(llike[k-1]-llike[k])/(1+abs(llike[k]))
    
    it = k-2
    print(paste("iteration=",it,sep=""))
    
    if(delta<tol || it>=itmax) break;
  }
  
  error = llike[k]-llike[k-1] # log-lik difference of last two iterations
  
  uncertainty = numeric() # records classification uncertainty of each obs.
  for(i in 1:n){
    rowmax = max(zhat[i,])
    uncertainty[i] = 1-rowmax
  }
  
  nu = G*p*(p+1)/2 + G*p + (G-1)
  bic.me = 2*loglikelihood - nu*log(n) # Calculates BIC for the model
  
  # edit output so it's basically consistent with meVVV() from MCLUST:
  out = list(modelname="VVV with est error", BIC=bic.me, n=n, d=p, G=G, z=zhat, parameters=parameters, uncertainty=uncertainty,
    loglik=loglikelihood, iteration=k-2, error=error, member=member, center=center, likvec=llike[1:k][-(1:2)])
  
  return(out)
}







#-----------------------E-step function--------------------------#

EstepVVV.err = function(param, data, err){
  n = nrow(data)
  p = ncol(data)
  G = ncol(param[[1]])
  
  muhat = param[[1]]
  tauhat = param[[2]]
  sigmahat = param[[3]] # here sigmahat is an array consisting of G p*p cov matrices
  inv.sum = param[[4]]
  
  # set initial z matrix:
  zhat = matrix(rep(0,n*G),nrow=n)
    
  denom = rep(0,n)
  for (i in 1:n){
    for (k in 1:G){
      L = chol(inv.sum[,,k,i])
      temp = as.numeric(data[i,]-muhat[,k])
      denom[i] = denom[i] + tauhat[k]*(det(sigmahat[,,k]+err[,,i]))^(-1/2)*exp(-1/2*crossprod(L%*%temp))
    }
  }
  
  for(i in 1:n){
    for(k in 1:G){
      L = chol(inv.sum[,,k,i])
      temp = as.numeric(data[i,]-muhat[,k])
      num = tauhat[k]*(det(sigmahat[,,k]+err[,,i]))^(-1/2)*exp(-1/2*crossprod(L%*%temp))
      zhat[i,k] = num/denom[i]
    }
  }
  
  parameters = list(muhat=muhat, tauhat=tauhat, sigmahat=sigmahat, inv.sum=inv.sum)
  loglikelihood = log.lik.vvv.err(parameters, data, err) # observed log likelihood
  out = list(zhat, parameters, loglikelihood)
  
  return(out) # output estimated z's, parameter estimates and estimated log-likelihood
}






#-------------------------M-step function------------------------------#

### M-step function for VVV with estimation error
MstepVVV.err = function(z, data, err, ini.par, lb){
  # argument "ini.par" defines a user-specified initial value for L-BFGS-B algorithm
  n = nrow(data)
  p = ncol(data)
  G = ncol(z)
  
  ## obtain mixing proportion estimate:
  clustcount = colSums(z) + .Machine$double.xmin # n_k
  phat = clustcount/n # n_k/n
  
  # set lower bounds for parameters:
  #library(gdata)
  lower.bound = array(-Inf, dim=c(p, p, G))
  lower = numeric()
    for(k in 1:G){
      diag(lower.bound[,,k]) = lb
      lower = c(lower, gdata::lowerTriangle(lower.bound[,,k], diag=TRUE)) # set lower bound for lower triangle elements
    }
  
  ## find arg max for objective function w.r.t. the covariance matrix
  est.res = optim(par=ini.par, obj.fun.VVV.err, gr=NULL, z=z, data=data, err=err,
     lower=lower, method="L-BFGS-B")
   param.est = est.res$par
  
  ## transform back to covariance matrices:
  m = p*(p+1)/2
  sigmahat = array(0, dim=c(p, p, G))
    for(k in 1:G){     
      gdata::lowerTriangle(sigmahat[,,k], diag=TRUE) = param.est[(m*k-m+1):(m*k)]
      sigmahat[,,k] = tcrossprod(sigmahat[,,k])
    }
  
  ## obtain mean estimate:
  inv.sum = inv_sum(z,data,err,sigmahat)
  #print(inv.sum)
  muhat = matrix(0, p, G)
  temp1 = temp2 = 0
      for(k in 1:G){
        for(i in 1:n){         
          tem = inv.sum[,,k,i]
          temp1 = temp1 + z[i,k]*tem
          temp2 = temp2 + z[i,k]*tem%*%data[i,]
        }
        muhat[,k] = chol2inv(chol(temp1))%*%temp2
        temp1 = temp2 = 0
      }
  
  ## put estimates together into a list:
  parameters = list(muhat=muhat, phat=phat, sigma=sigmahat, inv.sum=inv.sum)
  
  return(parameters)
}






#----------------Objective function for M-step-------------------#

### objective function for VVV case, with estimation error
obj.fun.VVV.err =  function(param,z,data,err){
  ### Input:
  ### param --- vector of elements in lower triangular matrix of length Gp(p+1)/2 
  ### z --- matrix of membership probabilities. dimension = n*G. z[i,k]=P(obs. i is in cluster k)
  ### data --- matrix of data. dimension = n*p
  ### err --- array of estimation error matrices, of dimension (p, p, n), one for each observation
  
    n=nrow(data)
    p=ncol(data)
    G=ncol(z)
    length(param) = G*p*(p+1)/2
    piconst = n*p*log(2*pi)/2
  
  ### Construct covariance matrix:
  ### First convert param into lower triangular matrices, using lowerTriangle() from package GDATA
    # library(gdata)
  
    m = p*(p+1)/2
  
    L = array(0, dim=c(p, p, G))
      for(k in 1:G){
        gdata::lowerTriangle(L[,,k],diag=TRUE) = param[(m*k-m+1):(m*k)]
      }
    
  
  ### Then obtain the cov matrices:
    cov.mat = array(0, dim=c(p, p, G))
      for(k in 1:G){
        cov.mat[,,k] = tcrossprod(L[,,k])
      }
  
  ### After rewritting the log-likelihood for complete data, the objective function
  ### becomes:
  ###    F(Sigma) = sum_i sum_k z_ik*t(x_i-muhat_k)*(Sigma_k+Sigma_i)^(-1)*(x_i-muhat_k) + sum_i sum_k z_ik*log(det(Sigma_k+Sigma_i))
  ###
  ### Notice that in this case, MLE for mean is no longer trivial. Instead, the expression for muhat_k is:
  ###    muhat_k = (sum_i z_ik(Sigma_k+Sigma_i)^(-1))^(-1) * sum_i z_ik(Sigma_k+Sigma_i)^(-1)x_i
    
    muhat = matrix(0, p, G)
    inv.sum = array(0, dim=c(p,p,G,n))
    temp1 = matrix(0, p, p)
    temp2 = rep(0, p)
      for(k in 1:G){
        for(i in 1:n){
          #inv.sum[,,k,i] = chol2inv(chol(cov.mat[,,k]+err[,,i]))
          inv.sum[,,k,i] = solve(cov.mat[,,k]+err[,,i])
          temp1 = temp1 + z[i,k]*inv.sum[,,k,i]
          temp2 = temp2 + z[i,k]*inv.sum[,,k,i]%*%data[i,]
        }
        muhat[,k] = chol2inv(chol(temp1)) %*% temp2
        temp1 = matrix(0, p, p)
        temp2 = rep(0, p)
      }
  
  
  ### From the last step, the objective function is thus:
    
    maxfun = 0
      for(i in 1:n){
        for(k in 1:G){
          temp3 = data[i,]-muhat[,k]
          L = chol(inv.sum[,,k,i])
          maxfun = maxfun + z[i,k] * log(det(cov.mat[,,k]+err[,,i])) + z[i,k] * crossprod(L%*%temp3)
        }
      }
  
  ## obtain mixing proportion estimate:
    clustcount = colSums(z) # n_k
    phat = clustcount/n # n_k/n
  
  ## maximizing the likelihood is minimizing its opposite:
    out = (1/2) * maxfun + piconst - sum(colSums(z)*log(phat))
    return(out)
}






#------------------Observed likelihood function---------------------#

## log likelihood of mvn mixture
log.lik.vvv.err = function(param, data, err){
    
  muhat = param[[1]]
  tauhat = param[[2]]
  sigmahat = param[[3]]
  inv.sum = param[[4]]
  
  n = nrow(data)
  p = ncol(data)
  G = ncol(param[[1]])
  
  piconst = -n*p*log(2*pi)/2
  
  temp = rep(0,n)
  lik = 0
  for(i in 1:n){
    for(k in 1:G){
      tem = data[i,]-muhat[,k]
      L = chol(inv.sum[,,k,i])
      temp[i] = temp[i] + tauhat[k] * (det(sigmahat[,,k]+err[,,i]))^(-1/2) * exp((-1/2) * 
          crossprod(L%*%tem))      
    }
    lik = lik + log(temp[i])
  }
  
  out = lik + piconst
  return(out)
}






#----------------------------Other functions----------------------------#

## within cluster scatter matrix:
wkmat = function(z, data){
 
  n = nrow(data)
  p = ncol(data)
  G = ncol(z)
  
  clustcount = colSums(z) + .Machine$double.xmin # prevent a/0 situation
  
  m = p*G
  muhat = matrix(rep(0,m),nrow=p)
    for(k in 1:G){
      muhat[,k] = colSums(data*z[,k])/clustcount[k]
    }
  
  
  dim = c(p,p,G)
  W = array(0, dim=dim)
    
  for(k in 1:G){
    for(i in 1:n){
      temp = data[i,]-muhat[,k]
      W[,,k] = W[,,k] + z[i,k]*tcrossprod(temp)
    }
  }
  
  return(W)
}

#-------------------------------------------------------------------------#

### Function that outputs inverses of sum of covariances for each observation
inv_sum = function(z, data, err, sigmahat){
  
  p = ncol(data)
  n = nrow(data)
  G = ncol(z)
  
  ue = unique(err,MARGIN=3)
  n_unq = dim(ue)[3]
  unique_inv = array(dim=c(p,p,G,n_unq))
  for(k in 1:G){
    for(i in 1:n_unq){
      unique_inv[,,k,i] = chol2inv(chol(sigmahat[,,k]+ue[,,i]))
    }
  }
  
  dim(ue) = c(p*p,n_unq)
  dim(err) = c(p*p,n)
  vec = numeric(n)
  for(j in 1:n_unq){
    vec[colSums(err==ue[,j])==p*p]=j
  }
  
  out = array(dim=c(p,p,G,n))
  for(k in 1:G){
    for(i in 1:n){
      out[,,k,i] = unique_inv[,,k,(vec[i])]
    }
  }
  
  return(out)  
}

#-------------------------------------------------------------------------#

## Generate initial values for M-step:
# Identical errors:
ini.par.iderr = function(data,z,err){
  p = ncol(data)
  G = ncol(z)
  
  clustcount = colSums(z) + .Machine$double.xmin
  wk = wkmat(z, data)

  ini.mat = array(0,dim=c(p,p,G))
    for(k in 1:G){
      ini.mat[,,k] = wk[,,k]/clustcount[k] - err[,,1]
    }
  
  ini.par = numeric()
    for(k in 1:G){
      ini.par = c(ini.par, gdata::lowerTriangle(t(chol(ini.mat[,,k])), diag=TRUE))
    }  
  
  return(ini.par)
}


# Same errors within clusters:
# Caution: Here we assume the unique errors are in order 1,...,G
ini.par.clust = function(data,z,err){
  p = ncol(data)
  G = ncol(z)
  
  clustcount = colSums(z) + .Machine$double.xmin
  wk = wkmat(z, data)
  
  err.unique = unique(err, MARGIN=3)
  ini.mat = array(0,dim=c(p,p,G))
    for(k in 1:G){
      ini.mat[,,k] = wk[,,k]/clustcount[k] - err.unique[,,k]
    }
  
  ini.par = numeric()
    for(k in 1:G){
      ini.par = c(ini.par, gdata::lowerTriangle(t(chol(ini.mat[,,k])), diag=TRUE))
    }  
  
  return(ini.par)
}


# No constraints on error structure. Simply use identity matrix:
ini.par.no = function(data, z, d){
  p = ncol(data)
  G = ncol(z)
  ini.mat = matrix(0, p, p)
  diag(ini.mat) = d
  ini.par = rep(gdata::lowerTriangle(ini.mat, diag=TRUE), G)
  
  return(ini.par)
}