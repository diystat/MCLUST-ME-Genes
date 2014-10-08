
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
      ini.par = c(ini.par, lowerTriangle(t(chol(ini.mat[,,k])), diag=TRUE))
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
      ini.par = c(ini.par, lowerTriangle(t(chol(ini.mat[,,k])), diag=TRUE))
    }  
  
  return(ini.par)
}


# No constraints on error structure. Simply use identity matrix:
ini.par.no = function(data, z){
  p = ncol(data)
  G = ncol(z)
  
  ini.mat = matrix(0, p, p)
  diag(ini.mat) = 0.5
  ini.par = rep(lowerTriangle(ini.mat, diag=TRUE), G)
  
  return(ini.par)
}

