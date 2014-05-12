
### This is an attempt to defining a function for the function to maximize in the M-step.
### Here, we assume the data are 2-dimensional, there are 4 components(clusters),
### and we're given the data as well as the initial classification matrix c.

### We also assume that the covariance matrices of each component have the same shape,
### but potentially different volume and orientation. So the parametrization of Sigma_k is:
###                Sigma_k = lambda_k * D_k %*% A %*% t(D_k)
### according to Celeux and Govaert (1995).

 mix = function(mat,c,data){
               ### Here "mat" is an array of dimension(3,3,4).
               ### Matrices 1-4 contain D_1-D_4 as submatrices, 
   
                  #  D_k |    a     #
                  # _____|_________ #  
                  #      |          # = mat[,,k]
                  # a^T  | lambda_k #
              
               ### "a" is column vector of diagonal elements in A,
               ### "lambda_k" is the k-th sclaing factor.
   
  n=nrow(data)
  d=ncol(data)
  K=ncol(c)
  
  
  ### Define the orthogonal matrices of eigenvectors for each component:
  D = array(dim=c(d,d,K))
  for(k in 1:K){
    D[,,k] = mat[,,k][1:d,1:d]
  }
  
  
  ### Define the diagonal matrix of scaled eigenvalues:
  avec = mat[,,1][1:d,(d+1)]
  A = diag(avec)
  
  
  ### Obtain the scaling factor lambda_k's:
  lam = rep(0,K)
  for (k in 1:K){
    lam[k] = diag(mat[,,k])[d+1]
  }
  
  
  ### Define the covariance matrices for each component:
  Sigma = array(dim=c(d,d,K))
  for (m in 1:K){
    Sigma[,,m] = lam[m] * D[,,m] %*% A %*% t(D[,,m])
  }
  
  
  
  ### Estimate parameters for mean(muhat) and mixing proportion(phat), using formulae
  ### (13)-(15) in Celeux and Govaert (1995):
  clustcount = colSums(c)
  phat = clustcount/n
  muhat = matrix(rep(0,2*K),nrow=2)
  for(j in 1:K){
    muhat[,j] = colSums(data*c[,j])/clustcount[j]
  }
  
  
  ### After rewritting the log-likelihood for complete data, the objective function
  ### simplifies to two components, call them Vsum and Wsum.
  ### Here we parametrize both terms as sum of all elements of two matrices, V and W,
  ### defined as below:
  w = V = matrix(rep(0,n*K),nrow=n)
  onevec = rep(1,K)  
  
  for(i in 1:n){
    for(k in 1:K){
      w[i,k] = c[i,k]*t(data[i,]-muhat[,k])%*%solve(Sigma[,,k])%*%(data[i,]-muhat[,k])
    }
  }
  Wsum = sum(W %*% onevec)
    
  for(i in 1:n){
    for(k in 1:K){
      V[i,k] = c[i,k]*log(det(Sigma[,,k]))
    }
  }
  Vsum = sum(V %*% onevec)
  
  
  ### From the last step, the objective function is thus:
  maxfun = Wsum - Vsum
  
  
  ### Here we specify the constraints on the decomposed matrices D and A.
  ### Since D is supposed to be orthogonal, we require it has orthonormal set
  ### of columns, that is, columns that have unit length and are orthogonal
  ### to one another.
  ### Also, since sigma_k is a cov matrix, we require it to be positive definite,
  ### and therefore all its eigenvalues be positive.
  ### If the conditions are not met, the function is then stopped from running.
   
  if(solve(Sigma[,,1])==t(Sigma[,,1]) & solve(Sigma[,,2])==t(Sigma[,,2]) &
     solve(Sigma[,,3])==t(Sigma[,,3]) & solve(Sigma[,,4])==t(Sigma[,,4]) &
     lam[1]>0 & lam[2]>0 & lam[3]>0 & lam[4]>0 & avec[1]>0 & avec[2]>0) 
    return(maxfun)
  else
    stop()
}

 
 
 
 
 

