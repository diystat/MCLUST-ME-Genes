
### Objective function for M-step.

### Here, we assume the data are 2-dimensional, there are 4 components(clusters),
### and we're given the data as well as the initial classification matrix z.

### We also assume that the covariance matrices of each component have the same shape,
### same volume and same orientation. So the parametrization of Sigma_k is:
###                Sigma_k = lambda * D %*% A %*% t(D)
### according to Celeux and Govaert (1995).

obj.fun = function(param,z,data){
  ### Input:
  ### param --- vector of elements in lower triangular matrix of length p(p+1)/2 
  ### z --- matrix of membership probabilities. dimension = n*K. z[i,k]=P(obs. i is in cluster k)
  ### data --- matrix of data. dimension = n*d.
  
    n=nrow(data)
    p=ncol(data)
    G=ncol(z)
    length(param) = p*(p+1)/2
  
  ### Construct covariance matrix:
  ### First convert param into a lower triangular matrix, using lowerTriangle() from package GDATA
    library(gdata)
    L = matrix(0, p, p)
    lowerTriangle(L,diag=TRUE) = param
    # print(L)
  
  ### Then obtain the cov matrix:
    cov.mat = L %*% t(L)
    # print(cov.mat)
  
  
  ### After rewritting the log-likelihood for complete data, the objective function
  ### becomes:
  ###    F(Sigma) = trace(W%*%Sigma^-1) + n*log(det(Sigma))
  ### where W = sum_i sum_k z(i,k)(x_i-mu_k)(x_i-mu_k)^T
    
    W = wmat(z, data) # see function wmat in "W matrix.R"
  
  
  ### From the last step, the objective function is thus:
    maxfun = sum(diag(W%*%solve(cov.mat))) + n*log(det(cov.mat))
  
    return(maxfun)
    
}



#-------------------------------------------------------------------------------------------#



## debugging:
obj.fun.test = function(){
  
    m = 100;
    n = 3;
    G = 3;
    
    piconst = n*m*log(2*pi)/2

    set.seed(999);
    x = rnorm(m*n);
    dim(x) = c(m, n);

    z = matrix(0, m, G);
    z[,1] = 1;

    ## Test 1
    ## covariance matrix is diagonal:

    ## Use your code
    mat = diag(1:n);
    param = lowerTriangle(mat, diag=TRUE)
    temp = obj.fun(param, z, x)
    loglik = -temp/2 - piconst
    loglik


    ## Use dmvnorm
    mu =  colMeans(x);
    library(mvtnorm);
    d1 = sum(dmvnorm(x, mu, mat^2, log=TRUE))
  
    d1 == loglik # results match

  
    ## transform x:
    set.seed(999);
    A = matrix(rnorm(n*n), n, n);
    A = A / (abs(det(A))^(1/n));
    print(A);
    x1 = t(A %*% t(x))
    mu1 = t(A %*% mu);
    mat1 = A %*% (mat^2) %*% t(A);

    ## use my own code:   
    param1 = lowerTriangle(t(chol(mat1)), diag=TRUE)
    temp1 = obj.fun(param1, z, x1);
    loglik1 = -temp1/2 - piconst
    loglik1

    loglik1 == loglik
    loglik1 == d1 # results match

 
}


