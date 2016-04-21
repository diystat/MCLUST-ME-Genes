


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

