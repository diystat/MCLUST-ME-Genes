

## convert zhat matrix into classification matrix, and compute misclassification rate

MCR = function(z,z.true){
  n = nrow(z)
  G = ncol(z)
  
  for(i in 1:n){
    rowmax = max(z[i,])
    for(k in 1:G){
      z[i,k] = ifelse(z[i,k]==rowmax,1,0)
    }  
  }
  
  diff = z - z.true
  
  out = (sum(rowSums(diff^2))/2)/n
  
  return(out)  
}
