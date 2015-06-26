
library(combinat)

ff = function(n){
  
  rand = function(n,x){
    temp1 = dim(combn(x,2))[2] + dim(combn((n-x),2))[2]
    temp2 = dim(combn(n,2))[2]
  
    out = temp1/temp2
    return(out)
  }

  xx = seq(3,n-3,1)
  yy = sapply(xx,rand,n=n)
  
  plot(xx,yy,type="l",xlab="n1",ylab="Rand")
}

ff(100)







