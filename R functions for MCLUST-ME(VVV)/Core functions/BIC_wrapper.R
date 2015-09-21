

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
  res.opt = res[n.opt]
  
  ## Output BIC values for each choice of number of clusters
  tbl = cbind(nc,bic.vec)

  out = list(G=nc[n.opt], bic=max(bic.vec), result=res.opt, BIC=tbl,
    nc=nc, bicvec=bic.vec)
  return(out)
}  
 

## Plot BIC for each choice of number of clusters
plotbic = function(nc.result){
  nc = nc.result$nc
  bic.vec = nc.result$bicvec
  
  # get the range for the x and y axis
  xrange = range(nc)
  yrange = range(bic.vec)
  # set up the plot
  plot(xrange, yrange, type="n", xlab="number of components",
     ylab="BIC", xaxt="n", main="BIC, larger is better")
  # add lines
  lines(nc, bic.vec, type="b", lwd=1.5)
  axis(1, at=nc)
}  
  