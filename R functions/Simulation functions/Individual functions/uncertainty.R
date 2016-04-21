

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