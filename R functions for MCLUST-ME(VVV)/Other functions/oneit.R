
oneit = function(data,z,err){
    ini.par = ini.par.no(data, z, 1)    
    thetahat = MstepVVV.err(z, data, err, ini.par, 1e-3) # M-step    
    temp = EstepVVV.err(thetahat, data, err) # E-step
}





