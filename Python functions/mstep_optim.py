# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 01:12:22 2015

@author: wzhang
"""
import numpy as np
import scipy.optimize as opt
import rpy2.robjects as robjects
import numpy.random as rnd
import timeit

robjects.r['load']("RNASeq.RData")        


obs = np.array(robjects.r['obs'])
err = np.array(robjects.r['errary'])
       
G = 3       
z = np.zeros(shape=(1000,G))
ind = rnd.random_integers(0,G-1,(1000))
for i in range(1000):
    z[i,ind[i]] = 1
    
p = 5
m = int(G*p*(p+1)/2)        
iniPar = np.ones(m)   
np.seterr(all='ignore')

start = timeit.default_timer()     
opt.fmin_slsqp(objectiveFunction,x0=iniPar,f_ieqcons=cons,args=(z,obs,err))
stop = timeit.default_timer()
print(stop-start)



    