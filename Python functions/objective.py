# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 00:18:36 2015

@author: wzhang
"""
import numpy as np
import math
import scipy as sp
import numpy.linalg as nl

def objectiveFunction(param,z,data,err):
    """
    Note on input: param is a (G,m) array, with each row 
    representing lower-triangular elements (including diagonal)
    of covariance matrix for each component
    """
    param = np.asarray(param)    
    n = np.shape(data)[0]
    p = np.shape(data)[1]
    G = np.shape(z)[1]
    m = int(p*(p+1)/2)
    par2D = np.reshape(param,(G,m))
    piConstant = n*p*math.log(2*math.pi)/2
    
    """
    Construct covariance matrix
    """    
    L = np.zeros(shape=(G,p,p))
    covMat = np.zeros(shape=(G,p,p))
    for g in range(G):
        for k in range(p):
            lower = int(k*p-k**2/2+k/2)
            upper = int(lower+p-k)
            L[g,range(k,p),range(0,p-k)] = par2D[g,lower:upper]
        covMat[g] = np.dot(L[g],L[g].T)  
    
    """
    Compute mean estimate
    """
    muhat = np.zeros(shape=(p,G))
    invSum = np.zeros(shape=(n,G,p,p))
    temp1 = np.zeros(shape=(p,p))
    temp2 = np.zeros(p)
    for g in range(G):
        for i in range(n):
            invSum[i,g] = nl.inv(covMat[g]+err[:,:,i])
            temp1 = temp1 + z[i,g]*invSum[i,g]
            temp2 = temp2 + z[i,g]*np.dot(invSum[i,g],data[i,:].T)
        muhat[:,g] = np.dot(nl.inv(temp1),temp2)
        temp1 = np.zeros(shape=(p,p))
        temp2 = np.zeros(p)
        
    """
    Objective function
    """
    maxFun = 0.0
    for g in range(G):
        for i in range(n):
            temp3 = data[i,:].T-muhat[:,g]
            L = np.dot(np.dot(temp3.T,invSum[i,g]),temp3)
            maxFun = maxFun + z[i,g]*math.log(nl.det(covMat[g]+err[:,:,i])) + z[i,g]*L
           
    """
    Mixing proportion estimate
    """
    clusterCount = np.sum(z,axis=0)
    phat = clusterCount/n+10**(-6)
    
    out = maxFun/2 + piConstant - np.sum(clusterCount*np.log(phat))
    print(out)
    return out


"""
Returns lower bounds on parameters
"""
def cons(param,z,data,err):
    G = np.shape(z)[1]
    p = np.shape(data)[1]
    m = int(p*(p+1)/2)
    par2D = np.reshape(np.asarray(param),(G,m))[:,:p]
    vec = par2D.ravel()-10**(-5)
    return vec