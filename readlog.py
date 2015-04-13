# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 21:54:09 2015

@author: mark
"""

import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
import scipy.stats as sps
import scipy.linalg as spl
from scipy import interpolate
import time
import GPep
import GPd
from tools import *
import EntropyPredict
import DIRECT
import copy as cp
import pickle
%pylab inline

# %%


    result = pickle.load(open(fname, 'rb'))
    xmintrue = result['xmintrue']
    xminest=[abs(i[0]-xmintrue) for i in result['minest']]
    plt.semilogy(xminest)
    
# %%
def plotxminerr(res):
    nr=len(res)
    
    xminerr=[]
    plt.figure()
    for r in res:
        mt=r['xmintrue']
        
        xm=[abs(i[0]-mt) for i in r['minest']]
        xminerr.append(xm)
        plt.semilogy(xm,'b')
    np=len(xminerr[0])
    
    m = sp.zeros(np)
    m2 = sp.zeros(np)
    v = sp.zeros(np)
    
    for i in xrange(np):
        for j in xrange(nr):
            m[i]+=xminerr[j][i]
            m2[i]+=(xminerr[j][i])**2
        v[i]=sp.sqrt(m[i]**2-m2[i])
    plt.figure()
    plt.semilogy(m,'b')
    plt.semilogy(m+v,'r')
    plt.semilogy(m-v,'r')
    return
    
    
    
def readfiles(fnames):
    res = []
    for f in fnames:
        res.append(pickle.load(open(f, 'rb')))
    return res

# %%


fnames=[]
for i in xrange(12):
    fnames.append('res/test'+str(i)+'.obj')
    
r=readfiles(fnames)
plotxminerr(r)