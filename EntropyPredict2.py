# -*- coding: utf-8 -*-
"""
Created on Thu May  7 16:30:38 2015
Rework of EntropyPredict to use multiprocessing and hopefully be cleaner
@author: mark
"""
import scipy as sp
from matplotlib import pyplot as plt
import scipy.stats as sps
import scipy.linalg as spl
from scipy.stats import norm as norm
import GPep
import GPd
import GPset
# from scipy.optimize import minimize as mnz
import DIRECT
from functools import partial
# import time
from multiprocessing import pool
from tools import *
import sys

class EntPredictor():
    def __init__(self, D, lower, upper, kfgen, kfprior, para):
        self.D = D
        self.dim = len(lower)
        self.kfgen = kfgen
        self.kfprior = kfprior
        
        self.lb=lower
        self.ub=upper
         
        self.nHYPsamples = para[0]
        return
        
    def setupEP(self):
        self.searchMLEHYP()
        self.drawHYPsamples()
        return
    
    def searchMLEHYP(self):
        print 'seeking MLE phyperparameters'
    
        dist = dist_obj(squaresllk, self.D, self.kfprior, self.kfgen)
        searchn = 500
        
        ee = lambda hyp, y: (-dist.loglike(hyp), 0)
        [xmintrue, miny, ierror] = DIRECT.solve(ee, [-3, -3], [3, 3], user_data=[], algmethod=1, maxf=searchn, logfilename='/dev/null')

        print 'MLEhyperparameters: '+str([10**i for i in xmintrue])
        self.logMLEHYPVal = xmintrue
        self.MLEHYPFn = self.kfgen([10**i for i in xmintrue])
        return
        
    def drawHYPsamples(self):
        n=self.nHYPsamples
        print 'Drawing '+str(n)+' hyperparameter samples'
        w_0 = self.logMLEHYPVal
        sigma = 0.05*sp.ones(self.dim+1)
        dist = dist_obj(squaresllk, self.D, self.kfprior, self.kfgen)
        
        samples = slice_sample(dist, w_0, iters=n, sigma=sigma)
        # print samples
        kfSam = []
        hySam = []
        for i in xrange(n):
            kfSam.append(self.kfgen([10**samples[0][i], 10**samples[1][i]]))
            hySam.append([10**samples[0][i], 10**samples[1][i]])
        
        self.HYPsampleVals=hySam
        self.HYPsampleFns=kfSam
        
        return
    
    def plotHYPsamples(self, d0=0, d1=1):
        #plts hyperparameter samples. defaults to first two values but can specify others
        f = plt.figure()
        a = f.add_subplot(111)
        for hs in self.HYPsampleVals:
            a.plot(hs[d0],hs[d1],'bx')
        a.set_yscale('log')
        a.set_xscale('log')
        a.set_xlabel(str(d0)+'st hyperparameter')
        a.set_ylabel(str(d1)+'st hyperparameter')
        return [f,a]
        
        