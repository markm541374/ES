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
import traceback

class EntPredictor():
    def __init__(self, D, lower, upper, kfgen, kfprior, para):
        self.D = D
        self.dim = len(lower)
        self.kfgen = kfgen
        self.kfprior = kfprior
        
        self.lb=lower
        self.ub=upper
         
        self.nHYPsamples = para[0]
        self.HYPsearchLow = para[1]
        self.HYPsearchHigh = para[2]
        self.HYPMLEsearchn = para[3]
        self.HYPsamSigma = para[4]
        self.HYPsamBurn = para[5]
        self.ENTnsam = para[6]
        return
    
    def __del__(self):
        try:
            print 'closing FBinfer'
            self.FBInfer.close()
        except:
            print traceback.format_exc()
        return
        
    def setupEP(self):
        self.searchMLEHYP()
        self.drawHYPsamples()
        self.initFBInfer()
        return
    
    def searchMLEHYP(self):
        print 'seeking MLE phyperparameters'
    
        dist = dist_obj(squaresllk, self.D, self.kfprior, self.kfgen)
        searchn = self.HYPMLEsearchn
        
        ee = lambda hyp, y: (-dist.loglike(hyp), 0)
        [xmintrue, miny, ierror] = DIRECT.solve(ee, self.HYPsearchLow, self.HYPsearchHigh, user_data=[], algmethod=1, maxf=searchn, logfilename='/dev/null')

        print 'MLEhyperparameters: '+str([10**i for i in xmintrue])
        self.logMLEHYPVal = xmintrue
        self.MLEHYPFn = self.kfgen([10**i for i in xmintrue])
        self.MLEInfer = GPd.GPcore(self.D[0], self.D[1], self.D[2], self.D[3], self.MLEHYPFn)
        return
        
    def drawHYPsamples(self):
        n=self.nHYPsamples
        print 'Drawing '+str(n)+' hyperparameter samples'
        w_0 = self.logMLEHYPVal
        sigma = self.HYPsamSigma*sp.ones(self.dim+1)
        dist = dist_obj(squaresllk, self.D, self.kfprior, self.kfgen)
        
        samples = slice_sample(dist, w_0, iters=n, sigma=sigma, burn=self.HYPsamBurn)
        
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
        try:
            a.plot(10**self.logMLEHYPVal[d0], 10**self.logMLEHYPVal[d1],'bo')
        except:
            pass
        return [f,a]
        
    def initFBInfer(self):
        print 'seting up FB inference'
        self.FBInfer = GPset.multiGP()
        for kf in self.HYPsampleFns:
            self.FBInfer.addGPd(self.D[0], self.D[1], self.D[2], self.D[3], kf)
        status = self.FBInfer.status()
        ns = len([i for i in status if i[0]==0])
        nk = self.nHYPsamples
        print str(ns)+' of '+str(nk)+' kernel draws inited sucessfuly'
        return
    
    def inferMLEpost(self,X_s,D_s):
        
        m,v = self.MLEInfer.infer_diag(X_s,D_s)
        return [m,v]
    
    def plotMLEpost(self,axis=0,point='None',np=100,obstype=[[sp.NaN]]):
        print 'plotting MLEpost'
        X=[]
        if point=='None':
            point=sp.zeros(self.dim)
       
        x_r = sp.linspace(self.lb[axis],self.ub[axis],np)
        
        for i in x_r:
            pi = point.copy()
            pi[axis]=i
            X.append(pi)
        [m,v] = self.inferMLEpost(sp.matrix(X),obstype*np)
        
        u_b=sp.zeros(np)
        l_b=sp.zeros(np)
        for j in xrange(np):
            u_b[j]=m[j]+2*sp.sqrt(v[j])
            l_b[j]=m[j]-2*sp.sqrt(v[j])
        f = plt.figure()
        a = f.add_subplot(111)
        a.plot(x_r,m,'b')
        a.fill_between(x_r, l_b, u_b, facecolor='lightskyblue', alpha=0.5)
        a.set_title('MLE')
        return [f,a]
        
    def inferFBpost(self,X_s, D_s):
        #np=len(X_s)
        infr = self.FBInfer.infer_diag(sp.matrix(X_s).T,D_s)
        
        n_hyp = len(self.FBInfer.processes)
        clean = sp.zeros(n_hyp)
        y_s = [[] for i in xrange(n_hyp)]
        v_s = [[] for i in xrange(n_hyp)]
        for j,r in enumerate(infr):
            if r[0]==0:
                y_s[j]=sp.array(r[1][0]).flatten()
                v_s[j]=sp.array(r[1][1]).flatten()
            else:
                
                clean[j]=-1
        
        
        np=len(y_s[0])
        y_r=sp.zeros(np)
        v_r=sp.zeros(np)
        n_c = len([i for i in clean if i==0])
        for i,y in enumerate(y_s):
            
            if clean[i]==0:
                for j in xrange(np):
                    y_r[j]+=y[j]
                    v_r[j]+=v_s[i][j]
        
        y_r = [y/float(n_c) for y in y_r]
        v_r = [v/float(n_c) for v in v_r]
        
        return [y_r,v_r]
        
    def plotFBpost(self,axis=0,point='None',np=100,obstype=[[sp.NaN]]):
        print 'plotting FBpost'
        X=[]
        if point=='None':
            point=sp.zeros(self.dim)
        n_hyp = len(self.FBInfer.processes)
        clean = sp.zeros(n_hyp)
        x_r = sp.linspace(self.lb[axis],self.ub[axis],np)
        y_s = [[] for i in xrange(n_hyp)]
        v_s = [[] for i in xrange(n_hyp)]
        for i in x_r:
            pi = point.copy()
            pi[axis]=i
            X.append(pi)
        [m,v] = self.inferFBpost(sp.matrix(X).T,obstype*np)
        
        
        u_b=sp.zeros(np)
        l_b=sp.zeros(np)
        for j in xrange(np):
            u_b[j]=m[j]+2*sp.sqrt(v[j])
            l_b[j]=m[j]-2*sp.sqrt(v[j])
        f = plt.figure()
        a = f.add_subplot(111)
        a.plot(x_r,m,'b')
        a.fill_between(x_r, l_b, u_b, facecolor='lightskyblue', alpha=0.5)
        a.set_title('FB')
        return [f,a]
    
    def EIMLE(self, X_s):
        X_s=sp.matrix(X_s)
        np = X_s.shape[0]
        D_s = [[sp.NaN]]*np
        m,v = self.inferMLEpost(X_s, D_s)
        E=sp.zeros(np)
        best = self.D[1].min()
        
        for i in xrange(np):
            E[i] = EI(best, m[i],sp.sqrt(v[i]))
        return E

    def EIFB(self, X_s):
        X_s=sp.matrix(X_s)
        np = X_s.shape[0]
        D_s = [[sp.NaN]]*np
        m,v = self.inferFBpost(sp.matrix(X_s).T, D_s)
        E=sp.zeros(np)
        best = self.D[1].min()
        
        for i in xrange(np):
            E[i] = EI(best, m[i],sp.sqrt(v[i]))
        return E

    def drawmins(self):
        res = self.FBInfer.drawmins(self.ENTnsam,[self.lb,self.ub])
        return res