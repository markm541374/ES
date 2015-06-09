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
import dill as pickle
import EntropyPredict

class OptEval():
    def __init__(self,fname):
        self.O = EntropyPredict.restartOpt(fname)
        return
        
    def xerr(self):
        xerr=[]
        xmintrue = sp.matrix(self.O.para['xmintrue'])
        for s in self.O.states[1:]:
            xminest = sp.matrix(s['xminIR'])
            xerr.append(spl.norm(xmintrue-xminest, ord='fro'))
        return xerr
    
    def yIR(self):
        yIR=[]
        ymintrue = self.O.para['ymintrue']
        for s in self.O.states[1:]:
            xminest = sp.matrix(s['xminIR'])
            yIR.append(abs(self.O.f(xminest)[0,0]-ymintrue))
        
        return sp.array(yIR)
    
    def steptimes(self):
        t=[]
        for s in self.O.states[1:]:
            t.append(s['time'])
        return sp.array(t)
    
def plotset(Y,x,f,a,c,extras=[]):
    
    if f==None:
        f = plt.figure()
        a = f.add_subplot(111)
    for i in xrange(Y.shape[0]):
        a.plot(x,Y[i,:],c)
    for e in extras:
        if e['name']=='median':
            plt.plot(x,getmedian(Y),e['colorcode'])
        if e['name']=='mean':
            plt.plot(x,getmean(Y),e['colorcode'])
    return [f,a]
    
def getmedian(X):
    #median along each column
    m=[]
    for j in xrange(X.shape[1]):
        m.append(sp.median(X[:,j]))
    return m
    
def getmean(X):
    #median along each column
    m=[]
    for j in xrange(X.shape[1]):
        m.append(sp.mean(X[:,j]))
    return m