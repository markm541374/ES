import scipy as sp
from matplotlib import pyplot as plt
import scipy.stats as sps
import scipy.linalg as spl
from scipy.stats import norm as norm
import GPep
import GPdc as GPd
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
    def __init__(self,fname,lastinvalid=False):
        try:
            self.O = EntropyPredict.restartOpt(fname)
        except:
            print 'opt did not complete cleanly, loadding to penultimate step.'
            self.O = EntropyPredict.restartOpt(fname,lastinvalid=lastinvalid)
        return
        
    def xerr(self):
        xerr=[]
        xmintrue = sp.matrix(self.O.para['xmintrue'])
        
        for s in self.O.states[1:]:
            xminest = sp.matrix(s['xminIR'])
            xerr.append(spl.norm(xmintrue-xminest, ord='fro'))
        return xerr
    
    def yIR(self):
        print 'YYYY'
        yIR=[]
        ymintrue = self.O.para['ymintrue']
        for s in self.O.states[1:]:
             
            xminest = sp.matrix(s['xminIR'])
            print self.O.f([0,0])
            yIR.append(abs(self.O.f(sp.array(xminest).flatten())-ymintrue))
        
        return sp.array(yIR)


    
    def region(self):
        reg=[]
        for s in self.O.states[1:]:
            reg.append(s['region_radius'])
        return reg
        
    def d3f(self):
        reg=[]
        for s in self.O.states[1:]:
            reg.append(max(s['err_d3f']))
        return reg
        
    def global_hyp(self):
        reg=[]
        for s in self.O.states[1:]:
            reg.append(s['global_hyp'])
        return reg
    
    def aqu(self):
        reg=[]
        for s in self.O.states[1:]:
            reg.append(s['searchres'][4])
        return reg
        
    def dvdf(self):
        reg=[]
        for s in self.O.states[1:]:
            reg.append(max(s['err_dvdf']))
        return reg
        
    def steptimes(self):
        t=[]
        for s in self.O.states[1:]:
            t.append(s['time'])
        return sp.array(t)
        
    def s_chosen(self):
        reg=[]
        for s in self.O.states[1:]:
            reg.append(s['searchres'][2])
        return reg
    def x_chosen(self):
        reg=[]
        for s in self.O.states[1:]:
            reg.append(s['searchres'][0])
        return reg
    def y_chosen(self):
        reg=[]
        for s in self.O.states[1:]:
            reg.append(s['searchres'][1])
        return reg
    def u_chosen(self,s2u=None):
        reg=[]
        for s in self.O.states[1:]:
            try:
                reg.append(s['u'])
            except KeyError:
                reg.append(s2u(s['searchres'][2]))
        return reg
        
    def MLEhyp(self):
        h=[]
        for s in self.O.states[1:]:
            h.append([10**x for x in s['logHYPMLE']])
        return h
        
    def sprofile(self):
        
        S=[]
        for s in self.O.states[1:]:
            S.append(s['sprofdata'])
        
        return S
        
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
