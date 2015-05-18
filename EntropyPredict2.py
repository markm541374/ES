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
         
        self.nHYPsamples = para['nHYPsamples']
        self.HYPsearchLow = para['HYPsearchLow']
        self.HYPsearchHigh = para['HYPsearchHigh']
        self.HYPMLEsearchn = para['HYPMLEsearchn']
        self.HYPsamSigma = para['HYPsamSigma']
        self.HYPsamBurn = para['HYPsamBurn']
        self.ENTnsam = para['ENTnsam']
        self.ENTzeroprecision = para['ENTzeroprecision']
        self.ENTsearchn = para['ENTsearchn']
        return
    
    def __del__(self):
        try:
            print 'closing FBinfer'
            self.FBInfer.close()
        except:
            print traceback.format_exc()
        try:
            print 'closing EPinfer'
            self.EPInfer.close()
        except:
            print traceback.format_exc()
        return
        
    def setupEP(self):
        self.searchMLEHYP()
        self.drawHYPsamples()
        self.initFBInfer()
        self.drawmins()
        self.initEPInfer()
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
    
    def initEPInfer(self):
        Xo=self.D[0]
        Yo=self.D[1]
        So=self.D[2]
        Do=self.D[3]
        self.EPInfer = GPset.multiGP()
        self.mask = sp.zeros(self.nHYPsamples)
        for i in xrange(self.nHYPsamples):
            md = self.ENTmindraws[i]
            if not md[0]==0:
                self.mask[i]=-1
                #create a bad GP in the correct position so that others are alligned
                self.EPInfer.addGPd(-1,-1,-1,-1,-1)
                continue
            xs = md[1][0,0]
            Xg = sp.matrix([[xs]])
            Yg = sp.matrix([[0.]])
            Sg = sp.matrix([[self.ENTzeroprecision]])
            Dg = [[0]]
            
            [Xc, Yc, Sc, Dc] = GPep.catObs([[Xo, Yo, So, Do], [Xg, Yg, Sg, Dg]])
            
            Xz = sp.matrix([[xs], [xs]])
            Dz = [[sp.NaN], [0, 0]]
            # the inequality
            Iz = sp.matrix([[Yo[Yo.argmin(), :][0, 0]], [0.]])
            # sign of the inequality
            Gz = [0, 0]
            Nz = [So[Yo.argmin(), :], 0.] #!!!!!this value is important, should it e the sigma for hte min obs or the posterior at that piont??
        
            self.EPInfer.addGPep(Xc, Yc, Sc, Dc, Xz, Dz, Iz, Gz, Nz, self.HYPsampleFns[i])
            
        
        return
        
        
    def inferMLEpost(self,X_s,D_s):
        m,v = self.MLEInfer.infer_diag(X_s,D_s)
        return [m,v]
        
    def plotMinDraws(self):
        #explicitly 1d
        h = []
        for j in self.ENTmindraws:
            if not j[0]==0:
                continue
            h.append(j[1][0,0])
        f = plt.figure()
        a = f.add_subplot(111)
        a.hist(h,40)
        return [f,a]        
        
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
    
    def infer_both(self, X_s, D_s):
        infr_before = self.FBInfer.infer_diag(sp.matrix(X_s).T,D_s)
        infr_after = self.EPInfer.infer_diag(sp.matrix(X_s).T,D_s)
        return [infr_before, infr_after]
    
    def findMV(self, X_s, D_s):
        i1 = self.FBInfer.infer_diag(X_s, D_s)
        currentmask=sp.zeros(self.nHYPsamples)
        Vydx = sp.zeros(self.nHYPsamples)
        Mydx = sp.zeros(self.nHYPsamples)
        for i in xrange(self.nHYPsamples):
            if not i1[i][0]==0:
                currentmask[i]=-1
                continue
            Vydx[i] = i1[i][1][1][0,0]
            Mydx[i] = i1[i][1][0][0,0]
        
        Vydxxs = sp.zeros(self.nHYPsamples)
        Mydxxs = sp.zeros(self.nHYPsamples)
        Xmcs=[]
        Dmcs=[]
        for i in xrange(self.nHYPsamples):
            X = sp.vstack([X_s,self.ENTmindraws[i][1]])
            Xmcs.append(X)
            Dmcs.append(D_s+[[sp.NaN]])
        i2 = self.EPInfer.infer_full_var(Xmcs,Dmcs)
        for i in xrange(self.nHYPsamples):
            if not i2[i][0]==0:
                currentmask[i]=-2
                continue
           
            V = i2[i][1][1]
            
            m = i2[i][1][0]
            #print V
            #<magic>
            s = V[0, 0]+V[1, 1]-V[0, 1]-V[1, 0]
            if s<=10**-10:
                print 's<10**10'
                print s
            mu = m[1, 0]-m[0, 0] 
            alpha = - mu / sp.sqrt(s) # sign difference compared to the paper because I am minimizing
        
            beta = sp.exp(sps.norm.logpdf(alpha) - sps.norm.logcdf(alpha))
            coef = beta*(beta+alpha)*(1./s)*(V[0, 0]-V[0, 1])**2
            vnxxs = V[0, 0]-coef
            Vydxxs[i]=vnxxs
            Mydxxs[i]=m[0,0]
            #</magic>
        return [Mydx, Vydx, Mydxxs, Vydxxs, currentmask]
        
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
        y2_r=sp.zeros(np)
        mv_r=sp.zeros(np)
        
        n_c = len([i for i in clean if i==0])
        for i,y in enumerate(y_s):
            
            if clean[i]==0:
                for j in xrange(np):
                    y_r[j]+=y[j]
                    y2_r[j]+=y[j]**2
                    mv_r[j]+=v_s[i][j]
        #mean of the mean
        y_r = [y/float(n_c) for y in y_r]
        #mean of the variance
        mv_r = [v/float(n_c) for v in mv_r]
        #mean of the square of hte mean
        y2_r = [y/float(n_c) for y in y2_r]
        #variance of the mean
        vy_r = [ y2-y_r[i]**2 for i,y2 in enumerate(y2_r)]
        #posterior mean is mean of mean, posterior variance os mean of variance plus variance of mean
        v_r = [mv+vy_r[i] for i,mv in enumerate(mv_r)]
        return [y_r,v_r]
        
    def findENT(self, xs, ds, ss):
        n_hyp = self.nHYPsamples
        [m0,v0,m1,v1,mask] = self.findMV(sp.matrix(xs).T,[ds])
        H0 = sp.zeros(n_hyp)
        H1 = sp.zeros(n_hyp)
        for i in xrange(n_hyp):
            if not mask[i]==0:
                continue
            H0[i] = 0.5*sp.log(2*sp.pi*sp.e*(v0[i]+ss))
            H1[i] = 0.5*sp.log(2*sp.pi*sp.e*(v1[i]+ss))

        H=(sum(H0)-sum(H1))/float(sum([1 for i in mask if i==0]))
        return H
        
    def plotEPchanges(self, axis=0,point='None',np=100,obstype=[[sp.NaN]]):
        print 'plotting EPchanges'
        X=[]
        if point=='None':
            point=sp.zeros(self.dim)
        n_hyp = len(self.FBInfer.processes)
        clean = sp.zeros(n_hyp)
        x_r = sp.linspace(self.lb[axis],self.ub[axis],np)
        m_0 = [sp.zeros(np) for i in xrange(n_hyp)]
        v_0 = [sp.zeros(np) for i in xrange(n_hyp)]
        m_1 = [sp.zeros(np) for i in xrange(n_hyp)]
        v_1 = [sp.zeros(np) for i in xrange(n_hyp)]
        for i,x in enumerate(x_r):
            pi = point.copy()
            pi[axis]=x
            res = self.findMV(sp.matrix(pi),[[sp.NaN]])
            
            for j in xrange(n_hyp):
                if res[4][j]==0:
                    m_0[j][i] = res[0][j]
                    v_0[j][i] = res[1][j]
                    m_1[j][i] = res[2][j]
                    v_1[j][i] = res[3][j]
        split = int(sp.ceil(sp.sqrt(n_hyp)))
       
        [f,axs] = plt.subplots(split,split,sharex='col', sharey='row')
        i=0;j=0;k=0
        while i<n_hyp:
            axs[j][k].plot(x_r,m_0[i],'b')
            axs[j][k].plot(x_r,m_1[i],'r')
            axs[j][k].set_ylim([-3,3])
            i+=1
            j=i/split
            k=i%split
            
        
        return [f,axs]
    def searchENTs(self,ss, obstype=[sp.NaN]):
        print 'Searhing for MaxENT'
        def ee(x,y):
            global ENTsearchi
            ENTsearchi+=1
            ENTsearchi
            print '\rIter: %d    ' % ENTsearchi,
            return (-self.findENT(x, obstype,ss), 0)
        global ENTsearchi
        ENTsearchi=0
        [xmin, miny, ierror] = DIRECT.solve(ee, self.lb, self.ub, user_data=[], algmethod=1, maxf=self.ENTsearchn, logfilename='/dev/null')
        del(ENTsearchi)
        print 'maxENT '+str(xmin)
        return [xmin, miny]
            
    def plotENT(self,ss,axis=0,point='None',np=100,obstype=[[sp.NaN]]):
        print 'plotting predictive Entropy'
        [f,a] = self.plotFBpost(axis=axis, point=point,np=np,obstype=obstype)
        X=[]
        if point=='None':
            point=sp.zeros(self.dim)
        n_hyp = len(self.FBInfer.processes)
        
        x_r = sp.linspace(self.lb[axis],self.ub[axis],np)
        
        for i in x_r:
            pi = point.copy()
            pi[axis]=i
            X.append(pi)
        
        H=sp.zeros(np)
        for i in xrange(np):
            x=X[i]
            h=self.findENT(x,obstype[0],ss)
            
            H[i]=h
        a2 = a.twinx()
        a2.plot(x_r,H,'r')
        return [f,[a,a2]]
        
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
        #this bit is 1D
        xs = sp.array(self.D[0]).flatten()
        ys = sp.array(self.D[1]).flatten()
        a.plot(xs, ys, 'rx')
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
        self.ENTmindraws = res
        return res
        
class Optimizer():
    def __init__(self, f, kfGen, kfPrior, lb, ub, para):
        self.f = f
        self.kfGen = kfGen
        self.kfPrior = kfPrior
        self.lb = lb
        self.ub = ub
        self.para=para
        
        return
        
    def initrandobs(self, n_init, s_init):
        x = sp.random.uniform(self.lb, self.ub, n_init)
        y = map(self.f, x)+sp.random.normal(scale=sp.sqrt(s_init), size=n_init)
        self.Xo = sp.matrix(x).T
        self.Yo = sp.matrix(y).T
        self.So = sp.matrix([[s_init]]*n_init)
        self.Do = [[sp.NaN]]*n_init
        return

    def initspecobs(self, Xo, Yo, So, Do):
        self.Xo = Xo
        self.Yo = Yo
        self.So = So
        self.Do = Do
        return
        
    def setupEP(self):
        try:
            del(self.EP)
        except:
            pass
        self.EP = EntPredictor([self.Xo,self.Yo,self.So,self.Do], self.lb, self.ub, self.kfGen, self.kfPrior, self.para )
        self.EP.setupEP()
        return
        
    def plotstate(self):
        [f0,a0] = self.EP.plotHYPsamples(d0=0, d1=1)
        [f1,a1] = self.EP.plotFBpost()
        [f2,a2] = self.EP.plotMLEpost()
        [f3,a3] = self.EP.plotMinDraws()
        [f4,as4] = self.EP.plotENT(0.01,np=100)
        return
    
    def searchnextFixS(self,s,obstype=[sp.NaN]):
        [x_n, H_e] = self.EP.searchENTs(s,obstype=obstype)
        yIn = self.f(x_n)+sp.random.normal(scale=sp.sqrt(s))
        self.Xo = sp.vstack([self.Xo, x_n])
        self.Yo = sp.vstack([self.Yo, yIn])
        self.So = sp.vstack([self.So, s])
        self.Do.append(obstype)
        return