# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 14:27:45 2015

@author: markm
"""
import scipy as sp
from scipy import linalg as spl
from scipy import interpolate
import GPd
from scipy.stats import norm as norms
import scipy.stats as sps
import numpy as np
from functools import partial
import sys
from pyDOE import lhs as lhs
import traceback
class MJMError(Exception):
    pass

def EI(ER,mu,sigma):
        alpha=(-ER+mu)/sigma
        Z = norms.cdf(-alpha)
        if Z==0.0:
            return sp.matrix(0.0)
	#print "alpha: "+str(alpha)
	#print "Z: "+str(Z)
        E=-mu+norms.pdf(alpha)*sigma/Z+ER
        ei=Z*E
        if np.isfinite(ei):
            return sp.matrix(ei)
        else:
            return sp.matrix(0.0)

class fgen1d():
    def __init__(self, lower, upper, npoints, kf):
        self.npoints = npoints
        self.Xi = sp.matrix(sp.linspace(lower, upper, npoints)).T
        V_p = GPd.buildKsym_d(kf, self.Xi, [[sp.NaN]]*npoints)
        self.vpCho = spl.cholesky(V_p, lower=True)
        
        return

    def genfun(self):
        fn = self.vpCho*sp.matrix(sp.random.normal(size=self.npoints)).T
        f = interpolate.interp1d(sp.array(self.Xi).flatten(), sp.array(fn).flatten(),kind='cubic')
        return f


class fgennd():
    #assumes the hypercube [-1,1]^D
    def __init__(self, D, npoints, kf):
        self.npoints = npoints

        support = lhs(D,samples=npoints,criterion='centermaximin')
        self.Xi = sp.matrix([[y*2-1 for y in x] for x in support])
        V_p = GPd.buildKsym_d(kf, self.Xi, [[sp.NaN]]*npoints)
        self.vpCho = spl.cholesky(V_p, lower=True)
        self.kf = kf
        return

    def genfun(self):
        fn = self.vpCho*sp.matrix(sp.random.normal(size=self.npoints)).T
        g = GPd.GPcore(self.Xi,fn,sp.matrix([10e-6]*self.npoints).T,[[sp.NaN]]*self.npoints,self.kf)
        f = lambda x: sp.array(g.infer_m(sp.matrix(x),[[sp.NaN]])[0,0])
        return f

def slice_sample(dist, init, iters, sigma, step_out=True,burn=10):
    """
    from http://isaacslavitt.com/2013/12/30/metropolis-hastings-and-
    slice-sampling/
    with some changes
    based on http://homepages.inf.ed.ac.uk/imurray2/teaching/09mlss/
    """

    # dist = joint_dist()

    # set up empty sample holder
    D = len(init)
    samples = sp.zeros((D, iters+burn))
    
    # initialize
    xx = init.copy()

    for i in xrange(iters+burn):
        mn=i-burn+1
        print '\r Drawn %d    ' % mn,
        sys.stdout.flush()

        perm = range(D)
        sp.random.shuffle(perm)
        
        last_llh = dist.loglike(xx)

        for d in perm:
            llh0 = last_llh + sp.log(sp.random.rand())
            rr = sp.random.rand(1)
            x_l = xx.copy()
            x_l[d] = x_l[d] - rr * sigma[d]
            x_r = xx.copy()
            x_r[d] = x_r[d] + (1 - rr) * sigma[d]

            if step_out:
                llh_l = dist.loglike(x_l)
                while llh_l > llh0:
                    # print x_l
                    x_l[d] = x_l[d] - sigma[d]
                    llh_l = dist.loglike(x_l)
                llh_r = dist.loglike(x_r)
                while llh_r > llh0:
                    x_r[d] = x_r[d] + sigma[d]
                    llh_r = dist.loglike(x_r)

            x_cur = xx.copy()
            while True:
                xd = sp.random.rand() * (x_r[d] - x_l[d]) + x_l[d]
                x_cur[d] = xd.copy()
                last_llh = dist.loglike(x_cur)
                if last_llh > llh0:
                    xx[d] = xd.copy()
                    break
                elif xd > xx[d]:
                    x_r[d] = xd
                elif xd < xx[d]:
                    x_l[d] = xd
                else:
                    raise RuntimeError('Slice sampler shrank too far.')

        samples[:, i] = xx.copy().ravel()
    print ''
    return samples[:,burn:]

def squaresllk(xx, D, kfgen):
            hyp = [10**(i) for i in xx]
            kf = kfgen(hyp)
            try:
                g = GPd.GPcore(D[0], D[1], D[2], D[3], kf)
                lk = g.llk()
            except:
                traceback.print_exc()
                print 'warn: can''t get llk with: '+str(hyp)
                lk=-sp.Inf
            return lk


class dist_obj():
    def __init__(self, llk, D, prior, kfgen):
        self.llk = llk
        self.D = D
        self.prior = prior
        self.kfgen = kfgen
        return

    def loglike(self, xx):
        
        P = self.prior(xx)
        L = self.llk(xx, self.D, self.kfgen)
        
        # print str(L)+'  xx  '+str(P)
        # print L+P
        return L+P



def sqExpPrior(paras, xx):
    
    dim=len(paras)
    p = 1
    for i in xrange(dim):
        #do not allow lengthscales that will cause zero correlation numerical error between points separated by 0.1
        #this is 0.004 on dunvegan
        if i>0 and sp.exp(-10**(-2*xx[i])*0.01)==0:
            return -sp.Inf
        mul = sps.norm.pdf(xx[i], loc=paras[i][0], scale = sp.sqrt(paras[i][1]))
        p*=mul
        
    return sp.log(p)
    
def genSqExpPrior(para):
    # print 'genk: '+str(theta)
    k = partial(sqExpPrior,para)
    return k