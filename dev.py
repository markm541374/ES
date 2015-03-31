# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 14:15:14 2015

@author: markm
"""
# %%

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
%pylab inline

# %%

hyptrue = [1., 0.25]
kfGen = GPep.gen_sqexp_k_d
kf = GPep.gen_sqexp_k_d(hyptrue)
upper = 1
lower = -1
ff = fgen1d(lower, upper, 500, kf)

def sqexpprior(xx):
    p = sps.norm.pdf(((xx[0])-0.)/2.) * sps.norm.pdf(((xx[1])-0.)/2.)
    return sp.log(p)

# %%

f = ff.genfun()

plt.plot(sp.linspace(-1, 1, 100), map(f, sp.linspace(-1, 1, 100)))
ee = lambda x, y: (f(x), 0)
[xmintrue, miny, ierror] = DIRECT.solve(ee, lower, upper, user_data=[], algmethod=1, maxf=4000)
xmintrue = xmintrue[0]
print miny
plt.plot(xmintrue, miny, 'rx')

# %%

n_init = 12
x = sp.random.uniform(-1, 1, n_init)
y = map(f, x)+sp.random.normal(scale=0.01, size=n_init)
Xo = sp.matrix(x).T
Yo = sp.matrix(y).T
So = sp.matrix([[10**-16]]*n_init)
Do = [[sp.NaN]]*n_init
g1 = GPd.GPcore(Xo, Yo, So, Do, kf)
a = GPd.plot1(g1, [-1], [1])


# %%


# %%
o = EntropyPredict.Optimizer(f,kfGen,sqexpprior,[-1],[1])
o.initspecobs(Xo,Yo,So,Do)
#o.showinov(0.000001, obstype=[0])
o.searchats(10**-16, method='Ent')

# %%

for i in xrange(15):
    o.searchats(10**-16,method='EI')

# %%
#reload
#reload(GPd)
#reload(GPep)
#reload(EntropyPredict)
e = EntropyPredict.EntPredictor([Xo, 100*Yo, 100*So, Do], [-1], [1], kfGen, sqexpprior)
e.HYPsamplen = 25
e.drawHypSamples(25, plot='verbose')


# %%
reload(EntropyPredict)

def sqexpprior2(xx):
    p = sps.norm.pdf(((xx[0])-0.)/2.) * sps.norm.pdf(((xx[1])-0.)/2.)
    return sp.log(p)


o2 = EntropyPredict.Optimizer(f,kfGen,sqexpprior,[-1],[1])
o2.initspecobs(o.Xo,10*o.Yo,100*o.So,o.Do)
#o.showinov(0.000001, obstype=[0])
o2.searchats(10**-14, method='Ent')