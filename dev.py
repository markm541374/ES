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
ff = fgen1d(lower, upper, 800, kf)


# %%

f = ff.genfun()
plt.plot(sp.linspace(-1, 1, 100), map(f, sp.linspace(-1, 1, 100)))
ee = lambda x, y: (f(x), 0)
[xmintrue, miny, ierror] = DIRECT.solve(ee, lower, upper, user_data=[], algmethod=1, maxf=4000)
xmintrue = xmintrue[0]
print miny
plt.plot(xmintrue, miny, 'rx')

# %%

n_init = 5
x = sp.random.uniform(-1, 1, n_init)
y = map(f, x)+sp.random.normal(scale=0.01, size=n_init)
Xo = sp.matrix(x).T
Yo = sp.matrix(y).T
So = sp.matrix([[0.000001]]*n_init)
Do = [[sp.NaN]]*n_init
g1 = GPd.GPcore(Xo, Yo, So, Do, kf)
a = GPd.plot1(g1, [-1], [1])
# %%


def sqexpprior(xx):
    p = sps.norm.pdf(((xx[0])-0.)/2.) * sps.norm.pdf(((xx[1])-0.)/2.)
    return sp.log(p)



s = 0.000000001


# %%

for k in xrange(1):
    e = EntropyPredict.EntPredictor([Xo, Yo, So, Do], [-1], [1], kfGen, sqexpprior)
    e.MLEsearchn = 200
    e.HYPsamplen = 50
    e.ENTsearchn = 160
    
    n = 300
    Xi = sp.linspace(-1, 1, n)
    Si = sp.logspace(-4, 0, 10)
    a = e.showEntGraph(Xi, [s], plot='basic')
    xmin = e.searchAtS([-1], [1], [s])
    a.plot(xmin, [0], 'go')
    a.plot(xmin, f(xmin[0]), 'ro')
    xIn = xmin[0]
    yIn = f(xIn)+sp.random.normal(scale=s)

    Xo = sp.vstack([Xo, xIn])
    Yo = sp.vstack([Yo, yIn])
    So = sp.vstack([So, s])
    Do.append([sp.NaN])

    
    plt.show()
    #hh = EntropyPredict.inferHmulti(e.Pred, sp.matrix(Xi).T, [s], plot='verbose')


# %%
reload(GPep)
reload(EntropyPredict)

# %%
Xi = sp.linspace(0.1,0.5,150)
e = EntropyPredict.EntPredictor([Xo, Yo, So, Do], [-1], [1], kfGen, sqexpprior)
e.drawHypSamples(1,plot='verbose')

e.initPredictor(plot='verbose')
a=e.showEntGraph(Xi,[0.00000000001])
a.axis([0.15,0.5,-1.5,-0.5])

# %%


hh = EntropyPredict.inferHmulti(e.Pred, sp.matrix(Xi).T, [0.000001], plot='verbose')

# %%
kfm = GPep.gen_sqexp_k_d([5.5280215860924278, 2.1461063922860951])
g1 = GPd.GPcore(Xo, Yo, So, Do, kfm)
a = GPd.plot1(g1, [-1], [1])
