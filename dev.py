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

hyptrue = [1., 0.2]
kfGen = GPep.gen_sqexp_k_d
kf = GPep.gen_sqexp_k_d(hyptrue)
upper = 1
lower = -1
ff = fgen1d(lower, upper, 1001, kf)


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
So = sp.matrix([[0.001]]*n_init)
Do = [[sp.NaN]]*n_init
g1 = GPd.GPcore(Xo, Yo, So, Do, kf)
a = GPd.plot1(g1, [-1], [1])
# %%


def sqexpprior(xx):
    p = sps.norm.pdf(((xx[0])-0.)/2.) * sps.norm.pdf(((xx[1])-0.)/2.)
    return sp.log(p)

e = EntropyPredict.EntPredictor([Xo, Yo, So, Do], [-1], [1], kfGen, sqexpprior)

s = 0.001


# %%
for k in xrange(30):
    n = 51
    Xi = sp.linspace(-1, 1, n)
    Si = sp.logspace(-4, 0, 10)
    a = e.showEntGraph(Xi, [s])
    xmin = e.searchAtS([-1], [1], [s])
    a.plot(xmin, [0], 'go')
    a.plot(xmin, f(xmin[0]), 'ro')
    xIn = xmin[0]
    yIn = f(xIn)+sp.random.normal(scale=s)

    Xo = sp.vstack([Xo, xIn])
    Yo = sp.vstack([Yo, yIn])
    So = sp.vstack([So, s])
    Do.append([sp.NaN])

    e = EntropyPredict.EntPredictor([Xo, Yo, So, Do], [-1], [1], kfGen, sqexpprior)
    plt.show()

# %%


e = EntropyPredict.EntPredictor([Xo, Yo, So, Do], [-1], [1], kfGen, sqexpprior)
e.drawHypSamples(100,plot=True)
e.initPredictor()
e.showEntGraph(Xi,[0.01])
# %%
n = 251
Xi = sp.linspace(-1, 1, n)
H = e.predictGraph(Xi,[0.01])
print H
# %%
for j in e.Pred[1]:
    print j[0].infer_full(sp.matrix(0.439999),[[sp.NaN]])

#GPep.GPcore.infer_full()
