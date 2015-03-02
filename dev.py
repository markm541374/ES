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

# %%

hyptrue = [1., 0.2]
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

n_init = 21
x = sp.random.uniform(-1, 1, n_init)
y = map(f, x)+sp.random.normal(scale=0.1, size=n_init)
Xo = sp.matrix(x).T
Yo = sp.matrix(y).T
So = sp.matrix([[0.01]]*n_init)
Do = [[sp.NaN]]*n_init
g1 = GPd.GPcore(Xo, Yo, So, Do, kf)
a = GPd.plot1(g1, [-1], [1])

# %%

kfGen=GPep.gen_sqexp_k_d
def squaresllk(xx, D):
    hyp = [10**(i) for i in xx]
    # print
    # print xx
    # print hyp
    kf = kfGen(hyp)
    g = GPd.GPcore(D[0], D[1], D[2], D[3], kf)
    lk = g.llk()
    return lk


def sqexpprior(xx):
    p = sps.norm.pdf(((xx[0])-0.)/2.) * sps.norm.pdf(((xx[1])-0.)/2.)
    return sp.log(p)


dist = dist_obj(squaresllk, [Xo, Yo, So, Do], sqexpprior)
# %%

n = 40
M = sp.zeros(n*n)
A = sp.linspace(-4, 4, n)
L = sp.linspace(-4, 4, n)

for y in xrange(40):
    for x in xrange(40):
        a = A[x]
        l = L[y]
        try:
            M[n*y+x] = dist.loglike([a, l])
        except:
            print a,l

Mgrid = sp.vstack(sp.split(M, n))
plt.contour(A, L, Mgrid, 50)
plt.xlabel('log10 outscale')
plt.ylabel('log10 lengthscale')
# %%

ee = lambda hyp, y: (-dist.loglike(hyp), 0)
# find the map hyperparameters
[xmintrue, miny, ierror] = DIRECT.solve(ee, [-4, -4], [4, 4], user_data=[], algmethod=1, maxf=500)
print [10**i for i in xmintrue]

# %%

w_0 = xmintrue

n = 100
sigma = 0.1*np.ones(2)
samples = slice_sample(dist, w_0, iters=n, sigma=sigma)

plt.figure()
plt.xlabel('outputscale')
plt.ylabel('lengthscale')
kfPass = []

for i in xrange(len(samples[0])):
    plt.loglog(10**samples[0][i], 10**samples[1][i], 'rx')
    kfPass.append(kfGen([10**samples[0][i], 10**samples[1][i]]))

# %%

GG = [[], []]
i = 0
for i, k in enumerate(kfPass):
    i += 1
    if i % 10 == 0:
        print i
    try:
        g = EntropyPredict.makedraws([Xo, Yo, So, Do], k, nd=1)
        GG[0].append(g[0][0])
        GG[1].append(g[1][0])
    except:
        kfPass.pop(i)
        print "removed a kf"

# %%

n = 201
Xi = sp.linspace(-1, 1, n)
# Xi=[0,12]
Si = [0.01]

# import sys
# tmp = sys.stderr
# sys.stderr = open('/dev/null','w')

H = sp.zeros(n)
for i in xrange(n):
    if i % 25 == 0:
        print i
    H[i] = EntropyPredict.inferHmulti(GG, sp.matrix(Xi[i]).T, Si)
# sys.stderr=tmp
# %%


kftrue = GPep.gen_sqexp_k_d(hyptrue)
g1 = GPd.GPcore(Xo, Yo, So, Do, kftrue)

# Xi=[0,12]
Si = [0.01]
GG0 = [[], []]
for i in xrange(n):
    g = EntropyPredict.makedraws([Xo, Yo, So, Do], kfGen(hyptrue), nd=1)
    GG0[0].append(g[0][0])
    GG0[1].append(g[1][0])
H0 = sp.zeros(n)
for i in xrange(n):
    if i % 25 == 0:
        print i
    H0[i] = EntropyPredict.inferHmulti(GG0, sp.matrix(Xi[i]).T, Si)

# %%

plt.figure(figsize=(20, 20))
a = GPd.plot1(g1, [-1], [1])
a.plot(Xi, H, 'r')
a.plot(Xi, H0, 'g')