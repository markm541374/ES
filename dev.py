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

hyptrue = [10., 0.20]
kfGen = GPep.gen_sqexp_k_d
kf = GPep.gen_sqexp_k_d(hyptrue)
upper = 1
lower = -1
ff = fgen1d(lower, upper, 500, kf)

def sqexpprior(xx):
    p = sps.norm.pdf(((xx[0])-0.)/(2*2.)) * sps.norm.pdf(((xx[1])-0.)/(2*2.))
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
s2n = 1.
N = 15
Uset = sp.linspace(1,N,11)
Sset = map(lambda u:(N-u)*s2n+10**-10, Uset)
plt.plot(Uset,Sset)


# %%

n_init = 5
sconst=0.001
x = sp.random.uniform(-1, 1, n_init)
y = map(f, x)+sp.random.normal(scale=0.01, size=n_init)
Xo = sp.matrix(x).T
Yo = sp.matrix(y).T
So = sp.matrix([[sconst]]*n_init)
# So = sp.matrix([[Sset[-1]]]*n_init)
Do = [[sp.NaN]]*n_init
g1 = GPd.GPcore(Xo, Yo, So, Do, kf)
a = GPd.plot1(g1, [-1], [1])



# %%
reload(EntropyPredict)
o = EntropyPredict.Optimizer(f,kfGen,sqexpprior,[-1],[1])
o.initspecobs(Xo,Yo,So,Do)
o.HYPsamplen=25
##o.showinov(0.000001, obstype=[0])
##o.searchats(10**-14, method='Ent')
#
#
## %%
#reload(EntropyPredict)
#e = EntropyPredict.EntPredictor([Xo,Yo,So,Do],[-1],[1],kfGen,sqexpprior)
#e.HYPsamplen=22
#e.ENTsearchn=160
##e.searchAtS([-1],[1],[0.0001],)
#n = 140
#Xi = sp.linspace(-1, 1, n)
##a = e.showEntGraph(Xi, [0.1,0.0001,0.0000000001], plot='basic')
#
#a = e.showEntGraph(Xi, sp.logspace(2,-10,35), U=sp.logspace(0,2,35),plot='basic')
#[xmin,j] = e.searchOverS([-1],[1], sp.logspace(2,-10,25), sp.logspace(0,8,25))

# %%
mt='EI'
reload(EntropyPredict)
for i in xrange(15):
    # o.searchOvers(sp.logspace(1,-7,30),sp.logspace(0,8,30), method='Ent')
    #o.searchOvers(Sset,Uset, method=mt)
    o.searchats([sconst], method=mt)

# %%

optresult = dict()
optresult['X'] = o.Xo
optresult['Y'] = o.Yo
optresult['S'] = o.So
optresult['D'] = o.Do
optresult['U'] = o.Uo
optresult['minest'] = o.result
optresult['xmintrue'] = xmintrue
optresult['ymintrue'] = miny
optresult['ninit'] = n_init
optresult['method'] = mt
optresult['time'] = time.time()
import pickle
object = optresult
file_pi = open('res/test4.obj', 'wb') 
pickle.dump(object, file_pi)

# %%

r = o.result


ymine=[]
xmine=[]
uacc=0
U=[]
for i in xrange(max(len(r),len(r2))):
    try:
        xmine.append(abs(r[i][0]-xmintrue))
        ymine.append(abs(r[i][1]-miny))
        uacc+=o.Uo[i]
        U.append(uacc)
    except:
        pass
    


plt.loglog(U,xmine,'b')

plt.figure()
plt.loglog(U,ymine,'b')

# %%
for j in xrange(10):
    f = ff.genfun()

    plt.plot(sp.linspace(-1, 1, 100), map(f, sp.linspace(-1, 1, 100)))
    ee = lambda x, y: (f(x), 0)
    [xmintrue, miny, ierror] = DIRECT.solve(ee, lower, upper, user_data=[], algmethod=1, maxf=4000)
    xmintrue = xmintrue[0]
    print miny
    plt.plot(xmintrue, miny, 'rx')
    
    n_init = 5
    sconst=0.001
    x = sp.random.uniform(-1, 1, n_init)
    y = map(f, x)+sp.random.normal(scale=0.01, size=n_init)
    Xo = sp.matrix(x).T
    Yo = sp.matrix(y).T
    So = sp.matrix([[sconst]]*n_init)
    # So = sp.matrix([[Sset[-1]]]*n_init)
    Do = [[sp.NaN]]*n_init
    g1 = GPd.GPcore(Xo, Yo, So, Do, kf)
    a = GPd.plot1(g1, [-1], [1])
    
    o = EntropyPredict.Optimizer(f,kfGen,sqexpprior,[-1],[1])
    o.initspecobs(Xo,Yo,So,Do)
    o.HYPsamplen=25
    
    mt='EI'
    try:
        for i in xrange(20):
    # o.searchOvers(sp.logspace(1,-7,30),sp.logspace(0,8,30), method='Ent')
    #o.searchOvers(Sset,Uset, method=mt)
            o.searchats([sconst], method=mt)
    except:
        pass
    
    optresult = dict()
    optresult['X'] = o.Xo
    optresult['Y'] = o.Yo
    optresult['S'] = o.So
    optresult['D'] = o.Do
    optresult['U'] = o.Uo
    optresult['minest'] = o.result
    optresult['xmintrue'] = xmintrue
    optresult['ymintrue'] = miny
    optresult['ninit'] = n_init
    optresult['method'] = mt
    optresult['time'] = time.time()
    import pickle
    object = optresult
    fname = 'res/N0/test'+str(j)+'.obj'
    file_pi = open(fname, 'wb') 
    pickle.dump(object, file_pi)
# %%