# -*- coding: utf-8 -*-
"""
Created on Tue May  5 15:22:15 2015

@author: markm
"""

from multiprocessing import Process, Pipe, Event, active_children

import sys
import GPd
import GPep
import GPset
import time
import scipy as sp
from tools import *
import DIRECT
import EntropyPredict2
import scipy.stats as sps
%pylab inline

# %%
hyptrue = [1., 0.20]
kfGen = GPep.gen_sqexp_k_d
kftrue = GPep.gen_sqexp_k_d(hyptrue)
upper = [1]
lower = [-1]
ff = fgen1d(lower[0], upper[0], 500, kftrue)




f = ff.genfun()

plt.plot(sp.linspace(-1, 1, 100), map(f, sp.linspace(-1, 1, 100)))
ee = lambda x, y: (f(x), 0)
[xmintrue, miny, ierror] = DIRECT.solve(ee, lower, upper, user_data=[], algmethod=1, maxf=4000, logfilename='/dev/null')
xmintrue = xmintrue[0]
print miny

plt.plot(xmintrue, miny, 'rx')
# %%
n_init = 4
sconst=0.0001
x = sp.random.uniform(-1, 1, n_init)
y = map(f, x)+sp.random.normal(scale=0.01, size=n_init)
Xo = sp.matrix(x).T
Yo = sp.matrix(y).T
So = sp.matrix([[sconst]]*n_init)
# So = sp.matrix([[Sset[-1]]]*n_init)
Do = [[sp.NaN]]*n_init
g1 = GPd.GPcore(Xo, Yo, So, Do, kftrue)
a = GPd.plot1(g1, [-1], [1])


# %%

reload(EntropyPredict2)
reload(GPset)

PO=[]
for c in active_children():
    c.terminate()
#outputscaleLogmean
OSLM=0.
#outputscaleLogvar
OSLV=2.**2
#inputscalaLogmean
I1LM=0.
#inputscaleLogvar
I1LV=2.**2

kfprior = genSqExpPrior([[OSLM,OSLV],[I1LM,I1LV]])

nHYPsam=16
HYPsearchLow = [-3, -3]
HYPsearchHigh = [3, 3]
HYPMLEsearchn = 800
HYPsamSigma = 0.05
HYPsamBurn = 8
ENTnsam = 100
ENTzeroprecision = 10**-6
para = [nHYPsam, HYPsearchLow, HYPsearchHigh, HYPMLEsearchn, HYPsamSigma, HYPsamBurn, ENTnsam, ENTzeroprecision]
PO=EntropyPredict2.EntPredictor([Xo,Yo,So,Do], lower, upper, kfGen, kfprior, para )
# %%
PO.setupEP()
#print PO.HYPsampleVals
[f0,a0] = PO.plotHYPsamples(d0=0, d1=1)
[f1,a1] = PO.plotFBpost()
[f2,a2] = PO.plotMLEpost()
[f3,a3] = PO.plotMinDraws()
[f4,as4] = PO.plotEPchanges()
# %%


# %%
X_s=sp.matrix([0])
D_s=[[sp.NaN]]
Xmcs=[]
Dmcs=[]
for i in xrange(PO.nHYPsamples):
    X = sp.vstack([X_s,PO.ENTmindraws[i][1]])
    Xmcs.append(X)
    Dmcs.append(D_s+[[sp.NaN]])
i2 = PO.EPInfer.infer_full_var(Xmcs,Dmcs)