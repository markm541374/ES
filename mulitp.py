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
import readlog
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
para=dict()

para['xmintrue'] = xmintrue
para['ymintrue'] = miny

para['nHYPsamples']=8
para['HYPsearchLow'] = [-2, -2]
para['HYPsearchHigh'] = [2, 2]
para['HYPMLEsearchn'] = 800
para['HYPsamSigma'] = 0.05
para['HYPsamBurn'] = 16
para['ENTnsam'] = 100
para['ENTzeroprecision'] = 10**-6
para['ENTsearchn'] = 500
para['IRsearchn'] = 500
#para['searchmethod']='fixs'
#para['fixs'] = 0.0001
para['searchmethod']='EIMLE'
para['fixs'] = 0.0001

para['obstype'] = [sp.NaN]
#para = [nHYPsam, HYPsearchLow, HYPsearchHigh, HYPMLEsearchn, HYPsamSigma, HYPsamBurn, ENTnsam, ENTzeroprecision, ENTsearchn]
# %%

reload(EntropyPredict2)
reload(GPset)

O = EntropyPredict2.Optimizer(f,kfGen, kfprior, lower, upper, para)
O.initrandobs(5,para['fixs'])
#O.setupEP()
#O.plotstate()
# %%
for i in xrange(15):
    O.runopt(1)
    O.plotstate()
    plt.show()
# %%
O.savestate()

# %%
E = readlog.OptEval('states.obj')

# %%

plt.semilogy(E.xerr())