# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 14:27:46 2015

@author: mark
"""
#optimisation run parameters
import scipy as sp

optpara=dict()
optpara['D']=1
optpara['nHYPsamples']=6
optpara['nHYPmargin']=1.8
optpara['HYPsearchLow'] = [-2, -2]
optpara['HYPsearchHigh'] = [2, 2]
optpara['HYPMLEsearchn'] = 800
optpara['HYPsamSigma'] = 0.05
optpara['HYPsamBurn'] = 12
optpara['ENTnsam'] = 400
optpara['ENTsamQ'] = 'method2'
optpara['ENTzeroprecision'] = 10**-6
optpara['ENTsearchn'] = 500
optpara['IRsearchn'] = 1000
optpara['searchmethod']='discretes'
optpara['slist'] = sp.logspace(-8,2,30)
optpara['ulist'] = sp.logspace(-4,-8,30)
optpara['fixs'] = 0.001
from scipy import NaN as nan
optpara['obstype'] = [nan]
optpara['covtype'] = 'sqexp'
#outpuscallogmean
optpara['OSLM']=0.
#outputscalelogvar
optpara['OSLV']=2.**2
optpara['I1LM']=0.
optpara['I1LV']=2.**2
optpara['inittype']='rand'
optpara['nrand']=5
optpara['boundregion']=0.995



objf=dict()
objf['type']='drawfromcov'
objf['covgen']='sqexp'
objf['D']=1
objf['hyp']=[1.,0.2]
objf['lower']=[-1.]
objf['upper']=[1.]



runs=dict()
runs['nopts']=1
runs['nsteps']=4