# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 14:27:46 2015

@author: mark
"""
#optimisation run parameters

optpara=dict()
optpara['D']=1
optpara['nHYPsamples']=8
optpara['HYPsearchLow'] = [-2, -2]
optpara['HYPsearchHigh'] = [2, 2]
optpara['HYPMLEsearchn'] = 800
optpara['HYPsamSigma'] = 0.05
optpara['HYPsamBurn'] = 12
optpara['ENTnsam'] = 100
optpara['ENTzeroprecision'] = 10**-6
optpara['ENTsearchn'] = 500
optpara['IRsearchn'] = 500
optpara['searchmethod']='EIMLE'
optpara['fixs'] = 0.0001
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




objf=dict()
objf['type']='drawfromcov'
objf['covgen']='sqexp'
objf['D']=1
objf['hyp']=[1.,0.2]
objf['lower']=[-1.]
objf['upper']=[1.]



runs=dict()
runs['nopts']=4
runs['nsteps']=12
