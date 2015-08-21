# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 14:27:46 2015

@author: mark
"""
#optimisation run parameters
import scipy as sp

optpara=dict()
optpara['D']=3
optpara['nHYPsamples']=10
optpara['nHYPmargin']=1.8
optpara['HYPsearchLow'] = [-2, -2, -2, -2]
optpara['HYPsearchHigh'] = [2, 2, 2, 2]
optpara['HYPMLEsearchn'] = 1200
optpara['HYPsamSigma'] = 0.05
optpara['HYPsamBurn'] = 12
optpara['ENTnsam'] = 400
optpara['ENTsamQ'] = 'method2'
optpara['ENTzeroprecision'] = 10**-6
optpara['ENTsearchn'] = 2000
optpara['IRsearchn'] = 2000
optpara['searchmethod']=['EIMLE','EIFB','fixs','discretes'][3]
optpara['slist'] = sp.logspace(1,-6,20)
optpara['ulist'] = [1./(x**0.1) for x in optpara['slist']]
optpara['fixs'] = 0.001
from scipy import NaN as nan
optpara['obstype'] = [nan]
optpara['covtype'] = 'sqexp'
#outpuscallogmean
#optpara['OSLM']=0.
#outputscalelogvar
#optpara['OSLV']=2.**2
#optpara['I1LM']=0.
#optpara['I1LV']=2.**2
#prior over hyp, mean and var on logscale, output then axes
optpara['prior'] = [[0.,2.*2],[0., 2.**2], [0., 2.**2], [0., 2.**2]]
optpara['inittype']='rand'
optpara['nrand']=20
optpara['boundregion']=0.99



objf=dict()
objf['type']='drawfromcov'
objf['covgen']='sqexp'
objf['D']=optpara['D']
objf['hyp']=[1.,0.5, 0.4,0.6]
#always scale to [-1.1]^D
#objf['lower']=[-1.]
#objf['upper']=[1.]



runs=dict()
runs['nopts']=10
runs['nsteps']= 35
