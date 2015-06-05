# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 14:27:46 2015

@author: mark
"""

para=dict()



para['nHYPsamples']=8
para['HYPsearchLow'] = [-2, -2]
para['HYPsearchHigh'] = [2, 2]
para['HYPMLEsearchn'] = 800
para['HYPsamSigma'] = 0.05
para['HYPsamBurn'] = 12
para['ENTnsam'] = 100
para['ENTzeroprecision'] = 10**-6
para['ENTsearchn'] = 500
para['IRsearchn'] = 500

para['searchmethod']='EIMLE'
para['fixs'] = 0.0001
from scipy import NaN as nan
para['obstype'] = [nan]