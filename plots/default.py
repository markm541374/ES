# -*- coding: utf-8 -*-
"""
Created on Sat Jun  6 14:51:25 2015

@author: mark
"""

#ploting parameters

datasets=['EIMLE_20_401','EIFB_20_401','ENFS_20_40']
plots=[]

#xmintrueplot
p=dict()
p['name']='xerr'
p['dsets']=[True,True,True]
p['colorcodes']=['lightblue','lightgreen','lightcoral']
p['xscale']='log'
p['yscale']='log'
p['title']='xmin_error vs. steps'
ex1=[dict()]
ex1[0]['name']='median'
ex1[0]['colorcode']='b'
ex2=[dict()]
ex2[0]['name']='median'
ex2[0]['colorcode']='g'
ex3=[dict()]
ex3[0]['name']='median'
ex3[0]['colorcode']='r'
p['extras']=[ex1,ex2,ex3]
plots.append(p)

#IR plot
p=dict()
p['name']='IR'
p['dsets']=[True,True,True]
p['colorcodes']=['lightblue','lightgreen','lightcoral']
p['xscale']='log'
p['yscale']='log'
p['title']='IR vs. steps'

ex1=[dict()]
ex1[0]['name']='median'
ex1[0]['colorcode']='b'
ex2=[dict()]
ex2[0]['name']='median'
ex2[0]['colorcode']='g'
ex3=[dict()]
ex3[0]['name']='median'
ex3[0]['colorcode']='r'
p['extras']=[ex1,ex2,ex3]



plots.append(p)

#time plot
p=dict()
p['name']='times'
p['dsets']=[True,True,True]
p['colorcodes']=['lightblue','lightgreen','lightcoral']
p['xscale']='linear'
p['yscale']='linear'
p['title']='time vs. steps'

ex1=[dict()]
ex1[0]['name']='mean'
ex1[0]['colorcode']='b'
ex2=[dict()]
ex2[0]['name']='mean'
ex2[0]['colorcode']='g'
ex3=[dict()]
ex3[0]['name']='mean'
ex3[0]['colorcode']='r'
p['extras']=[ex1,ex2,ex3]
plots.append(p)

#MLEhyperparameters
p=dict()
p['name']='MLEhyp'
p['dsets']=[True,True,True]
p['colorcodes']=['lightblue','lightgreen','lightcoral']
p['xscale']='linear'
p['yscale']='log'
p['title']='MLYHYP vs. steps'
p['extras']=[[],[],[]]
plots.append(p)
