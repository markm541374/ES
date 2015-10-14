# -*- coding: utf-8 -*-
"""
Created on Sat Jun  6 14:51:25 2015

@author: mark
"""

#ploting parameters


datasets=['default']#,'m3EIMLE']
plots=[]

s2u = lambda x:1./(x**0.25)
#xmintrueplot
p=dict()
p['name']='xerr'
p['dsets']=[True,True,True]
p['colorcodes']=['lightblue','lightgreen','lightcoral']
p['xscale']='linear'
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


#xerru plot
p=dict()
p['name']='xerru'
p['dsets']=[True,True,True]
p['colorcodes']=['lightblue','lightgreen','lightcoral']
p['xscale']='linear'
p['yscale']='log'
p['title']='xerr vs. uacc'
x1=[dict()]
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

#dxerru plot
p=dict()
p['name']='dxerru'
p['dsets']=[True,True,True]
p['colorcodes']=['lightblue','lightgreen','lightcoral']
p['xscale']='linear'
p['yscale']='linear'
p['title']='dxerr vs. u'
p['extras']=[[],[],[]]
#plots.append(p)


#s plot
p=dict()
p['name']='schosen'
p['dsets']=[True,True,True]
p['colorcodes']=['lightblue','lightgreen','lightcoral']
p['xscale']='linear'
p['yscale']='log'
p['title']='s vs. steps'

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

#u plot
p=dict()
p['name']='uchosen'
p['dsets']=[True,True,True]
p['colorcodes']=['lightblue','lightgreen','lightcoral']
p['xscale']='linear'
p['yscale']='log'
p['title']='u vs. steps'

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

#sprofdata
p=dict()
p['name']='sprofdata'
p['dsets']=[True,True,True]
p['colorcodes']=['lightblue','lightgreen','lightcoral']
p['xscale']='log'
p['yscale']='log'
p['title']='H vs s'
p['extras']=[[],[],[]]
#plots.append(p)

#region
p=dict()
p['name']='region'
p['dsets']=[True,True,True]
p['colorcodes']=['lightblue','lightgreen','lightcoral']
p['xscale']='log'
p['yscale']='log'
p['title']='region vs. steps'
p['extras']=[[],[],[]]
#plots.append(p)

#d3f
p=dict()
p['name']='d3f'
p['dsets']=[True,True,True]
p['colorcodes']=['lightblue','lightgreen','lightcoral']
p['xscale']='linear'
p['yscale']='log'
p['title']='d3f vs. steps'
p['extras']=[[],[],[]]
#plots.append(p)

#aqu
p=dict()
p['name']='aqu'
p['dsets']=[True,True,True]
p['colorcodes']=['lightblue','lightgreen','lightcoral']
p['xscale']='linear'
p['yscale']='linear'
p['title']='aqu'
p['extras']=[[],[],[]]
plots.append(p)

#dvdf
p=dict()
p['name']='dvdf'
p['dsets']=[True,True,True]
p['colorcodes']=['lightblue','lightgreen','lightcoral']
p['xscale']='linear'
p['yscale']='log'
p['title']='dvdf vs. steps'
p['extras']=[[],[],[]]
#plots.append(p)

#global hyp
p=dict()
p['name']='global_hyp'
p['dsets']=[True,True,True]
p['colorcodes']=['lightblue','lightgreen','lightcoral']
p['xscale']='linear'
p['yscale']='linear'
p['title']='global_hyp vs. steps'
p['extras']=[[],[],[]]
#plots.append(p)

#IR plot
p=dict()
p['name']='IR'
p['dsets']=[True,True,True]
p['colorcodes']=['lightblue','lightgreen','lightcoral']
p['xscale']='linear'
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

#IRu plot
p=dict()
p['name']='IRu'
p['dsets']=[True,True,True]
p['colorcodes']=['lightblue','lightgreen','lightcoral']
p['xscale']='linear'
p['yscale']='log'
p['title']='IR vs. uacc'

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

