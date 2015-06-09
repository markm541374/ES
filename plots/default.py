# -*- coding: utf-8 -*-
"""
Created on Sat Jun  6 14:51:25 2015

@author: mark
"""

#ploting parameters

datasets=['EIMLE_20_40','EIFB_20_401']
plots=[]

#xmintrueplot
p=dict()
p['name']='xerr'
p['dsets']=[True,True]
p['colorcodes']=['lightblue','lightgreen']
p['xscale']='linear'
p['yscale']='log'
p['title']='xmin_error vs. steps'
ex1=[dict()]
ex1[0]['name']='median'
ex1[0]['colorcode']='r'
ex2=[dict()]
ex2[0]['name']='median'
ex2[0]['colorcode']='g'

p['extras']=[ex1,ex2]
plots.append(p)

#IR plot
p=dict()
p['name']='IR'
p['dsets']=[True,True]
p['colorcodes']=['lightblue','lightgreen']
p['xscale']='linear'
p['yscale']='log'
p['title']='IR vs. steps'

ex1=[dict()]
ex1[0]['name']='median'
ex1[0]['colorcode']='r'
ex2=[dict()]
ex2[0]['name']='median'
ex2[0]['colorcode']='g'

p['extras']=[ex1,ex2]



plots.append(p)

#time plot
p=dict()
p['name']='times'
p['dsets']=[True,True]
p['colorcodes']=['lightblue','lightgreen']
p['xscale']='linear'
p['yscale']='linear'
p['title']='time vs. steps'

ex1=[dict()]
ex1[0]['name']='mean'
ex1[0]['colorcode']=['r']
ex2=[dict()]
ex2[0]['name']='mean'
ex2[0]['colorcode']=['g']

p['extras']=[ex1,ex2]
plots.append(p)
