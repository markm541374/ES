# -*- coding: utf-8 -*-
"""
Created on Sat Jun  6 14:51:25 2015

@author: mark
"""

#ploting parameters

datasets=['default']
plots=[]

#xmintrueplot
p=dict()
p['name']='xerr'
p['colorcodes']=['lightblue']
p['xscale']='linear'
p['yscale']='log'
p['title']='xmin_error vs. steps'
p['extras']=['median']
plots.append(p)

#IR plot
p=dict()
p['name']='IR'
p['colorcodes']=['lightblue']
p['xscale']='linear'
p['yscale']='log'
p['title']='IR vs. steps'
p['extras']=['median']
plots.append(p)

#time plot
p=dict()
p['name']='times'
p['colorcodes']=['lightblue']
p['xscale']='linear'
p['yscale']='linear'
p['title']='time vs. steps'
p['extras']=['mean']
plots.append(p)