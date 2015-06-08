# -*- coding: utf-8 -*-
"""
Created on Sat Jun  6 14:51:25 2015

@author: mark
"""

#ploting parameters

datasets=['EIMLE_20_40']
plots=[]

#xmintrueplot
p=dict()
p['name']='xerr'
p['colorcodes']=['lightblue']
p['xscale']='log'
p['yscale']='log'
p['title']='xmin_error vs. steps'
p['extras']=['median']
plots.append(p)

#IR plot
p=dict()
p['name']='IR'
p['colorcodes']=['lightblue']
p['xscale']='log'
p['yscale']='log'
p['title']='IR vs. steps'
p['extras']=['median']
plots.append(p)