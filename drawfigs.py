#! /usr/bin/env python


# -*- coding: utf-8 -*-
"""
Created on Sat Jun  6 14:46:46 2015
Drawfigs.py
@author: mark
"""

import os
import sys
import argparse
import time
import pprint
from tools import *
import GPep
import GPd
import DIRECT
import EntropyPredict
import logging
import readlog
from matplotlib import pyplot as plt
#commandline input parser

parser=argparse.ArgumentParser(prog='drawfigs')
parser.add_argument('-p','--para', nargs='?', default='default')
parser.add_argument('-n','--name', nargs='?', default='default')
args = parser.parse_args()

sys.path.append('plots')
para = __import__(args.para)

#create directory for results
if args.name=='default':
    rpath = 'plots/default'
    from shutil import rmtree
    if os.path.exists(rpath):
        rmtree(rpath)
    os.mkdir(rpath)
else:
    rpath = os.path.join('plots',args.name)
    if os.path.exists(rpath):
        i=1
        rpath=rpath+str(i)
        while os.path.exists(rpath):
            i+=1
            rpath=rpath[:-1]+str(i)
        
    os.mkdir(rpath)
    pass

dpaths=[os.path.join('results',p) for p in para.datasets]
data=[]
for dp in dpaths:
    data.append([readlog.OptEval(os.path.join(dp,f)) for f in os.listdir(dp) if f[-4:]=='.obj'])

for plot in para.plots:
    if plot['name']=='xerr':
        ydata=[[r.xerr() for r in d] for d in data]
        xdata=[[range(len(x)) for x in xs] for xs in ydata]
        [f,a] = readlog.plotany(xdata,ydata,plot)
        f.savefig(os.path.join(rpath,'xerr.png'))