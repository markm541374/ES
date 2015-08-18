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
import GPset
import DIRECT
import EntropyPredict
import logging
import readlog
from matplotlib import pyplot as plt
import traceback
import imp
import seaborn as sns
#commandline input parser

parser=argparse.ArgumentParser(prog='drawfigs')
parser.add_argument('-p','--para', nargs='?', default='default')
parser.add_argument('-n','--name', nargs='?', default='default')
args = parser.parse_args()

sys.path.append('plots')
para = imp.load_source(args.para,'plots/'+args.para+'.py')

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
    try:
        data.append([readlog.OptEval(os.path.join(dp,f)) for f in os.listdir(dp) if f[-4:]=='.obj'])
    except:
        print traceback.format_exc()
        
        print data
        print dp
        print f
        raise MJMError('couldn\'t read data')
for plot in para.plots:
    f = plt.figure()
    
    for i,d in enumerate([d for i,d in enumerate(data) if plot['dsets'][i]]):
        if plot['name']=='xerr':
            a = f.add_subplot(111)
            ydata=sp.vstack([r.xerr() for r in d])
            xdata=range(len(ydata[0]))
            readlog.plotset(ydata,xdata,f,a,plot['colorcodes'][i],extras=plot['extras'][i])
        
        if plot['name']=='IR':        
            a = f.add_subplot(111)
            ydata=sp.vstack([r.yIR() for r in d])
            xdata=range(len(ydata[0]))
            readlog.plotset(ydata,xdata,f,a,plot['colorcodes'][i],extras=plot['extras'][i])
            
        if plot['name']=='region':        
            a = f.add_subplot(111)
            ydata=sp.vstack([r.region() for r in d])
            xdata=range(len(ydata[0]))
            readlog.plotset(ydata,xdata,f,a,plot['colorcodes'][i],extras=plot['extras'][i])
            
        if plot['name']=='global_hyp':        
            a = f.add_subplot(111)
            ydata=sp.vstack([r.global_hyp() for r in d])
            xdata=range(len(ydata[0]))
            readlog.plotset(ydata,xdata,f,a,plot['colorcodes'][i],extras=plot['extras'][i])
            
        if plot['name']=='d3f':        
            a = f.add_subplot(111)
            ydata=sp.vstack([r.d3f() for r in d])
            xdata=range(len(ydata[0]))
            readlog.plotset(ydata,xdata,f,a,plot['colorcodes'][i],extras=plot['extras'][i])
            
        if plot['name']=='dvdf':        
            a = f.add_subplot(111)
            ydata=sp.vstack([r.dvdf() for r in d])
            xdata=range(len(ydata[0]))
            readlog.plotset(ydata,xdata,f,a,plot['colorcodes'][i],extras=plot['extras'][i])
        
        if plot['name']=='aqu':        
            a = f.add_subplot(111)
            ydata=sp.vstack([r.aqu() for r in d])
            xdata=range(len(ydata[0]))
            readlog.plotset(ydata,xdata,f,a,plot['colorcodes'][i],extras=plot['extras'][i])
            
        if plot['name']=='schosen':        
            a = f.add_subplot(111)
            ydata=sp.vstack([r.s_chosen() for r in d])
            xdata=range(len(ydata[0]))
            readlog.plotset(ydata,xdata,f,a,plot['colorcodes'][i],extras=plot['extras'][i])
            
        
        if plot['name']=='times':        
            a = f.add_subplot(111)
            ydata=sp.vstack([r.steptimes() for r in d])
            xdata=range(len(ydata[0]))
            readlog.plotset(ydata,xdata,f,a,plot['colorcodes'][i],extras=plot['extras'][i])
            
        if plot['name']=='MLEhyp':
            alldata = [r.MLEhyp() for r in d]
            ydata=[]
            D=len(alldata[0][0])
            
            for j in xrange(D):
                ydata.append(sp.vstack([[h[j] for h in x] for x in alldata]))
            xdata=range(len(alldata[0]))
            for j in xrange(D):
                a = f.add_subplot(D,1,j+1)
                readlog.plotset(ydata[j],xdata,f,a,plot['colorcodes'][i],extras=plot['extras'][i])
                a.set_yscale(plot['yscale'])
                a.set_xscale(plot['xscale'])
                a.set_title(plot['title'])
            
        
        if plot['name']=='sprofdata':
            alldata = [r.sprofile() for r in d]
            x=alldata[0][0][0]
            D = len(alldata[0])
            n = len(alldata)
            
            f.set_figwidth(8)
            f.set_figheight(3*D)
            for i in xrange(D):
                a = f.add_subplot(D,1,i+1)
                for j in xrange(n):
                    a.plot(x,alldata[j][i][1],'g')
                    a.plot(x,alldata[j][i][2],'b')
                    a.plot(x,alldata[j][i][3],'r')
                try:
                    a.set_yscale(plot['yscale'])
                except:
                    a.set_yscale('linear')
                a.set_xscale(plot['xscale'])
    try:
        a.set_yscale(plot['yscale'])
    except:
        a.set_yscale('linear')
    a.set_xscale(plot['xscale'])
    a.set_title(plot['title'])
    f.savefig(os.path.join(rpath,plot['name']+'.png'))
