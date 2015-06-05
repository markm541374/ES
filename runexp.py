#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 14:26:07 2015
Runexp
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

#commandline input parser
parser=argparse.ArgumentParser(prog='runGPGO')
parser.add_argument('-p','--paras', nargs='?', default='default')
parser.add_argument('-n','--name', nargs='?', default='default')
args = parser.parse_args()

sys.path.append('paras')
paras = __import__(args.paras)

#create directory for results
if args.name=='default':
    rpath = 'results/default'
    from shutil import rmtree
    if os.path.exists(rpath):
        rmtree(rpath)
    os.mkdir(rpath)
else:
    rpath = os.path.join('results',args.name)
    if os.path.exists(rpath):
        i=1
        rpath=rpath+str(i)
        while os.path.exists(rpath):
            i+=1
            rpath=rpath[:-1]+str(i)
        
    os.mkdir(rpath)
    pass

#print header to logfile
logf = open(os.path.join(rpath,'log.txt'),'w')
logf.write(parser.prog+ ' started at '+time.strftime('%H:%M:%S on %a %-d %B %Y') +' on ' +os.uname()[1]+'\n\n')

logf.write('INPUT from '+args.paras+':\n\nOptimisation parameters:\n')
pprint.pprint(paras.optpara,logf)
logf.write('\nObjective function parameters:\n')
pprint.pprint(paras.objf,logf)
logf.write('\n')

#make the generator function ofor the objective function

if paras.objf['type']=='drawfromcov':
    if paras.objf['covgen']=='sqexp':
        objkfGen = GPep.gen_sqexp_k_d
    hyptrue = paras.objf['hyp']
    kftrue = objkfGen(hyptrue)
upper = paras.objf['upper']
lower = paras.objf['lower']
functiongenerator = fgen1d(lower[0], upper[0], 300, kftrue)

#init kfgen and prior for covariance
if paras.optpara['covtype']=='sqexp':
    optkfGen = GPep.gen_sqexp_k_d
    optkfprior = genSqExpPrior([[paras.optpara['OSLM'],paras.optpara['OSLV']],[paras.optpara['I1LM'],paras.optpara['I1LV']]])

for i in xrange(paras.runs['nopts']):
    #draw an objective function
    f=functiongenerator.genfun()
    ee = lambda x, y: (f(x), 0)
    [xmintrue, miny, ierror] = DIRECT.solve(ee, lower, upper, user_data=[], algmethod=1, maxf=4000, logfilename='/dev/null')
    xmintrue = xmintrue[0]
    paras.optpara['xmintrue']=xmintrue
    paras.optpara['ymintrue']=miny
    O = EntropyPredict.Optimizer(f,optkfGen, optkfprior, lower, upper, paras.optpara)
    if paras.optpara['inittype']=='rand':
        O.initrandobs(paras.optpara['nrand'],paras.optpara['fixs'])
    O.runopt(paras.runs['nsteps'])
    O.savestate(os.path.join(rpath,'trace'+str(i)+'.obj'))
        
    
    
logf.write(parser.prog+ ' exited at '+time.strftime('%H:%M:%S on %a %-d %B %Y'))
logf.close()
