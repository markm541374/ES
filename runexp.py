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
import logging
import traceback
#import cProfile, pstats, StringIO
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
        rmtree(rpath, ignore_errors=True)
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
logger=logging.getLogger()
logger.setLevel(logging.DEBUG)
db=logging.FileHandler(os.path.join(rpath,'debug.log'))
db.setLevel(logging.DEBUG)
inf=logging.FileHandler(os.path.join(rpath,'info.log'))
inf.setLevel(logging.INFO)
formatter=logging.Formatter('%(asctime)s:%(levelname)s:%(module)s:%(message)s')
db.setFormatter(formatter)
inf.setFormatter(formatter)
logger.addHandler(db)
logger.addHandler(inf)

header=''
header+=parser.prog+ ' started at '+time.strftime('%H:%M:%S on %a %-d %B %Y') +' on ' +os.uname()[1]+'\n\n'
header+='INPUT from '+args.paras+':\n\nOptimisation parameters:\n'
logger.info(header)
detail=''
detail+=pprint.pformat(paras.optpara)
detail+='\n\nObjective function parameters:\n'
detail+=pprint.pformat(paras.objf)
detail += '\n'
logger.debug(detail)
#make the generator function ofor the objective function
#pr = cProfile.Profile()
#pr.enable()


if paras.objf['type']=='drawfromcov':
    if paras.objf['covgen']=='sqexp':
        objkfGen = GPep.gen_sqexp_k_d
    hyptrue = paras.objf['hyp']
    kftrue = objkfGen(hyptrue)
upper = paras.objf['upper']
lower = paras.objf['lower']
functiongenerator = fgennd(1, 150, kftrue)

#init kfgen and prior for covariance
if paras.optpara['covtype']=='sqexp':
    optkfGen = GPep.gen_sqexp_k_d
    optkfprior = genSqExpPrior(paras.optpara['prior'])

for i in xrange(paras.runs['nopts']):
    logger.info('starting run '+str(i)+'\n')
    #draw an objective function
    xmintrue=lower[0]
    while min(xmintrue-lower[0], upper[0]-xmintrue) < 0.025*(upper[0]-lower[0]):
        f=functiongenerator.genfun()
        ee = lambda x, y: (f(x), 0)
        [xmintrue, miny, ierror] = DIRECT.solve(ee, lower, upper, user_data=[], algmethod=1, maxf=4000, logfilename='/dev/null')
        xmintrue = xmintrue[0]
        print 'truemin: '+str(xmintrue)
    paras.optpara['xmintrue']=xmintrue
    paras.optpara['ymintrue']=miny
    O = EntropyPredict.Optimizer(f,optkfGen, optkfprior, lower, upper, paras.optpara)
    if paras.optpara['inittype']=='rand':
        O.initrandobs(paras.optpara['nrand'],paras.optpara['fixs'])
    try:
        O.runopt(paras.runs['nsteps'])
    except:
        logger.error('optimisation did not complete cleanly\n'+traceback.format_exc())
    O.savestate(os.path.join(rpath,'trace'+str(i)+'.obj'))
        
#pr.disable()
#s = StringIO.StringIO()
#sortby = 'cumulative'
#ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
#ps.print_stats()
#print s.getvalue()

    
logger.info(parser.prog+ ' exited at '+time.strftime('%H:%M:%S on %a %-d %B %Y'))
try:
    if not __IPYTHON__:
        exit()
except:
    pass