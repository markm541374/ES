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
import readlog
#import cProfile, pstats, StringIO
#commandline input parser

def extend(fpath,nsteps):
    O = EntropyPredict.restartOpt(fpath)
    O.runopt(nsteps)
    O.savestate(os.path.join(os.path.split(fpath)[0],'trace'+str(O.states[0]['id'])+'_'+str(len(O.states)-1)+'.obj'))
    return

def startrun(rpath,paras,id,nsteps):
#print header to logfile
    logger=logging.getLogger()
    

    header=''
    header+='started at '+time.strftime('%H:%M:%S on %a %-d %B %Y') +' on ' +os.uname()[1]+'\n\n'
    header+='INPUT :\n\nOptimisation parameters:\n'
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

    D = paras.objf['D']
    upper = [1.]*D
    lower = [-1]*D
    if paras.objf['type']=='drawfromcov':
        if paras.objf['covgen']=='sqexp':
            objkfGen = GPep.gen_sqexp_k_d
        elif paras.objf['covgen']=='mat32':
            objkfGen = GPep.gen_mat32_k_d
        else:
            raise MJMError('bad covtype')
        hyptrue = paras.objf['hyp']
        kftrue = objkfGen(hyptrue)
        functiongenerator = fgennd(D, 50, kftrue)
    elif paras.objf['type']=='branin':
        functiongenerator = fgenbranin()
    elif paras.objf['type']=='rosen':
        functiongenerator = fgenrosen()

    

    #init kfgen and prior for covariance
    if paras.optpara['covtype']=='sqexp' or paras.optpara['covtype']=='mat32':
        optkfGen = GPep.gen_sqexp_k_d
        optkfprior = genSqExpPrior(paras.optpara['prior'])


    logger.info('starting run \n')
    #draw an objective function
    xmintrue=lower
    first = True
    while any([j>0.99 or j<-0.99 for j in xmintrue]) and (not paras.objf['type']=='branin' or paras.objf['type']=='rosen' or first):
        first = False
        f=functiongenerator.genfun()
        ee = lambda x, y: (f(x), 0)
        [xmintrue, miny, ierror] = DIRECT.solve(ee, lower, upper, user_data=[], algmethod=1, maxf=4000, logfilename='/dev/null')
        print 'truemin: '+str(xmintrue)
        logger.debug('truemin: '+str([xmintrue,miny]))
    paras.optpara['xmintrue']=xmintrue
    paras.optpara['ymintrue']=miny
    O = EntropyPredict.Optimizer(f,optkfGen, optkfprior, D, paras.optpara)
    if paras.optpara['inittype']=='rand':
        O.initrandobs(paras.optpara['nrand'],paras.optpara['fixs'])
    try:
        O.runopt(nsteps)
    except:
        logger.error('optimisation did not complete cleanly\n'+traceback.format_exc())
    O.states[0]['id'] = id
    O.savestate(os.path.join(rpath,'trace'+str(O.states[0]['id'])+'_'+str(len(O.states)-1)+'.obj'))
    del(O)
    return

