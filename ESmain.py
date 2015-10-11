#! /usr/bin/env python
__author__ = 'markm'

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
import runexp
import traceback
import shutil

parser=argparse.ArgumentParser(prog='runGPGO')
parser.add_argument('-p','--paras', nargs='?', default='default3') #paras, only used for new dsets
parser.add_argument('-n','--name', nargs='?', default='default') #name for dset toplevel folder
parser.add_argument('-a','--action', nargs='?', default='new') # new/extend/widen
parser.add_argument('-e', '--extendnumber', nargs='?',default='1')
args = parser.parse_args()




if args.action is None:
    print 'No action specified'
    exit()



if args.action == 'new':
    sys.path.append('paras')
    paras = __import__(args.paras)

    #create directory for results
    if args.name=='default':
        rpath = 'results/default'
        from shutil import rmtree
        if os.path.exists(rpath):
            rmtree(rpath, ignore_errors=True)
        os.mkdir(rpath)
        os.mkdir(os.path.join(rpath,'cache'))
    else:
        rpath = os.path.join('results',args.name)
        if os.path.exists(rpath):
            i=1
            rpath=rpath+str(i)
            while os.path.exists(rpath):
                i+=1
                rpath=rpath[:-1]+str(i)

        os.mkdir(rpath)
        os.mkdir(os.path.join(rpath,'cache'))
        pass

    shutil.copy(os.path.join('paras',args.paras+'.py'),os.path.join(rpath,'paras.py'))
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
    for i in xrange(paras.runs['nopts']):
        runexp.startrun(rpath,paras,i,paras.runs['nsteps'])

elif args.action == 'extend':
    rpath = os.path.join('results',args.name)
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
    fns = [os.path.join(rpath,f) for f in os.listdir(rpath) if f[-4:]=='.obj']
    en = int(args.extendnumber)
    for i in xrange(len(fns)):
        runexp.extend(fns[i],en)
        h,t = os.path.split(fns[i])
        shutil.move(fns[i],os.path.join(h,'cache',t))


elif args.action == 'widen':
    rpath = os.path.join('results',args.name)
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
    #paras = __import__(os.path.join(rpath,args.paras))
    import imp
    paras = imp.load_source('paras', os.path.join(rpath,'paras.py'))
    fns = [os.path.join(rpath,f) for f in os.listdir(rpath) if f[-4:]=='.obj']
    existing = len(fns)
    nsteps = len(EntropyPredict.restartOpt(fns[0]).states)-1
    for i in xrange(existing ,existing+int(args.extendnumber) ):
        runexp.startrun(rpath,paras,i,nsteps)


else:
    print 'invalid action'
    exit()


logger.info(parser.prog+ ' exited at '+time.strftime('%H:%M:%S on %a %-d %B %Y'))
try:
    if not __IPYTHON__:
        exit()
except:
    pass
