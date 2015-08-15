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

sys.path.append('paras')
paras = __import__('default3')

#print header to logfile
rpath = os.path.join('results','default')
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

if paras.objf['type']=='drawfromcov':
    if paras.objf['covgen']=='sqexp':
        objkfGen = GPep.gen_sqexp_k_d
    hyptrue = paras.objf['hyp']
    kftrue = objkfGen(hyptrue)
D = paras.objf['D']
upper = [1.]*D
lower = [-1]*D
functiongenerator = fgennd(D, 50, kftrue)

#init kfgen and prior for covariance
if paras.optpara['covtype']=='sqexp':
    optkfGen = GPep.gen_sqexp_k_d
    optkfprior = genSqExpPrior(paras.optpara['prior'])

#logger.info('starting run '+str(i)+'\n')
#draw an objective function
xmintrue=lower
while any([i>0.99 or i<-0.99 for i in xmintrue]):
    f=functiongenerator.genfun()
    ee = lambda x, y: (f(x), 0)
    [xmintrue, miny, ierror] = DIRECT.solve(ee, lower, upper, user_data=[], algmethod=1, maxf=1000, logfilename='/dev/null')
    print 'truemin: '+str(xmintrue)
paras.optpara['xmintrue']=xmintrue[0]
paras.optpara['ymintrue']=miny

e = EntropyPredict.Optimizer(f,optkfGen,optkfprior,D,paras.optpara)
e.initrandobs(paras.optpara['nrand'],paras.optpara['fixs'])

e.runopt(2)
del e
print traceback.print_exc()