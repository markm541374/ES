__author__ = 'mark'

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

parser=argparse.ArgumentParser(prog='drawfigs')
#parser.add_argument('-p','--para', nargs='?', default='default')
parser.add_argument('-n','--name', nargs='?', default='default')
args = parser.parse_args()


rpath = os.path.join('results',args.name)

out = open(os.path.join('results',args.name,'results.txt'),'w+')

out.write('RESULTS\n\n')


data=[]


try:
    data=[readlog.OptEval(os.path.join(rpath,f)) for f in os.listdir(rpath) if f[-4:]=='.obj']
except:
    print traceback.format_exc()

    print data
    print dp
    print f
    raise MJMError('couldn\'t read data')

out.write('{0:d} optimisations\n\n'.format(len(data)))
#terminal xerr
try:
    xerr=[r.xerr()[-1] for r in data]
    out.write('Terminal xerr data\n\n')
    out.write(statstring(xerr))
except:
    out.write('no xerr data')
    #raise

#terminal IR
try:
    v=[r.yIR()[-1] for r in data]
    out.write('\n\nTerminal IR data\n\n')
    out.write(statstring(v))
except:
    out.write('XXXXX no IR data')
    #raise

#terminal S
try:
    v=[r.s_chosen()[-1] for r in data]
    out.write('\n\nTerminal s data\n\n')
    out.write(statstring(v))
except:
    out.write('XXXXX no s data')
    raise

#initial S
try:
    v=[r.s_chosen()[0] for r in data]
    out.write('\n\nInitial s data\n\n')
    out.write(statstring(v))
except:
    out.write('XXXXX no s data')
    raise

out.close()

