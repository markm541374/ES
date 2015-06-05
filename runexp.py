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

parser=argparse.ArgumentParser(prog='runGPGO')
parser.add_argument('-p','--paras', nargs='?', default='default')
parser.add_argument('-n','--name', nargs='?', default='default')
args = parser.parse_args()

sys.path.append('paras')
paras = __import__(args.paras)

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

logf = open(os.path.join(rpath,'log.txt'),'w')
logf.write(parser.prog+ ' started at '+time.strftime('%H:%M:%S on %a %-d %B %Y') +' on ' +os.uname()[1]+'\n\n')

logf.write('INPUT from '+args.paras+':\n')
pprint.pprint(paras.para,logf)
logf.write('\n')




logf.write(parser.prog+ ' exited at '+time.strftime('%H:%M:%S on %a %-d %B %Y'))
logf.close()
