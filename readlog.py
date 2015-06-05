import scipy as sp
from matplotlib import pyplot as plt
import scipy.stats as sps
import scipy.linalg as spl
from scipy.stats import norm as norm
import GPep
import GPd
import GPset
# from scipy.optimize import minimize as mnz
import DIRECT
from functools import partial
# import time
from multiprocessing import pool
from tools import *
import sys
import traceback
import dill as pickle
import EntropyPredict

class OptEval():
    def __init__(self,fname):
        self.O = EntropyPredict.restartOpt(fname)
        return
        
    def xerr(self):
        xerr=[]
        xmintrue = sp.matrix(self.O.para['xmintrue'])
        for s in self.O.states[1:]:
            xminest = sp.matrix(s['xminIR'])
            xerr.append(spl.norm(xmintrue-xminest, ord='fro'))
        return xerr