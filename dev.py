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
sns.set_palette("Reds")


C=readlog.OptEval('results/default/trace1.obj',lastinvalid=True).O
print len(C.states)
C.gotostate(len(C.states)-1)

try:
    C.setupEP()
except:
    print '!!!'

C.EP.searchMLEHYP()
C.EP.drawHYPsamples()
s=C.EP.initFBInfer()
C.EP.FBInfer.drawmins(150,[2,'method2'])