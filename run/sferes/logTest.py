#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from optparse import OptionParser                                                   
sys.path.append("../../src")
from Wrap import TYMaze
from Models import *
import cPickle as pickle
from matplotlib import *
from pylab import *

parameters0 = {'length':0,'beta':0.34951,'gamma':0.99999999,'eta':0.437553999996}

model = VMWM(parameters0)
opt = TYMaze(model)
with open('/home/viejo/Fos/data/B28_52trials.pickle') as f:
	data = pickle.load(f)
llh = opt.sferes(data)
print llh
log0 = opt.cumloglike

parameters3 = {'length':3,'beta':0.761214,'gamma':0.99999999,'eta':0.868842999991}

model = VMWM(parameters3)
opt = TYMaze(model)
llh = opt.sferes(data)
print llh
log3 = opt.cumloglike

figure()
plot(log0, label = 'n=0')
plot(log3, label = 'n=3')
legend()
show()