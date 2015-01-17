#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from optparse import OptionParser                                                   
sys.path.append("../../src/")
from Wrap import TYMaze
import numpy as np
from Models import VMWM
import cPickle as pickle
from matplotlib import *
from pylab import *
# parser = OptionParser()                                                             
# parser.add_option("-b", "--beta", action="store", type = 'float')
# parser.add_option("-g", "--gamma", action="store", type = 'float')
# parser.add_option("-e", "--eta", action="store", type = 'float')
# parser.add_option("-n", "--length", action="store", type = 'float')
# (options, args) = parser.parse_args()
with open("../../data/B74_84trials.pickle") as f:
	data = pickle.load(f)
# with open("../latency.pickle", 'rb') as f:
# 	latency = pickle.load(f)
# model = VMWM()

# parameters = vars(options)
# for p in parameters.iterkeys():                                                                                                                   
#     if parameters[p] is not None:
#         parameters[p] = model.bounds[p][0]+parameters[p]*(model.bounds[p][1]-model.bounds[p][0])
# parameters['length'] = 3
epsilon = 0.0		
tmp = []
gamma = [0.01, 0.5, 0.99999999999]
beta = [0.1, 1.0, 4.40855, 8.0]
for b in beta :
	parameters = {'beta':b,
				'eta':0.0947603999991,
				'gamma':0.9999999999 ,
				'length':3}		
	model = VMWM(parameters)
	opt = TYMaze(model)
	llh = opt.sferes(data, epsilon)
	tmp.append(llh)
	# print llh

tmp = np.array(tmp)
print tmp.sum(1)
figure()

for i in xrange(len(tmp)):
	plot(tmp[i], 'o-', label=beta[i])
	legend()

show()