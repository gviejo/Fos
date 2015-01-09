#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from optparse import OptionParser                                                   
sys.path.append("../../src/")
from Wrap import TYMaze
from Models import VMWM
import cPickle as pickle
# parser = OptionParser()                                                             
# parser.add_option("-b", "--beta", action="store", type = 'float')
# parser.add_option("-g", "--gamma", action="store", type = 'float')
# parser.add_option("-e", "--eta", action="store", type = 'float')
# parser.add_option("-n", "--length", action="store", type = 'float')
# (options, args) = parser.parse_args()
with open("../../data/B139_68trials.pickle") as f:
	data = pickle.load(f)
# model = VMWM()
# parameters = vars(options)
# for p in parameters.iterkeys():                                                                                                                   
#     if parameters[p] is not None:
#         parameters[p] = model.bounds[p][0]+parameters[p]*(model.bounds[p][1]-model.bounds[p][0])
# parameters['length'] = 3
parameters = {'beta':0.124936,
			'eta':0.722747999993,
			'gamma':0.0,
			'length':0}			
model = VMWM(parameters)
opt = TYMaze(model)
llh = opt.sferes(data)
print llh