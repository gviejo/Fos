#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from optparse import OptionParser                                                   
sys.path.append("/home/viejo/Fos/src/")
from Wrap import TYMaze
from Models import Graph
import cPickle as pickle
parser = OptionParser()                                                             
parser.add_option("-b", "--beta", action="store", type = 'float')
parser.add_option("-g", "--gamma", action="store", type = 'float')
parser.add_option("-e", "--eta", action="store", type = 'float')
(options, args) = parser.parse_args()
parameters = vars(options)
# model = Graph()
# for p in parameters.iterkeys():                                                                                                              
#     if parameters[p] is not None:
#         parameters[p] = model.bounds[p][0]+parameters[p]*(model.bounds[p][1]-model.bounds[p][0])
parameters = {'beta': 10.405800000000001,
 'eta': 0.9,
 'gamma': 0.9}

model = Graph(parameters)
opt = TYMaze(model)
with open('/home/viejo/Fos/data/B137_52trials.pickle') as f:
	data = pickle.load(f)
llh = opt.sferes(data)
print llh