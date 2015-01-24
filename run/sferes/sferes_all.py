#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from optparse import OptionParser                                                   
sys.path.append("/home/viejo/Fos/src/")
from Wrap import TYMaze
from Models import VMWM
import cPickle as pickle
parser = OptionParser()                                                             
parser.add_option("-b", "--beta", action="store", type = 'float')
parser.add_option("-g", "--gamma", action="store", type = 'float')
parser.add_option("-e", "--eta", action="store", type = 'float')
parser.add_option("-n", "--length", action="store", type='float')
(options, args) = parser.parse_args()
parameters = vars(options)
model = VMWM()
for p in parameters.iterkeys():                                                                                                              
    if parameters[p] is not None:
        parameters[p] = model.bounds[p][0]+parameters[p]*(model.bounds[p][1]-model.bounds[p][0])
# parameters['length'] = int(parameters['length'])
parameters['length'] = 3
model = VMWM(parameters)
opt = TYMaze(model)

list_mice = ['B154_68trials',
 'B155_52trials',
 'B28_52trials',
 'B166_52trials',
 'B62_68trials',
 'B163_52trials',
 'B74_84trials',
 'B152_84trials',
 'B150_52trials',
 'B84_68trials',
 'B86_84trials',
 'B137_52trials',
 'B61_52trials',
 'B139_68trials']

llh = 0.0
for s in list_mice:
	with open('/home/viejo/Fos/data/'+s+'trials.pickle') as f:
		data = pickle.load(f)
	llh+=opt.sferes(data,1.0)
print llh