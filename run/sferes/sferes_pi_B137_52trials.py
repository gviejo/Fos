#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from optparse import OptionParser                                                   
sys.path.append("../../src")
from Wrap import TYMaze
from Models import PI
import cPickle as pickle

parameters = {'beta': 13.2636,
 'eta': 0.138274,
 'gamma': 0.183356}
   
# print parameters
model = PI(parameters)
opt = TYMaze(model)
with open('../../data/B163_52trials.pickle') as f:
	data = pickle.load(f)
llh = opt.sferes(data)
print llh