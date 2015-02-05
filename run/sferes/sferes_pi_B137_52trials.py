#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from optparse import OptionParser                                                   
sys.path.append("../../src")
from Wrap import TYMaze
from Models import PI
import cPickle as pickle

parameters = {'beta': 100.0,
 'eta': 0.109,
 'gamma': 0.1}
# print parameters
model = PI(parameters)
opt = TYMaze(model)
with open('../../data/B137_52trials.pickle') as f:
	data = pickle.load(f)
llh = opt.sferes(data)
# print llh