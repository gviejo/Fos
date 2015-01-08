#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

load and and plot multi objective results from Sferes 2 optimisation 


Copyright (c) 2013 Guillaume VIEJO. All rights reserved.
"""

import sys
import os
import cPickle as pickle
from optparse import OptionParser
import numpy as np
sys.path.append("../../src")
from matplotlib import *
from pylab import *
from Sferes import pareto
from itertools import *

# -----------------------------------
# ARGUMENT MANAGER
# -----------------------------------
if not sys.argv[1:]:
    sys.stdout.write("Sorry: you must specify at least 1 argument")
    sys.stdout.write("More help avalaible with -h or --help option")
    sys.exit(0)
parser = OptionParser()
parser.add_option("-i", "--input", action="store", help="The name of the directory to load", default=False)
parser.add_option("-m", "--model", action="store", help="The name of the model to test \n If none is provided, all files are loaded", default=False)
parser.add_option("-o", "--output", action="store", help="The output file of best parameters to test", default=False)
(options, args) = parser.parse_args()
# -----------------------------------

# -----------------------------------
# LOADING DATA
# -----------------------------------
front = pareto(options.input)
front.extractFrontLimits()
front.plotEvolution()

front.extractBestParameters()

with open("parameters.pickle", 'wb') as f:
	pickle.dump(front.best, f)
with open("parameters.txt", 'w') as f:
	m = 'VMWM'
	for s in front.best['VMWM'].keys():
		line="mouse="+s+"\t"+" \t ".join([k+"="+str(front.best[m][s][k]) for k in ['beta','eta','gamma']])+"\t\tloglikelihood = "+str(front.best_log[m][s])+"\n"		
		f.write(line)

# all_mice = np.array([i.split(".")[0] for i in os.listdir("../../data/")])
# few_mice = np.array(['B28','B61','B137','B155','B163','B166','B150','B62','B84','B139','B154','B74','B86','B152'])
# list_mice = np.array([i for i in all_mice if i.split("_")[0] in few_mice])


sys.exit()

fig = figure()
m = 'VMWM'
for i in xrange(len(front.p_order[m])):
	ax = fig.add_subplot(3,1,i+1)
	p = front.p_order[m][i]
	for j in xrange(len(front.data[m].keys())):
		s = front.data[m].keys()[j]
		ax.plot(j, front.best[m][s][p], 'o')
		ax.set_ylim(front.bounds[m][p][0], front.bounds[m][p][1])
		ax.set_title(p)

# show()
