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
#from matplotlib import *
#from pylab import *
from Sferes import *
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
front = maxlikelihood(options.input)
# front.extractFrontLimits()
# front.plotEvolution()
front.extractBestLog()
front.write(options.input)


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
