#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys, os
sys.path.append("../../src/")
import numpy as np
from Wrap import TYMaze
from Models import VMWM
import cPickle as pickle
from optparse import OptionParser

# -----------------------------------
# ARGUMENT MANAGER
# -----------------------------------
if not sys.argv[1:]:
    sys.stdout.write("Sorry: you must specify at least 1 argument")
    sys.stdout.write("More help avalaible with -h or --help option")
    sys.exit(0)
parser = OptionParser()
parser.add_option("-i", "--input", action="store", help="The name of the directory to load", default=False)
(options, args) = parser.parse_args()

with open("../latency.pickle", 'rb') as handle:
  latency = pickle.load(handle)

with open(options.input, 'rb') as handle:
	parameters = pickle.load(handle)
m = 'VMWM'
nb_exp = 50
wrap = TYMaze(VMWM())
data = {}
for s in parameters[m].keys():
	print m, s
	if 'length' not in parameters[m][s].keys():
		print "hello"
		parameters[m][s].update({'length':3})
	data[s] = wrap.test(parameters[m][s], nb_exp, int(s.split("_")[1]))
	# wrap.plot(data[s], s, parameters[m][s], latency[s.split("_")[0]], "../test/"+options.input.split("_")[0]+"_"+s+".pdf")

wrap.plotall(data, latency, "../test/"+options.input.split("_")[0]+"_all_test.pdf")

