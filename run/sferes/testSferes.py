#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys, os
sys.path.append("../../src/")
import numpy as np
from Wrap import TYMaze
from Models import *
import cPickle as pickle
from Sferes import *
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

models = {'VMWM':VMWM(),
		'Graph':Graph(),
		'PI':PI()}

nb_exp = 5
data = {}
for m in parameters.keys():
	wrap = TYMaze(models[m])
	data[m] = {}
	for s in parameters[m].keys():
		print m, s
		if m == 'VMWM' and 'length' not in parameters[m][s].keys():			
			parameters[m][s].update({'length':3})	
		data[m][s] = wrap.test(parameters[m][s], nb_exp, int(s.split("_")[1]))	
		
		wrap.plot(data[m][s], s, parameters[m][s], latency[s.split("_")[0]], "../test/"+options.input.split("_")[0]+"_"+s+".pdf")


with open("data_tmp2", 'rb') as handle:
	data = pickle.load(handle)
wrap.plotall(data, latency, "../test/"+options.input.split("_")[0])

os.system("evince ../test/SFERES9_group_test_all_models.pdf")

