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
from multiprocessing import Pool, Process
# -----------------------------------
# Functions 
# -----------------------------------
def pool_test(arg):	
	parameters = arg[0]
	nb_exp = arg[1]				
	data = {}
	model = PI()
	wrap = TYMaze(model)
	
	for s in parameters.keys():
		data[s] = wrap.test(parameters[s], nb_exp, int(s.split("_")[1]))

	return data


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

nb_exp = 100
data = {}
# for m in ['VMWM', 'Graph']:
for m in ['VMWM']:
	wrap = TYMaze(models[m])
	data[m] = {}
	for s in parameters[m].keys():
		print m, s
		if m == 'VMWM' and 'length' not in parameters[m][s].keys():
			parameters[m][s].update({'length':0})	
			print "Setting length to ", parameters[m][s]['length']			
			
		data[m][s] = wrap.test(parameters[m][s], nb_exp, int(s.split("_")[1]))	
		
		# wrap.plot(data[m][s], s, parameters[m][s], latency[s.split("_")[0]], "../test/"+options.input.split("_")[0]+"_"+s+".pdf")

# m = 'PI'
# wrap = TYMaze(models[m])
# data[m] = {}
# pool = Pool(6)
# tmp = np.array(parameters[m].keys()).reshape(6,5)
# args = [({s:parameters[m][s] for s in tmp[i]},nb_exp) for i in xrange(len(tmp))]

# tmp = pool.map(pool_test, args)
# for i in xrange(len(tmp)):
# 	for s in tmp[i].keys():
# 		data[m][s] = tmp[i][s]

# with open("data_tmp", 'rb') as handle:
# 	data = pickle.load(handle)
wrap.plotall(data, latency, "../test/"+options.input.split("_")[0])



# os.system("evince ../test/_group_test_all_models.pdf")

