#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys, os
sys.path.append("../../src/")
import numpy as np
from manip import Manip
from matplotlib import *
from pylab import *
import cPickle as pickle

def test(name, parameters, mouse_latency):
	print "Testing "+mouse	
	nb_trials = int(mouse.split("_")[1][0:2])
	manip = Manip(nb_mice, nb_trials, "TYM","VMWM",parameters,"guidage")
	mean_time, std_time = manip.test()
	figure()
	plot(np.arange(nb_trials), mean_time, 'o-', label = 'Model')
	fill_between(np.arange(nb_trials), mean_time-(std_time/2.), mean_time+(std_time/2.), alpha = 0.5)
	plot(mouse_latency, 'o-', color = 'red', label = 'Mouse')
	desc = ", ".join([k+"="+str(np.round(parameters[k],2)) for k in parameters.keys()])
	ylabel("Latency")
	title(name+" | "+desc, fontsize = 11)
	grid()
	legend()
	savefig("sferes_norm/test_latency_"+name+".pdf")	
	# os.system("evince test_latency_"+name+".pdf")


all_mice = np.array([i.split(".")[0] for i in os.listdir("../data/")])
few_mice = np.array(['B28','B61','B137','B155','B163','B166','B150','B62','B84','B139','B154','B74','B86','B152'])
list_mice = np.array([i for i in all_mice if i.split("_")[0] in few_mice])
nb_mice = 50
# LATENCY
with open("../latency/latency.pickle", 'rb') as handle:
	latency = pickle.load(handle)

# # GRID SEARCH PARAMETERS
# bounds = dict({'gamma':[0.01, 0.99], 'beta':[1.0, 50.0], 'eta':[0.01, 0.99] })
# step = 15
# grid_step = dict({ 'gamma':np.linspace(bounds['gamma'][0], bounds['gamma'][1], step), 'beta':np.linspace(bounds['beta'][0], bounds['beta'][1], 9), 'eta':np.linspace(bounds['eta'][0], bounds['eta'][1], step) })
# for mouse in list_mice:	
# 	grid_data = np.load("../grid/grid6/grid_search_"+mouse+".npy")	
# 	ind_best = np.where(np.max(grid_data) == grid_data)
# 	parameters = dict([('WMsize',3)]+[(k,grid_step[k][ind_best[i][0]]) for k,i in zip(['eta','gamma','beta'],xrange(3))])
# 	test(mouse, parameters, latency[mouse.split("_")[0]])

# SFERES PARAMETERS
with open("../sferes/parameters.pickle", 'rb') as handle:
	parameters_sferes = pickle.load(handle)
for mouse in list_mice:		
	parameters = parameters_sferes['actorcritic'][mouse]
	parameters['WMsize'] = 3
	test(mouse, parameters, latency[mouse.split("_")[0]])







# f = open("parameters_grid.txt", 'w')

# for mouse in list_mice:	
# 	grid_data = np.load("../grid/grid6/grid_search_"+mouse+".npy")	
# 	ind_best = np.where(np.max(grid_data) == grid_data)
# 	parameters = dict([('WMsize',3)]+[(k,grid_step[k][ind_best[i][0]]) for k,i in zip(['eta','gamma','beta'],xrange(3))])
# 	line="mouse="+mouse.split("_")[0]+"\t"+" \t ".join([k+"="+str(np.round(parameters[k],3)) for k in ['beta','eta','gamma']])+"\tloglikelihood = "+str(-(10**np.abs(np.max(grid_data))))+"\n"
# 	f.writelines(line)
# f.close()





