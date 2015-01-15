#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys
from optparse import OptionParser                                                   
sys.path.append("../../src/")
from Wrap import TYMaze
import numpy as np
from Models import VMWM
import cPickle as pickle
import numpy as np
from matplotlib import *
from pylab import *

directory = '../sferes/SFERES5'
model = VMWM()
p_order = ['beta', 'gamma', 'eta', 'length']
model_in_folders = os.listdir(directory)
data = {}
m = 'VMWM'
lrun = os.listdir(directory+"/"+m)
order = p_order
scale = model.bounds
for r in lrun:
    s = "_".join(r.split("_")[2:4])[:-6]
    n = int(r.split("_")[4].split(".")[0])
    print("\t\t Mouse : "+s)
    n = int(r.split("_")[4].split(".")[0])
    if s in data.keys():
        data[s][n] = np.genfromtxt(directory+"/"+m+"/"+r)
    else :
        data[s] = dict()
        data[s][n] = np.genfromtxt(directory+"/"+m+"/"+r)                                
    for p in order:
        data[s][n][:,order.index(p)+4] = scale[p][0]+data[s][n][:,order.index(p)+4]*(scale[p][1]-scale[p][0])

parameters = dict()

for s in data.iterkeys():
	parameters[s] = dict()
	for n in data[s].iterkeys():		
		best_ind = np.argmax(data[s][n][:,2])
		parameters[s][n] = {p:data[s][n][best_ind,4+p_order.index(p)] for p in p_order}


with open("../sferes/SFERES5_parameters.pickle", 'wb') as f:
    pickle.dump(parameters, f)

# TEST
with open("../latency.pickle", 'rb') as handle:
  latency = pickle.load(handle)
nb_exp = 30
wrap = TYMaze(VMWM())
test = {}
cost = {}
for s in parameters.keys():
	test[s] = dict()
	cost[s] = []
	for n in parameters[s].iterkeys():
		print s, n
		test[s][n] = wrap.test(parameters[s][n], nb_exp, int(s.split("_")[1]))
		cost[s].append(np.sum(np.power(latency[s.split("_")[0]]-np.mean(test[s][n],0), 2)))
		# wrap.plot(data[s], s, parameters[m][s], latency[s.split("_")[0]], "../test/"+options.input.split("_")[0]+"_"+s+".pdf")


biais = np.linspace(0,10,20)


figure()
subplot(5,1,1)
for s in cost.iterkeys():
	plot(biais, np.array(cost[s]), 'o--')
	meancost = np.mean(np.array([cost[s] for s in cost.keys()]), 0)
	plot(biais, meancost, '-', linewidth = 4)

for p, i in zip(order,xrange(len(order))):
	subplot(5,1,i+2)
	for s in parameters.keys():
		plot(biais, np.array([parameters[s][j][p] for j in xrange(len(biais))]), 'o-')
		ylim(scale[p])
		ylabel(p)




show()

