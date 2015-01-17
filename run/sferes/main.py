#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from optparse import OptionParser                                                   
sys.path.append("../../src/")
from Wrap import TYMaze
import numpy as np
from Models import VMWM
import cPickle as pickle
from matplotlib import *
from pylab import *

with open("SFERES9_parameters.pickle", 'rb') as handle:
	parameters = pickle.load(handle)['VMWM']



for s in parameters.iterkeys():
	with open("../../data/"+s+"trials.pickle") as f:
		data = pickle.load(f)
	epsilon = 0.0		
	tmp = []
	gamma = [0.01, 0.5,  parameters[s]['gamma']]
	beta = [0.1, 1.0, parameters[s]['beta'], 8.0]
	parameters[s]['length'] = 3
	p_test = parameters[s]
	# for b in beta :
	for g in gamma:		
		# p_test['beta'] = b
		p_test['gamma'] = g
		# p_test['length'] = 3

		model = VMWM(parameters[s])
		opt = TYMaze(model)
		llh = opt.sferes(data, epsilon)
		# print llh

		tmp.append(llh)
		
	tmp = np.array(tmp)
	# Plot steps vs mean loglikelihood
	x = np.unique(tmp[0,:,0])
	y = np.array([[np.mean(tmp[i,:,1][tmp[i,:,0]==j]) for j in x] for i in xrange(len(tmp))])
	v = np.array([[np.std(tmp[i,:,1][tmp[i,:,0]==j]) for j in x] for i in xrange(len(tmp))])

	print s, tmp[:,:,1].sum(1)
	figure()
	
	subplot(211)
	[plot(tmp[i][:,1], 'o-', label=gamma[i]) for i in xrange(len(tmp))]
	legend(loc = 'best')
	subplot(212)
	[errorbar(x, y[i], v[i], label=gamma[i]) for i in xrange(len(tmp))]
	legend(loc = 'best')

	show()