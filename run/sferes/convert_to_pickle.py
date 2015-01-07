#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import os, sys
import cPickle as pickle
sys.path.append("../../src/")
from world import World

pos2pos = np.array(['undefined','1a','1b','I','10','II','9b','9a','8','III','2','V','3b','3a','4','IV'])
pos_to_state = dict({	'1a':'U',
						'1b':'I',
						'I' :'Y',
						'10':'I',
						'II':'Y',
						'8' :'U',
					   'III':'U',
					   	'9a':'U',
					   	'9b':'I',
					   	'2' :'I',
					   	'V' :'Y',
					   	'3b':'I',
					   	'3a':'U',
					   	'4' :'I',
					   	'IV':'U'})
actions = ['F', 'L', 'U', 'R']
files = os.listdir("/home/viejo/Benedicte/data")

for f in files:
	data = np.loadtxt("/home/viejo/Benedicte/data/"+f)
	nb_trials = len(np.unique(data[:,0]))
	nb_point = 0
	pick = {}
	for i in xrange(int(nb_trials)):
		pick[i] = {}
		trial = data[data[:,0] == i+1]
		pos_sequence = pos2pos[trial[:,1].astype(int)]
		state_sequence = np.array([pos_to_state[p] for p in pos_sequence])
		action_sequence = trial[0:-1,2].astype(int)
		reward_sequence = np.zeros(len(action_sequence)).astype(int)
		if pos_sequence[-1] == 'III':
			reward_sequence[-1] = 1
		pick[i]['state'] = state_sequence
		pick[i]['action'] = action_sequence
		pick[i]['reward'] = reward_sequence
		pick[i]['pos'] = pos_sequence
		pick[i]['ind'] = np.arange(len(action_sequence))+nb_point # position dans le vecteur log-likelihood
		# on instancie un world pour aller plus vite et examiner les possible a chaque etape
		world = World("TYM")
		world.startingPos()
		possible = [world.readPathwaysALaLouche()]
		for j in xrange(len(action_sequence)):
			world.moveALaLouche(action_sequence[j])
			possible.append(world.readPathwaysALaLouche())
		possible = np.array(possible)		
		nb_point+=(len(action_sequence))
		pick[i]['possible'] = possible
	pick['info'] = {'nb_point':nb_point,
					'nb_trial':nb_trials}
	with open("../../data/"+f.split(".")[0]+".pickle", 'wb') as handle:
		pickle.dump(pick, handle)
		
