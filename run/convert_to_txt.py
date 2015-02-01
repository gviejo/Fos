#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import os, sys
import cPickle as pickle
sys.path.append("../src/")
from world import World


# convert pos et state to int
states = ['I', 'Y', 'U']
pos = {'1a':0,'1b':1,'I':2,'10':3,'2':4,'V':8,'3b':9,'3a':10,'4':12,'IV':14,'II':7,'9a':5,'9b':6,'8':11,'III':13}


list_files = os.listdir("../data/")

list_files.remove("latency.pickle")

for f in list_files:
	with open("../data/"+f, 'rb') as handle:
		data = pickle.load(handle)
	mouse = f.split(".")[0][:-6]
	os.system("mkdir ccall/data_txt/"+mouse)
	tmp = []
	length_trials = ""
	for i in xrange(data['info']['nb_trial']):
		length_trials+=(str(data[i]['state'].shape[0])+" ")
		for j in xrange(len(data[i]['reward'])):
			line = [pos[data[i]['pos'][j]], states.index(data[i]['state'][j]), data[i]['action'][j], data[i]['reward'][j]]
			line += list((data[i]['possible'][j]).astype(int))
			tmp.append(line)		
		tmp.append([pos[data[i]['pos'][-1]], states.index(data[i]['state'][-1]), -1, -1, 0,0,0,0])
	tmp = np.array(tmp)
	with open("ccall/data_txt/"+mouse+"/possarpossible.txt", 'w') as f:
		for i in xrange(len(tmp)):
			line = " ".join(tmp[i].astype(str))
			f.write(line+"\n")
	with open("ccall/data_txt/"+mouse+"/info.txt", 'w') as f:
		f.write(length_trials)
	

