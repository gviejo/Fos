#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
import sys
sys.path.append("../src")
from Wrap import TYMaze
from Models import Graph
from matplotlib import *
from pylab import *

def action(a):		
	print "Position", world.mousePos
	print "Possible", world.readPathwaysALaLouche()
	model.computeValue(wrap.pos_to_state[world.mousePos],a,world.readPathwaysALaLouche())	
	wrap.move(wrap.actions[a])
	state = wrap.pos_to_state[world.mousePos]
	reward = world.readRew()	
	model.updateValue(reward, state)
	
	print "Delta ", model.delta
	# print "Possible", np.array(wrap.actions)[np.array(world.readPathwaysALaLouche())==1]

	print "Reward", reward

parameters = {'beta': 13.405800000000001,
 'eta': 0.92488999999075105,
 'gamma': 0.96420399990357952}
 



model = Graph(parameters)

wrap = TYMaze(model)
# world = wrap.world
# data = wrap.test(parameters, 1, 52)

model.startExp()
wrap.world.startingPos()

# # for i in [0,0,1,0,3,0]:
# # 	action(i)
# # 	sys.stdin.readline()
