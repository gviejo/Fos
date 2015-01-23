#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
import sys
sys.path.append("../src")
from Wrap import TYMaze
from Models import Graph
from world import World
from matplotlib import *
from pylab import *

count = []


def rec(n, previous):		
	world = World("TYM")
	world.startingPos()
	for a in previous:
		world.moveALaLouche(a)
	if n == 57:
		count.append(1)
		return
	elif world.readRew():
		count.append(1)
		return
	else:		
		possible = world.readPathwaysALaLouche()		
		for i in np.where(possible)[0]:
			rec(n+1, previous+[i])

parameters = {'beta': 20.405800000000001,
 'eta': 0.9,
 'gamma': 0.96420399990357952}
 
model = Graph(parameters)

wrap = TYMaze(model)

wrap.guidage()

for i in xrange(1):
	model.startTrial()
	wrap.world.startingPos()
	state = wrap.pos_to_state[wrap.world.mousePos]	
	for j in xrange(10):	
	# for s,a in zip(wrap.state_guiding, wrap.action_guiding):
		possible = wrap.world.readPathwaysALaLouche()
		print "state=",state
		action = model.chooseAction(state, wrap.world.readPathwaysALaLouche())
		# p_a = model.computeValue(s, a, possible)
		print "node=",model.current_node
		print "qvalues=",model.q_values[possible==1]
		print "action=",action
		wrap.move(action)
		state = wrap.pos_to_state[wrap.world.mousePos]
		reward = wrap.world.readRew()
		# print model.nodes
		# print model.edges		
		# print model.T
		# print model.V
		wrap.model.updateValue(reward, state)	
		if reward:
			print "reward found"
			break
		sys.stdin.readline()
# rec(0, [])
