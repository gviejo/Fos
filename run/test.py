#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
import sys
sys.path.append("../src")
from Wrap import TYMaze
from Models import VMWM

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

parameters = dict({'length':2,
					'beta':2.0,
					'eta':0.1,
					'gamma':0.9})

model = VMWM(parameters)

wrap = TYMaze(model)
world = wrap.world
wrap.test(parameters, 1, 52)
# model.startExp()
# wrap.world.startingPos()

# # for i in [0,0,1,0,3,0]:
# # 	action(i)
# # 	sys.stdin.readline()
