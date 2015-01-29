#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
import sys
sys.path.append("../src")
from Wrap import TYMaze
from Models import *
from world import World
from matplotlib import *
from pylab import *



parameters = {'beta': 20.405800000000001,
 'eta': 0.9,
 'gamma': 0.96420399990357952}
 
model = PI(parameters)

wrap = TYMaze(model)
world = wrap.world
# wrap.guidage()
# sys.exit()
# for i in xrange(1):
# 	model.startTrial()
# 	wrap.world.startingPos()
# 	for j in xrange(10):	
# 		possible = wrap.world.readPathwaysALaLouche()
# 		print "position=",world.mousePos," Possible=",possible
# 		action = model.chooseAction(world.mousePos, wrap.pos_to_state[wrap.world.mousePos], wrap.world.readPathwaysALaLouche())
# 		# # p_a = model.computeValue(s, a, possible)		

# 		# print "qvalues=",model.q_values[possible==1]
# 		# print "action=",action
# 		wrap.move(action)		
# 		reward = wrap.world.readRew()
# 		model.updateValue(reward, wrap.pos_to_state[wrap.world.mousePos])
# 		print "V(pos)=", model.varPos
# 		print "V(goal)=",model.varGoal
# 		sys.stdin.readline()



# transition = set()

# world.startingPos()
# old_position = world.mousePos
# world.moveALaLouche(0)
# for a in [0,1,0,1,0,2,0,1,0,2,0,1,0,1,0,1,0,2,0,1,0,2,0,1,0,1,0,2]:
# 	position = world.mousePos
# 	if np.sum(world.readPathwaysALaLouche())>1:
# 		transition.add((old_position,position,a))		
# 	world.moveALaLouche(a)
# 	old_position = position	

# world.startingPos()
# old_position = world.mousePos
# world.moveALaLouche(0)
# for a in [0,3,0,3,0,2,0,3,0,2,0,3,0,3,0,3,0,2,0,3,0,2,0,3,0,3,0,2]:
# 	position = world.mousePos
# 	if np.sum(world.readPathwaysALaLouche())>1:
# 		transition.add((old_position,position,a))
# 	world.moveALaLouche(a)
# 	old_position = position	

# tmp = dict.fromkeys(transition)
# alpha = np.linspace(0,2*np.pi, 61)

# world.startingPos()
# old_position = world.mousePos
# world.moveALaLouche(0)
# for a in [0,1,0,1,0,2,0,1,0,2,0,1,0,1,0,1,0,2,0,1,0,2,0,1,0,1,0,2]:
# 	position = world.mousePos
# 	# if (old_position,position,a) in tmp.keys():
# 	# 	tmp[(old_position,position,a)] = [world.mouseDir]
# 	world.moveALaLouche(a)	
# 	if (old_position,position,a) in tmp.keys():
# 		tmp[(old_position,position,a)] = world.mousePos
# 	old_position = position	


# world.startingPos()
# old_position = world.mousePos
# world.moveALaLouche(0)
# for a in [0,3,0,3,0,2,0,3,0,2,0,3,0,3,0,3,0,2,0,3,0,2,0,3,0,3,0,2]:
# 	position = world.mousePos
# 	# if (old_position,position,a) in tmp.keys():
# 	# 	tmp[(old_position,position,a)] = [world.mouseDir]
# 	world.moveALaLouche(a)
# 	if (old_position,position,a) in tmp.keys():
# 		tmp[(old_position,position,a)] = world.mousePos
# 	old_position = position	

