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

# parameters = {'beta': 13.405800000000001,
#  'eta': 0.92488999999075105,
#  'gamma': 0.96420399990357952}
 
# model = Graph(parameters)

# wrap = TYMaze(model)
# wrap.world.startingPos()
# branch = {}

# for a in ["F", "F", "R", "U", "L"]:
# 	state = wrap.pos_to_state[wrap.world.mousePos]
# 	branch[]
# 	wrap.move(a)

rec(0, [])
