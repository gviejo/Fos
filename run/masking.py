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




# sys.exit()

parameters = {'beta': 20.405800000000001,
 'eta': 0.9,
 'gamma': 0.96420399990357952}
 
model = PI(parameters)

wrap = TYMaze(model)
world = wrap.world

mask = {'Y':{p:np.zeros((31,31,31*31)) for p in ['I','II','V']}}

transitions = {'V':['4','2','3b'],
				'I':['2','10','1b'],
				'II':['8','9b','10']}

positions = model.grid.reshape(31*31,2)



p1 = 'V'
direction = np.array([model.computeAngle(p1, p) for p in transitions[p1]])
#angle du milieu entre deux directions
arc = [direction[i]+(direction[i+1]-direction[i])/2. for i in xrange(2)]
arc.append(np.mean([direction[0],direction[-1]-2*np.pi]))
# rangé dans l'ordre croissant
arc = np.sort(arc)
# coefficient directeur de la droite
coeff = np.array([np.sin(arc[i])/np.cos(arc[i]) for i in xrange(3)])
for i in xrange(31*31):
	p = positions[i]	
	# cadran 1 right
	mask['Y'][p1][:,:,i][((model.grid[:,:,1]-(model.grid[:,:,0]*coeff[0]+(p[1]-coeff[0]*p[0])))<0) * ((model.grid[:,:,1]-(model.grid[:,:,0]*coeff[2]+(p[1]-coeff[2]*p[0])))>0)] = 1.0
	# cadran 2 upper
	mask['Y'][p1][:,:,i][((model.grid[:,:,1]-(model.grid[:,:,0]*coeff[1]+(p[1]-coeff[1]*p[0])))>0) * ((model.grid[:,:,1]-(model.grid[:,:,0]*coeff[0]+(p[1]-coeff[0]*p[0])))>0)] = 2.0
	# cadran 3 lower right
	mask['Y'][p1][:,:,i][((model.grid[:,:,1]-(model.grid[:,:,0]*coeff[1]+(p[1]-coeff[1]*p[0])))<0) * ((model.grid[:,:,1]-(model.grid[:,:,0]*coeff[2]+(p[1]-coeff[2]*p[0])))<0)] = 3.0


p1 = 'I'
direction = np.array([model.computeAngle(p1, p) for p in transitions[p1]])
#angle du milieu entre deux directions
arc = [direction[i]+(direction[i+1]-direction[i])/2. for i in xrange(2)]
arc.append(np.mean([direction[0],direction[-1]-2*np.pi])+2*np.pi)
# rangé dans l'ordre croissant
arc = np.sort(arc)
# coefficient directeur de la droite
coeff = np.array([np.sin(arc[i])/np.cos(arc[i]) for i in xrange(3)])
# careful ceoff[0] est une droite vertical
for i in xrange(31*31):
	p = positions[i]	
	# cadran 1 right
	mask['Y'][p1][:,:,i][(model.grid[:,:,0]>p[0]) * (((model.grid[:,:,1]-(model.grid[:,:,0]*coeff[2]+(p[1]-coeff[2]*p[0])))>0))] = 1.0
	# cadran 2 left
	mask['Y'][p1][:,:,i][(model.grid[:,:,0]<p[0]) * (((model.grid[:,:,1]-(model.grid[:,:,0]*coeff[1]+(p[1]-coeff[1]*p[0])))>0))] = 2.0
	# cadran 3 lower
	mask['Y'][p1][:,:,i][((model.grid[:,:,1]-(model.grid[:,:,0]*coeff[1]+(p[1]-coeff[1]*p[0])))<0) * ((model.grid[:,:,1]-(model.grid[:,:,0]*coeff[2]+(p[1]-coeff[2]*p[0])))<0)] = 3.0


p1 = 'II'
direction = np.array([model.computeAngle(p1, p) for p in transitions[p1]])
#angle du milieu entre deux directions
arc = [direction[i]+(direction[i+1]-direction[i])/2. for i in xrange(2)]
arc.append(np.mean([direction[0],direction[-1]-2*np.pi])+2*np.pi)
# rangé dans l'ordre croissant
arc = np.sort(arc)
# coefficient directeur de la droite
coeff = np.array([np.sin(arc[i])/np.cos(arc[i]) for i in xrange(3)])
# careful ceoff[0] est une droite vertical
for i in xrange(31*31):
	p = positions[i]	
	# cadran 1 lower right
	mask['Y'][p1][:,:,i][((model.grid[:,:,1]-(model.grid[:,:,0]*coeff[1]+(p[1]-coeff[1]*p[0])))<0) * ((model.grid[:,:,1]-(model.grid[:,:,0]*coeff[2]+(p[1]-coeff[2]*p[0])))<0)] = 1.0
	# cadran 2 upper
	mask['Y'][p1][:,:,i][((model.grid[:,:,1]-(model.grid[:,:,0]*coeff[2]+(p[1]-coeff[2]*p[0])))>0) * ((model.grid[:,:,1]-(model.grid[:,:,0]*coeff[0]+(p[1]-coeff[0]*p[0])))>0)] = 2.0
	# cadran 3 left
	mask['Y'][p1][:,:,i][((model.grid[:,:,1]-(model.grid[:,:,0]*coeff[0]+(p[1]-coeff[0]*p[0])))<0) * ((model.grid[:,:,1]-(model.grid[:,:,0]*coeff[1]+(p[1]-coeff[1]*p[0])))>0)] = 3.0



# figure()
# for i,j in zip(np.arange(25),np.arange(0,31*31,25)):
# 	subplot(5,5,i+1)
# 	imshow(mask['Y'][p1][:,:,j], origin='lower')
# show()


# for i in xrange(31*31):
# 	imshow(mask['Y'][p1][:,:,i], origin = 'lower')
# 	print i
# 	show()
# figure()
# subplot(211)
# # [plot([0,np.cos(direction)[i]], [0,np.sin(direction)[i]], '-') for i in xrange(3)]
# [plot([0,np.cos(arc)[i]], [0,np.sin(arc)[i]], '-') for i in xrange(3)]

# subplot(212)
# print p


# imshow(cadran, origin = 'lower')





# x = np.arange(-0.2,0.2,0.01)
# [plot(x, x*coeff[i], 'o') for i in xrange(3)]

# show()
