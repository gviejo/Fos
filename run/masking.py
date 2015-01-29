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
import cPickle as pickle



# sys.exit()

parameters = {'beta': 20.405800000000001,
 'eta': 0.9,
 'gamma': 0.96420399990357952}
 
model = PI(parameters)

wrap = TYMaze(model)
world = wrap.world

mask = {'Y':{p:{i:np.zeros((31,31,31*31)) for i in [1,2,3]} for p in ['I','II','V']}}
mask['I'] = {p:{i:np.zeros((31,31,31*31)) for i in [1,2]} for p in ['1b','10','2','3b','4','9b','8']}

next_states = {'V':['4','2','3b'],
				'I':['2','10','1b'],
				'II':['8','9b','10'],
				'1b':['I','1a'],
				'10':['II','I'],
				'2':['I','V'],
				'9b':['II','9a'],
				'8':['III','II'],
				'4':['IV','V'],
				'3b':['V','3a']}

xy = model.grid.reshape(31*31,2)

p1 = '1b'
for i in xrange(31*31):
	p = xy[i]
	# cadran 1 upper
	mask['I'][p1][1][:,:,i][model.grid[:,:,1]>p[1]] = 1.0
	# cadran 2 upper 
	mask['I'][p1][2][:,:,i][model.grid[:,:,1]<p[1]] = 2.0

for p1 in ['10', '2', '3b', '4', '9b', '8']:
	direction = model.computeAngle(p1, next_states[p1][0])
	arc = (direction+(np.pi/2.))%(2*np.pi)
	coeff = np.sin(arc)/np.cos(arc)
	for i in xrange(31*31):
		p = xy[i]
		mask['I'][p1][1][:,:,i][(model.grid[:,:,1]-(model.grid[:,:,0]*coeff+(p[1]-coeff*p[0])))>0] = 1.0
		mask['I'][p1][2][:,:,i][(model.grid[:,:,1]-(model.grid[:,:,0]*coeff+(p[1]-coeff*p[0])))<0] = 2.0

p1 = 'V'
direction = np.array([model.computeAngle(p1, p) for p in next_states[p1]])
#angle du milieu entre deux directions
arc = [direction[i]+(direction[i+1]-direction[i])/2. for i in xrange(2)]
arc.append(np.mean([direction[0],direction[-1]-2*np.pi]))
# rangé dans l'ordre croissant
arc = np.sort(arc)
# coefficient directeur de la droite
coeff = np.array([np.sin(arc[i])/np.cos(arc[i]) for i in xrange(3)])
for i in xrange(31*31):
	p = xy[i]	
	# cadran 1 right
	mask['Y'][p1][1][:,:,i][((model.grid[:,:,1]-(model.grid[:,:,0]*coeff[0]+(p[1]-coeff[0]*p[0])))<0) * ((model.grid[:,:,1]-(model.grid[:,:,0]*coeff[2]+(p[1]-coeff[2]*p[0])))>0)] = 1.0
	# cadran 2 upper
	mask['Y'][p1][2][:,:,i][((model.grid[:,:,1]-(model.grid[:,:,0]*coeff[1]+(p[1]-coeff[1]*p[0])))>0) * ((model.grid[:,:,1]-(model.grid[:,:,0]*coeff[0]+(p[1]-coeff[0]*p[0])))>0)] = 2.0
	# cadran 3 lower right
	mask['Y'][p1][3][:,:,i][((model.grid[:,:,1]-(model.grid[:,:,0]*coeff[1]+(p[1]-coeff[1]*p[0])))<0) * ((model.grid[:,:,1]-(model.grid[:,:,0]*coeff[2]+(p[1]-coeff[2]*p[0])))<0)] = 3.0


p1 = 'I'
direction = np.array([model.computeAngle(p1, p) for p in next_states[p1]])
#angle du milieu entre deux directions
arc = [direction[i]+(direction[i+1]-direction[i])/2. for i in xrange(2)]
arc.append(np.mean([direction[0],direction[-1]-2*np.pi])+2*np.pi)
# rangé dans l'ordre croissant
arc = np.sort(arc)
# coefficient directeur de la droite
coeff = np.array([np.sin(arc[i])/np.cos(arc[i]) for i in xrange(3)])
# careful ceoff[0] est une droite vertical
for i in xrange(31*31):
	p = xy[i]	
	# cadran 1 right
	mask['Y'][p1][1][:,:,i][(model.grid[:,:,0]>p[0]) * (((model.grid[:,:,1]-(model.grid[:,:,0]*coeff[2]+(p[1]-coeff[2]*p[0])))>0))] = 1.0
	# cadran 2 left
	mask['Y'][p1][2][:,:,i][(model.grid[:,:,0]<p[0]) * (((model.grid[:,:,1]-(model.grid[:,:,0]*coeff[1]+(p[1]-coeff[1]*p[0])))>0))] = 2.0
	# cadran 3 lower
	mask['Y'][p1][3][:,:,i][((model.grid[:,:,1]-(model.grid[:,:,0]*coeff[1]+(p[1]-coeff[1]*p[0])))<0) * ((model.grid[:,:,1]-(model.grid[:,:,0]*coeff[2]+(p[1]-coeff[2]*p[0])))<0)] = 3.0


p1 = 'II'
direction = np.array([model.computeAngle(p1, p) for p in next_states[p1]])
#angle du milieu entre deux directions
arc = [direction[i]+(direction[i+1]-direction[i])/2. for i in xrange(2)]
arc.append(np.mean([direction[0],direction[-1]-2*np.pi])+2*np.pi)
# rangé dans l'ordre croissant
arc = np.sort(arc)
# coefficient directeur de la droite
coeff = np.array([np.sin(arc[i])/np.cos(arc[i]) for i in xrange(3)])
# careful ceoff[0] est une droite vertical
for i in xrange(31*31):
	p = xy[i]	
	# cadran 1 lower right
	mask['Y'][p1][1][:,:,i][((model.grid[:,:,1]-(model.grid[:,:,0]*coeff[1]+(p[1]-coeff[1]*p[0])))<0) * ((model.grid[:,:,1]-(model.grid[:,:,0]*coeff[2]+(p[1]-coeff[2]*p[0])))<0)] = 1.0
	# cadran 2 upper
	mask['Y'][p1][2][:,:,i][((model.grid[:,:,1]-(model.grid[:,:,0]*coeff[2]+(p[1]-coeff[2]*p[0])))>0) * ((model.grid[:,:,1]-(model.grid[:,:,0]*coeff[0]+(p[1]-coeff[0]*p[0])))>0)] = 2.0
	# cadran 3 left
	mask['Y'][p1][3][:,:,i][((model.grid[:,:,1]-(model.grid[:,:,0]*coeff[0]+(p[1]-coeff[0]*p[0])))<0) * ((model.grid[:,:,1]-(model.grid[:,:,0]*coeff[1]+(p[1]-coeff[1]*p[0])))>0)] = 3.0



# with open("mask.pickle", 'wb') as handle:
# 	pickle.dump(mask, handle)


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
