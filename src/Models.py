#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Models.py

model-free reinforcement + working memory = VMWM
"""


import numpy as np
import sys, os
import scipy.stats
states = np.array(['I', 'Y', 'U'])
actions = np.array(['F', 'L', 'U', 'R'])


class VMWM():
	""" from B. Babayan thesis
	"""

	def __init__(self, parameters = {'length':3}):
		# Parameters
		self.parameters = parameters
		self.parameters['length'] = int(self.parameters['length'])
		self.n_action = len(actions)
		self.n_state = len(states)
		self.n_possible = self.n_state*np.sum(map(lambda x:self.n_action**x, range(self.parameters['length']+1)))
		self.bounds = dict({"length":[0,3.9],
							"gamma":[0.0, 0.9999999999],
							"beta":[0.0, 200.0],
							"eta":[0.0, 0.99999999999]})
		self.delta = 0.0		
		# Values initialization		
		# self.critic	= np.ones(self.n_possible)*0.1
		# self.actor = np.ones((self.n_possible, self.n_action))*0.1
		self.critic	= np.zeros(self.n_possible)
		self.actor = np.zeros((self.n_possible, self.n_action))
		#Various Init
		self.current_state = None
		self.current_action = None
		self.last_actions = np.array(['']*self.parameters['length'])
		self.ind = np.arange(self.n_action)
		# Conversion Init				
		self.convert = dict()
		self.keys = []		
		for s in states: self.initConversion(s, self.parameters['length'])			
		self.keys = np.array(self.keys)		
		assert len(self.keys) == self.n_possible
		for i in xrange(self.n_possible):
			self.convert[self.keys[i]] = i

	def initConversion(self, key, n):
		self.keys.append(key)
		if n == 0:
			return
		else:
			for a in actions:
				self.initConversion(key+a, n-1)

	def setParameters(self, name, value):            
		if value < self.bounds[name][0]:
			self.parameters[name] = self.bounds[name][0]
		elif value > self.bounds[name][1]:
			self.parameters[name] = self.bounds[name][1]
		else:
			self.parameters[name] = value                

	def setAllParameters(self, parameters):
		for i in parameters.iterkeys():
			if i in self.bounds.keys():
				self.setParameters(i, parameters[i])

	def startExp(self):
		# self.critic	= np.ones(self.n_possible)*0.1
		# self.actor = np.ones((self.n_possible, self.n_action))*0.1
		self.critic	= np.zeros(self.n_possible)
		self.actor = np.zeros((self.n_possible, self.n_action))

	def startTrial(self):
		self.last_actions = np.array(['']*self.parameters['length'])

	def softMax(self, values):
		tmp = np.exp(values*float(self.parameters['beta']))
		if np.isinf(tmp).sum() == 1:
			tmp = np.isinf(tmp)*(1.0-len(tmp)*1e-8) + 1e-8
		p_a = tmp/float(np.sum(tmp))
		if np.sum(p_a==0.0):
			p_a += 1e-8
			p_a = p_a/float(np.sum(p_a))			
		return p_a

	def sampleSoftMax(self, values):
		tmp = np.exp(values*float(self.parameters['beta']))
		if np.isinf(tmp).sum() == 1:
			tmp = np.isinf(tmp)*(1.0-len(tmp)*1e-8) + 1e-8
		p_a = tmp/float(np.sum(tmp))
		if np.sum(p_a==0.0):
			p_a += 1e-8
			p_a = p_a/float(np.sum(p_a))	
		tmp = [np.sum(p_a[0:i]) for i in range(len(p_a))]
		return np.sum(np.array(tmp) < np.random.rand())-1 

	def computeValue(self, pos, state, a, possible):
		# Very tricky : state in [U,Y,I] and a in [0,1,2,3]
		self.current_state = self.convert[state+"".join(self.last_actions)]
		self.current_action = a
		ind = self.ind[possible==1]
		q_values = self.actor[self.current_state][possible==1]
		p_a = np.zeros(self.n_action)		
		p_a[ind] = self.softMax(q_values)			
		return p_a[self.current_action]

	def chooseAction(self, pos, state, possible):
		self.current_state = self.convert[state+"".join(self.last_actions)]
		ind = self.ind[possible==1]
		q_values = self.actor[self.current_state][possible==1]		
		self.current_action = ind[self.sampleSoftMax(q_values)]
		return actions[self.current_action]		

	def updateValue(self, reward, next_state):
		r = (reward==0)*0.0+(reward==1)*1.0+(reward==-1)*0.0		        
		# Update list of previous action
		if self.parameters['length']:
			self.last_actions[1:] = self.last_actions[:-1]
			self.last_actions[0] = actions[self.current_action]
		# Delta update
		n_s = self.convert[next_state+"".join(self.last_actions)]
		self.delta = r + self.parameters['gamma']*self.critic[n_s] - self.critic[self.current_state]
		self.critic[self.current_state] += (self.parameters['eta']*self.delta)
		self.actor[self.current_state, self.current_action] += (self.parameters['eta']*self.delta)

class Graph():
	""" from B.Babayan thesis
	"""
	def __init__(self, parameters={}):
		# Parameters
		self.parameters = parameters	
		self.n_action = len(actions)
		self.n_state = len(states)
		self.bounds = dict({"gamma":[0.0, 0.9999999999],
							"beta":[0.0, 200.0],
							"eta":[0.0, 0.99999999999]})		
		# Graph initialization		
		self.nodes = {0:np.zeros(self.n_action,dtype=int)} # node : [node1,0,node4,0]
		self.states = {0:'U'}
		self.edges = {0:[]} # node:[a1, a2, a3]
		self.current_node = 0
		self.T = {0:np.zeros(self.n_action,dtype=float)}
		self.V = {0:0.0}
		self.R = {0:0.0}
		self.back = {}
		#Various Init
		self.ind = np.arange(self.n_action)
		self.q_values = np.zeros(self.n_action)

	def setParameters(self, name, value):            
		if value < self.bounds[name][0]:
			self.parameters[name] = self.bounds[name][0]
		elif value > self.bounds[name][1]:
			self.parameters[name] = self.bounds[name][1]
		else:
			self.parameters[name] = value                

	def setAllParameters(self, parameters):
		for i in parameters.iterkeys():
			if i in self.bounds.keys():
				self.setParameters(i, parameters[i])

	def startExp(self):
		# Graph initialization		
		self.nodes = {0:np.zeros(self.n_action,dtype=int)} # node : [node1,0,node4,0]
		self.states = {0:'U'}
		self.edges = {0:[]} # node:[a1, a2, a3]
		self.current_node = 0
		self.T = {0:np.zeros(self.n_action,dtype=float)}
		self.V = {0:0.0}
		self.R = {0:0.0}
		self.back = {}
		#Various Init
		self.ind = np.arange(self.n_action)
		self.q_values = np.zeros(self.n_action)		

	def startTrial(self):
		self.current_node = 0

	def softMax(self, values):
		tmp = np.exp(values*float(self.parameters['beta']))
		if np.isinf(tmp).sum() == 1:
			tmp = np.isinf(tmp)*(1.0-len(tmp)*1e-8) + 1e-8
		p_a = tmp/float(np.sum(tmp))
		if np.sum(p_a==0.0):
			p_a += 1e-8
			p_a = p_a/float(np.sum(p_a))			
		return p_a

	def sampleSoftMax(self, values):
		tmp = np.exp(values*float(self.parameters['beta']))
		if np.isinf(tmp).sum() == 1:
			tmp = np.isinf(tmp)*(1.0-len(tmp)*1e-8) + 1e-8
		p_a = tmp/float(np.sum(tmp))
		if np.sum(p_a==0.0):
			p_a += 1e-8
			p_a = p_a/float(np.sum(p_a))	
		tmp = [np.sum(p_a[0:i]) for i in range(len(p_a))]
		return np.sum(np.array(tmp) < np.random.rand())-1 

	def retropropagate(self, node):		
		self.V[node] = np.max([self.R[node],self.V[node],np.max([self.parameters['gamma']*self.T[node][a]*self.V[self.nodes[node][a]] for a in self.edges[node]])])
		if self.back[node]:
			self.retropropagate(self.back[node])
		else:
			return

	def computeValue(self, pos, state, a, possible):
		# Very tricky : state in [U,Y,I] and a in [0,1,2,3]
		self.current_state = self.states[self.current_node]
		if self.current_state != state:
			print "problem" 
			sys.exit()		
		self.current_action = a
		ind = self.ind[possible==1]
		self.q_values *= 0.0
		for a in self.edges[self.current_node]:
			self.q_values[a] = self.parameters['gamma']*self.T[self.current_node][a]*self.V[self.nodes[self.current_node][a]]		
		p_a = np.zeros(self.n_action)		
		p_a[ind] = self.softMax(self.q_values[possible==1])			
		return p_a[self.current_action]

	def chooseAction(self, pos, state, possible):
		self.current_state = self.states[self.current_node]		
		if self.current_state != state:
			print "problem"
			sys.exit()
		ind = self.ind[possible==1]
		self.q_values *= 0.0
		for a in self.edges[self.current_node]:
			self.q_values[a] = self.parameters['gamma']*self.T[self.current_node][a]*self.V[self.nodes[self.current_node][a]]	
		self.current_action = ind[self.sampleSoftMax(self.q_values[possible==1])]
		return actions[self.current_action]		

	def updateValue(self, reward, next_state):
		r = (reward==0)*0.0+(reward==1)*1.0+(reward==-1)*0.0		        		
		if self.current_action in self.edges[self.current_node]:
			new_node = self.nodes[self.current_node][self.current_action]
			self.T[self.current_node][self.current_action] += (self.parameters['eta']*(1.0-self.T[self.current_node][self.current_action]))			
		else:
			new_node=np.max(self.nodes.keys())+1
			self.nodes[self.current_node][self.current_action] = new_node			
			self.edges[self.current_node].append(self.current_action)
			self.T[self.current_node][self.current_action] = self.parameters['eta']
			self.nodes[new_node] = np.zeros(self.n_action, dtype=int)
			self.states[new_node] = next_state
			self.edges[new_node] = []						
			self.T[new_node] = np.zeros(self.n_action, dtype=float)
			self.R[new_node] = r
			self.V[new_node] = r			
			self.back[new_node] = self.current_node
		self.current_node = new_node
		self.current_state = next_state
		if r:
			self.retropropagate(self.back[self.current_node])



			

		# Delta update
		# n_s = self.convert[next_state+"".join(self.last_actions)]
		# self.delta = r + self.parameters['gamma']*self.critic[n_s] - self.critic[self.current_state]
		# self.critic[self.current_state] += (self.parameters['eta']*self.delta)
		# self.actor[self.current_state, self.current_action] += (self.parameters['eta']*self.delta)

class PI():
	""" from B. Babayan thesis
	"""

	def __init__(self, parameters = {}):
		# Parameters
		self.parameters = parameters		
		self.n_action = len(actions)
		self.n_state = len(states)		
		self.bounds = dict({"gamma":[0.0, 0.9999999999],
							"beta":[0.0, 200.0],
							"eta":[0.0, 0.99999999999]})		
		# Values initialization
		self.transition = {('1a','1b',0):1,
							('1a','1b',2):2,
							('9a','9b',0):1,
							('9a','9b',2):2,
							('3a','3b',0):1,
							('3a','3b',2):2,
							('III','8',0):2,
							('III','8',2):1,
							('IV','4',0):2,
							('IV','4',2):1,							
							('I','1b',0):2,
							('I','1b',2):1,
							('I','10',0):1,
							('I','10',2):2,
							('I','2',0):1,
							('I','2',2):2,
							('II','9b',0):2,
							('II','9b',2):1,
							('II','8',0):1,
							('II','8',2):2,
							('II','10',0):2,
							('II','10',2):1,
							('V','2',0):2,
							('V','2',2):1,
							('V','3b',0):2,
							('V','3b',2):1,
							('V','4',0):1,
							('V','4',2):2,
							('1b','I',3):1,
							('1b','I',1):2,
							('1b','I',2):3,
							('10','I',1):1,
							('10','I',2):2,
							('10','I',3):3,
							('2','I',2):1,
							('2','I',3):2,
							('2','I',1):3,							
							('10','II',2):1,
							('10','II',3):2,
							('10','II',1):3,
							('9b','II',3):1,
							('9b','II',1):2,
							('9b','II',2):3,
							('8','II',1):1,
							('8','II',2):2,
							('8','II',3):3,							
							('3b','V',2):1,
							('3b','V',3):2,
							('3b','V',1):3,
							('2','V',3):1,
							('2','V',1):2,
							('2','V',2):3,
							('4','V',1):1,
							('4','V',2):2,
							('4','V',3):3}		
		self.positions = {'10': [-0.19, 0.608],
						 '1a': [0, 0],
						 '1b': [0.0, 0.235],
						 '2': [0.19, 0.608],
						 '3a': [0.827, 0.601],
						 '3b': [0.604, 0.674],
						 '4': [0.308, 0.97],
						 '8': [-0.308, 0.97],
						 '9a': [-0.827, 0.601],
						 '9b': [-0.604, 0.674],
						 'I': [0.0, 0.47],
						 'II': [-0.38, 0.746],
						 'III': [-0.235, 1.193],
						 'IV': [0.235, 1.193],
						 'V': [0.38, 0.746]}
		self.varPos = 0.0
		self.varGoal = 0.0				
		# Matrix init
		self.n_case = 30
		self.grain = 6./self.n_case	
		self.grid = np.dstack(np.meshgrid(np.linspace(-3,3,self.n_case+1),np.linspace(-3,3,self.n_case+1)))
		self.Pgoal = np.zeros((self.n_case+1, self.n_case+1))
		self.Ppos = np.zeros((self.n_case+1, self.n_case+1))
		self.xy = self.grid.reshape(31*31,2)
		# Mask init
		self.mask = {'Y':{p:{i:np.zeros((31,31,31*31)) for i in [1,2,3]} for p in ['I','II','V']}}
		self.mask['I'] = {p:{i:np.zeros((31,31,31*31)) for i in [1,2]} for p in ['1b','10','2','3b','4','9b','8']}
		self.next_states = {'V':['4','2','3b'],
							'I':['2','10','1b'],
							'II':['8','9b','10'],
							'1b':['I','1a'],
							'10':['II','I'],
							'2':['I','V'],
							'9b':['II','9a'],
							'8':['III','II'],
							'4':['IV','V'],
							'3b':['V','3a']}
		self.fillMaskArray()
		#Various Init
		self.q_values = None
		self.current_position = '1a'
		self.previous_position = '1a'
		self.current_action = None		
		self.ind = np.arange(self.n_action)
		self.reward_position = self.positions['III']		
		self.reward_found = False

	def fillMaskArray(self):
		p1 = '1b'
		for i in xrange(31*31):
			p = self.xy[i]
			# cadran 1 upper
			self.mask['I'][p1][1][:,:,i][self.grid[:,:,1]>p[1]] = 1.0
			# cadran 2 upper 
			self.mask['I'][p1][2][:,:,i][self.grid[:,:,1]<p[1]] = 2.0

		for p1 in ['10', '2', '3b', '4', '9b', '8']:
			direction = self.computeAngle(p1, self.next_states[p1][0])
			arc = (direction+(np.pi/2.))%(2*np.pi)
			coeff = np.sin(arc)/np.cos(arc)
			for i in xrange(31*31):
				p = self.xy[i]
				self.mask['I'][p1][1][:,:,i][(self.grid[:,:,1]-(self.grid[:,:,0]*coeff+(p[1]-coeff*p[0])))>0] = 1.0
				self.mask['I'][p1][2][:,:,i][(self.grid[:,:,1]-(self.grid[:,:,0]*coeff+(p[1]-coeff*p[0])))<0] = 2.0

		p1 = 'V'
		direction = np.array([self.computeAngle(p1, p) for p in self.next_states[p1]])
		#angle du milieu entre deux directions
		arc = [direction[i]+(direction[i+1]-direction[i])/2. for i in xrange(2)]
		arc.append(np.mean([direction[0],direction[-1]-2*np.pi]))
		# rangé dans l'ordre croissant
		arc = np.sort(arc)
		# coefficient directeur de la droite
		coeff = np.array([np.sin(arc[i])/np.cos(arc[i]) for i in xrange(3)])
		for i in xrange(31*31):
			p = self.xy[i]	
			# cadran 1 right
			self.mask['Y'][p1][1][:,:,i][((self.grid[:,:,1]-(self.grid[:,:,0]*coeff[0]+(p[1]-coeff[0]*p[0])))<0) * ((self.grid[:,:,1]-(self.grid[:,:,0]*coeff[2]+(p[1]-coeff[2]*p[0])))>0)] = 1.0
			# cadran 2 upper
			self.mask['Y'][p1][2][:,:,i][((self.grid[:,:,1]-(self.grid[:,:,0]*coeff[1]+(p[1]-coeff[1]*p[0])))>0) * ((self.grid[:,:,1]-(self.grid[:,:,0]*coeff[0]+(p[1]-coeff[0]*p[0])))>0)] = 2.0
			# cadran 3 lower right
			self.mask['Y'][p1][3][:,:,i][((self.grid[:,:,1]-(self.grid[:,:,0]*coeff[1]+(p[1]-coeff[1]*p[0])))<0) * ((self.grid[:,:,1]-(self.grid[:,:,0]*coeff[2]+(p[1]-coeff[2]*p[0])))<0)] = 3.0


		p1 = 'I'
		direction = np.array([self.computeAngle(p1, p) for p in self.next_states[p1]])
		#angle du milieu entre deux directions
		arc = [direction[i]+(direction[i+1]-direction[i])/2. for i in xrange(2)]
		arc.append(np.mean([direction[0],direction[-1]-2*np.pi])+2*np.pi)
		# rangé dans l'ordre croissant
		arc = np.sort(arc)
		# coefficient directeur de la droite
		coeff = np.array([np.sin(arc[i])/np.cos(arc[i]) for i in xrange(3)])
		# careful ceoff[0] est une droite vertical
		for i in xrange(31*31):
			p = self.xy[i]	
			# cadran 1 right
			self.mask['Y'][p1][1][:,:,i][(self.grid[:,:,0]>p[0]) * (((self.grid[:,:,1]-(self.grid[:,:,0]*coeff[2]+(p[1]-coeff[2]*p[0])))>0))] = 1.0
			# cadran 2 left
			self.mask['Y'][p1][2][:,:,i][(self.grid[:,:,0]<p[0]) * (((self.grid[:,:,1]-(self.grid[:,:,0]*coeff[1]+(p[1]-coeff[1]*p[0])))>0))] = 2.0
			# cadran 3 lower
			self.mask['Y'][p1][3][:,:,i][((self.grid[:,:,1]-(self.grid[:,:,0]*coeff[1]+(p[1]-coeff[1]*p[0])))<0) * ((self.grid[:,:,1]-(self.grid[:,:,0]*coeff[2]+(p[1]-coeff[2]*p[0])))<0)] = 3.0


		p1 = 'II'
		direction = np.array([self.computeAngle(p1, p) for p in self.next_states[p1]])
		#angle du milieu entre deux directions
		arc = [direction[i]+(direction[i+1]-direction[i])/2. for i in xrange(2)]
		arc.append(np.mean([direction[0],direction[-1]-2*np.pi])+2*np.pi)
		# rangé dans l'ordre croissant
		arc = np.sort(arc)
		# coefficient directeur de la droite
		coeff = np.array([np.sin(arc[i])/np.cos(arc[i]) for i in xrange(3)])		
		for i in xrange(31*31):
			p = self.xy[i]	
			# cadran 1 lower right
			self.mask['Y'][p1][1][:,:,i][((self.grid[:,:,1]-(self.grid[:,:,0]*coeff[1]+(p[1]-coeff[1]*p[0])))<0) * ((self.grid[:,:,1]-(self.grid[:,:,0]*coeff[2]+(p[1]-coeff[2]*p[0])))<0)] = 1.0
			# cadran 2 upper
			self.mask['Y'][p1][2][:,:,i][((self.grid[:,:,1]-(self.grid[:,:,0]*coeff[2]+(p[1]-coeff[2]*p[0])))>0) * ((self.grid[:,:,1]-(self.grid[:,:,0]*coeff[0]+(p[1]-coeff[0]*p[0])))>0)] = 2.0
			# cadran 3 left
			self.mask['Y'][p1][3][:,:,i][((self.grid[:,:,1]-(self.grid[:,:,0]*coeff[0]+(p[1]-coeff[0]*p[0])))<0) * ((self.grid[:,:,1]-(self.grid[:,:,0]*coeff[1]+(p[1]-coeff[1]*p[0])))>0)] = 3.0

	def setParameters(self, name, value):            
		if value < self.bounds[name][0]:
			self.parameters[name] = self.bounds[name][0]
		elif value > self.bounds[name][1]:
			self.parameters[name] = self.bounds[name][1]
		else:
			self.parameters[name] = value                

	def setAllParameters(self, parameters):
		for i in parameters.iterkeys():
			if i in self.bounds.keys():
				self.setParameters(i, parameters[i])

	def startExp(self):
		self.varPos = self.parameters['gamma']
		self.varGoal = 0.0
		self.reward_found = False

	def startTrial(self):
		self.current_position = '1a'
		self.previous_position = '1a'
		self.varPos = self.parameters['gamma']
		self.current_action = None				

	def softMax(self, values):
		tmp = np.exp(values*float(self.parameters['beta']))
		if np.isinf(tmp).sum() == 1:
			tmp = np.isinf(tmp)*(1.0-len(tmp)*1e-8) + 1e-8
		p_a = tmp/float(np.sum(tmp))
		if np.sum(p_a==0.0):
			p_a += 1e-8
			p_a = p_a/float(np.sum(p_a))			
		return p_a

	def sampleSoftMax(self, values):
		tmp = np.exp(values*float(self.parameters['beta']))
		if np.isinf(tmp).sum() == 1:
			tmp = np.isinf(tmp)*(1.0-len(tmp)*1e-8) + 1e-8
		p_a = tmp/float(np.sum(tmp))
		if np.sum(p_a==0.0):
			p_a += 1e-8
			p_a = p_a/float(np.sum(p_a))	
		tmp = [np.sum(p_a[0:i]) for i in range(len(p_a))]
		return np.sum(np.array(tmp) < np.random.rand())-1 

	def computeAngle(self, s1, s2):
		x = self.positions[s2][0]-self.positions[s1][0]
		y = self.positions[s2][1]-self.positions[s1][1]
		return (np.arctan2(y, x)+2*np.pi)%(2*np.pi)
		# return np.arctan2(y,x)

	def cdf_multi(self, lower, upper, mu, cov):
	    # lower = [x0,y0]
	    # upper= [x1,y1]
	    # mu = [xm, ym]
	    # cov = sigma
	    stdev = np.sqrt(np.diag(np.eye(2)*cov))
	    lower = (lower-mu)/stdev
	    upper = (upper-mu)/stdev
	    correl = np.zeros(1)
	    infin = 2.0*np.ones(2)
	    error, cdfvalue, inform = scipy.stats.kde.mvn.mvndst(lower,upper,infin,correl)
	    return cdfvalue    

	def fill_PGoal(self):
		for y in xrange(self.n_case):
			for x in xrange(self.n_case):
				self.Pgoal[y,x] = self.cdf_multi(self.grid[y,x],self.grid[y+1,x+1], self.reward_position, self.varGoal)

	def fill_PPos(self):
		for y in xrange(self.n_case):
			for x in xrange(self.n_case):
				self.Ppos[y,x] = self.cdf_multi(self.grid[y,x],self.grid[y+1,x+1], self.positions[self.current_position], self.varPos)

	def computeValue(self, position, state, a, possible):
		# Very tricky : state in [U,Y,I] and a in [0,1,2,3]
		self.current_position = position
		self.fill_PPos()
		self.current_action = a
		ind = self.ind[possible==1]
		self.q_values = np.ones(len(ind))
		if self.reward_found and len(ind) > 1:
			for i in xrange(len(ind)):				
				k = self.transition[(self.previous_position, self.current_position,ind[i])]
				tmp = np.sum(np.sum(self.mask[state][self.current_position][k]*np.atleast_3d(self.Pgoal),1), 0)				
				self.q_values[i] = np.sum(tmp * self.Ppos.flatten())
				# self.q_values[i] = np.sum((np.sum(self.mask[state][self.current_position][k]*np.atleast_3d(self.Pgoal), (0,1)))*(self.Ppos.flatten()))				
		p_a = np.zeros(self.n_action)		
		p_a[ind] = self.softMax(self.q_values)
		return p_a[self.current_action]

	def chooseAction(self, position, state, possible):		
		self.current_position = position
		self.fill_PPos()
		ind = self.ind[possible==1]
		self.q_values = np.ones(len(ind))
		if self.reward_found and len(ind) > 1:			
			for i in xrange(len(ind)):
				k = self.transition[(self.previous_position, self.current_position,ind[i])]
				self.q_values[i] = np.sum((np.sum(self.mask[state][self.current_position][k]*np.atleast_3d(self.Pgoal), (0,1)))*(self.Ppos.flatten()))
		# print "In PI: ", self.q_values				
		self.current_action = ind[self.sampleSoftMax(self.q_values)]
		return actions[self.current_action]		

	def updateValue(self, reward, next_state):
		r = (reward==0)*0.0+(reward==1)*1.0+(reward==-1)*0.0		        
		self.varPos += self.parameters['gamma']
		self.previous_position = self.current_position
		if r:
			self.reward_found = True
			self.varGoal = (1.0-self.parameters['eta'])*self.varGoal + self.parameters['eta']*self.varPos
			self.fill_PGoal()
