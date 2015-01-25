#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Models.py

model-free reinforcement + working memory = VMWM
"""


import numpy as np
import sys, os

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

	def computeValue(self, state, a, possible):
		# Very tricky : state in [U,Y,I] and a in [0,1,2,3]
		self.current_state = self.convert[state+"".join(self.last_actions)]
		self.current_action = a
		ind = self.ind[possible==1]
		q_values = self.actor[self.current_state][possible==1]
		p_a = np.zeros(self.n_action)		
		p_a[ind] = self.softMax(q_values)			
		return p_a[self.current_action]

	def chooseAction(self, state, possible):
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

	def computeValue(self, state, a, possible):
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

	def chooseAction(self, state, possible):
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
		self.transition = {('10', 'I', 1): [39, 51],
							('10', 'I', 3): [39, 30],
							('10', 'II', 1): [9, 18],
							('10', 'II', 3): [9, 57],
							('1a', '1b', 0): [0, 0],
							('1b', 'I', 1): [0, 9],
							('1b', 'I', 3): [0, 51],
							('2', 'I', 1): [21, 30],
							('2', 'I', 3): [21, 9],
							('2', 'V', 1): [51, 3],
							('2', 'V', 3): [51, 42],
							('3a', '3b', 0): [12, 12],
							('3b', 'V', 1): [12, 21],
							('3b', 'V', 3): [12, 3],
							('4', 'V', 1): [33, 42],
							('4', 'V', 3): [33, 21],
							('8', 'II', 1): [27, 39],
							('8', 'II', 3): [27, 18],
							('9a', '9b', 0): [48, 48],
							('9b', 'II', 1): [48, 57],
							('9b', 'II', 3): [48, 39],
							('I', '10', 0): [9, 9],
							('I', '1b', 0): [30, 30],
							('I', '2', 0): [51, 51],
							('II', '10', 0): [39, 39],
							('II', '8', 0): [57, 57],
							('II', '9b', 0): [18, 18],
							('III', '8', 0): [27, 27],
							('IV', '4', 0): [33, 33],
							('V', '2', 0): [21, 21],
							('V', '3b', 0): [42, 42],
							('V', '4', 0): [3, 3]}
		self.alpha = np.linspace(0, 2*np.pi, 60)
		self.position = {'10': [-0.19, 0.608],
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
		#Various Init

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
		pass

	def startTrial(self):
		pass

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

	def computeValue(self, state, a, possible):
		# Very tricky : state in [U,Y,I] and a in [0,1,2,3]
		self.current_state = state
		self.current_action = a
		ind = self.ind[possible==1]
		q_values = self.actor[self.current_state][possible==1]
		p_a = np.zeros(self.n_action)		
		p_a[ind] = self.softMax(q_values)			
		return p_a[self.current_action]

	def chooseAction(self, state, possible):
		self.current_state = self.convert[state+"".join(self.last_actions)]
		ind = self.ind[possible==1]
		q_values = self.actor[self.current_state][possible==1]		
		self.current_action = ind[self.sampleSoftMax(q_values)]
		return actions[self.current_action]		

	def updateValue(self, reward, next_state):
		r = (reward==0)*0.0+(reward==1)*1.0+(reward==-1)*0.0		        
