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
		self.bounds = dict({"length":[0,3],
							"gamma":[0.0, 1.0],
							"beta":[0.0, 100.0],
							"eta":[0.0, 1.0]})
		self.delta = 0.0		
		# Values initialization		
		self.critic	= np.ones(self.n_possible)*0.1
		self.actor = np.ones((self.n_possible, self.n_action))*0.1
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
		self.critic	= np.ones(self.n_possible)*0.1
		self.actor = np.ones((self.n_possible, self.n_action))*0.1

	def softMax(self, values):
		tmp = np.exp(values*float(self.parameters['beta']))
		if np.isinf(tmp).sum():
			self.p_a = np.isinf(self.p_a)*0.9999996+0.0000001            
		else :
			self.p_a = tmp/np.sum(tmp)           
		return tmp/float(np.sum(tmp))

	def sampleSoftMax(self, values):
		tmp = np.exp(values*float(self.parameters['beta']))
		if np.isinf(tmp).sum():
			tmp = np.isinf(self.p_a)*0.9999996+0.0000001
		else :
			tmp = tmp/float(np.sum(tmp))
		tmp = [np.sum(tmp[0:i]) for i in range(len(tmp))]
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
		self.last_actions[1:] = self.last_actions[:-1]
		self.last_actions[0] = actions[self.current_action]		
		# Delta update
		n_s = self.convert[next_state+"".join(self.last_actions)]
		self.delta = r + self.parameters['gamma']*self.critic[n_s] - self.critic[self.current_state]
		self.critic[self.current_state] += (self.parameters['eta']*self.delta)
		self.actor[self.current_state, self.current_action] += (self.parameters['eta']*self.delta)

