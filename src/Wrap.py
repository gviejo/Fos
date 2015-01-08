#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Class of Maze that wraps around the model.
For guidage, reward delivery, maze structure, ...

"""

import numpy as np
import sys, os
from world import World

class TYMaze():
	""" from B. Babayan thesis
	"""

	def __init__(self, model):
		self.model = model
		self.states = ['I', 'Y', 'U']
		self.actions = ['F', 'L', 'U', 'R']
		# SPECIFIC TO THE MAZE
		self.world = World("TYM")
		self.current_pos = None
		self.reward_found = None
		self.nb_trials_max = 52
		self.nb_steps_max = 57
		self.state_guiding = np.array([2,0,1,0,1,0,2])
		self.action_guiding = np.array([0,0,1,0,3,0])
		self.reward_guiding = np.array([0,0,0,0,0,1])
		self.pos_to_state = dict({	'1a':'U',
									'1b':'I',
									'I' :'Y',
									'10':'I',
									'II':'Y',
									'8' :'U',
								   'III':'U',
								   	'9a':'U',
								   	'9b':'I',
								   	'2' :'I',
								   	'V' :'Y',
								   	'3b':'I',
								   	'3a':'U',
								   	'4' :'I',
								   	'IV':'U'})

	def move(self, action):
		# 0 : Forward, 1 : Left, 2 : Backward, 3 : Right
		self.world.moveALaLouche(self.actions.index(action))

	def guidage(self):
		for i in xrange(6):
			self.model.current_state = self.state_guiding[i]
			self.model.current_action = self.action_guiding[i]
			self.model.updateValue(self.reward_guiding[i], self.states[self.state_guiding[i+1]])

	def sferes(self, data):
		nb_point = data['info']['nb_point']
		nb_trial = data['info']['nb_trial']
		loglike = np.zeros(nb_point)
		self.model.startExp()
		self.world
		for i in xrange(nb_trial):
			for j in xrange(len(data[i]['action'])):
				pa = self.model.computeValue(data[i]['state'][j], data[i]['action'][j], data[i]['possible'][j])
				
				# print np.log(pa)
				
				self.model.updateValue(data[i]['reward'][j], data[i]['state'][j+1])
				loglike[data[i]['ind'][j]] = np.log(pa)

			if data[i]['reward'][-1] == 0:
				self.guidage()
		llh = np.sum(loglike)
		if int(llh)==0:
			return -100000
		else:
			return np.sum(loglike)
		
	def test(self, nb_exp):
		for n in xrange(nb_exp):
			self.model.startExp()
			for i in xrange(self.nb_trials_max):			
				self.reward_found = False		
				self.world.startingPos()
				state = self.pos_to_state[self.world.mousePos]				
				for j in xrange(self.nb_steps_max):								
					possible = self.world.readPathwaysALaLouche()									
					action = self.model.chooseAction(state, possible)					
					self.move(action)
					state = self.pos_to_state[self.world.mousePos]
					reward = self.world.readRew()
					self.model.updateValue(reward, state)
					if reward:
						self.reward_found = True
						break
				print j
				if not self.reward_found:
					self.guidage()
					

					
