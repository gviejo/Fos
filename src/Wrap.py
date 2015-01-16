#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Class of Maze that wraps around the model.
For guidage, reward delivery, maze structure, ...

"""

import numpy as np
import os,sys
from world import World
from matplotlib import *
if os.uname()[1] in ['atlantis', 'paradise']:
   from pylab import *    


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
		self.state_guiding = np.array(['U','I','Y','I','Y','I','U'])
		self.possible_guiding = np.array([[1,0,0,0],[1,0,1,0],[0,1,1,1],[1,0,1,0],[0,1,1,1],[1,0,1,0],[0,0,1,0]])
		self.action_guiding = np.array([0,0,1,0,3,0])
		self.reward_guiding = np.array([0,0,0,0,0,1])
		self.pos_to_state = dict({	'1a':'U',
									'1b':'I',
									'I' :'Y',
									'10':'I',
									'II':'Y',
									'8' :'I',
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
			self.model.computeValue(self.state_guiding[i],self.action_guiding[i],self.possible_guiding[i])
			self.model.updateValue(self.reward_guiding[i], self.state_guiding[i+1])

	def sferes(self, data, epsilon):
		np.seterr(all='ignore')
		nb_point = data['info']['nb_point']
		nb_trial = data['info']['nb_trial']
		loglike = np.zeros(nb_point)
		self.model.startExp()
		self.world.startingPos()
		for i in xrange(nb_trial):
			self.model.startTrial()
			# biais = np.exp(-epsilon*(float(len(data[i]['action']))-6.0))
			for j in xrange(len(data[i]['action'])):				
				pa = self.model.computeValue(data[i]['state'][j], data[i]['action'][j], data[i]['possible'][j])								
				self.model.updateValue(data[i]['reward'][j], data[i]['state'][j+1])
				# loglike[data[i]['ind'][j]] = np.log(pa)*biais
				# loglike[data[i]['ind'][j]] = np.log(pa)*(1.0/float(len(data[i]['action'])-5))
				loglike[data[i]['ind'][j]] = np.log(pa)

			if data[i]['reward'][-1] == 0:
				self.guidage()
		llh = np.sum(loglike)
		if llh==0 or np.isnan(llh) or np.isinf(llh):
			return -100000
		else:
			return np.sum(loglike)
		
	def sferes2(self, data, latency):
		np.seterr(all='ignore')
		nb_point = data['info']['nb_point']
		nb_trial = data['info']['nb_trial']
		loglike = np.zeros(nb_point)
		self.model.startExp()
		self.world.startingPos()
		for i in xrange(nb_trial):
			self.model.startTrial()
			for j in xrange(len(data[i]['action'])):
				pa = self.model.computeValue(data[i]['state'][j], data[i]['action'][j], data[i]['possible'][j])
				self.model.updateValue(data[i]['reward'][j], data[i]['state'][j+1])
				loglike[data[i]['ind'][j]] = np.log(pa)
			if data[i]['reward'][-1] == 0:
				self.guidage()
		llh = np.sum(loglike)
		# Least squares test
		data = self.test(self.model.parameters, 20, nb_trial)
		lrs = -np.sum(np.power(latency-np.mean(data,0), 2))

		if llh==0 or np.isnan(llh) or np.isinf(llh):
			return -100000, lrs
		elif np.isnan(lrs) or np.isinf(lrs):
			return llh, -100000
		else:
			return llh, lrs

	def test(self, parameters, nb_exp, nb_trials):
		self.model.__init__(parameters)
		data = np.zeros((nb_exp, nb_trials))
		for n in xrange(nb_exp):
			self.model.startExp()
			for i in xrange(nb_trials):
				self.model.startTrial()	
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
				data[n,i] = j
				if not self.reward_found:
					self.guidage()
		return data*1.05

	def plot(self, data, mouse, parameters, latency, file_name):
		figure() # for each model all subject            
		mean_time = np.mean(data, 0)
		std_time = np.std(data, 0)
		plot(np.arange(len(mean_time)), mean_time, 'o-', label='Model')		
		fill_between(np.arange(len(mean_time)), mean_time-(std_time/2.), mean_time+(std_time/2.), alpha=0.5)
		plot(latency, 'o-', color = 'red', label='Mouse')
		ylim(0,61)
		desc = ", ".join([k+"="+str(np.round(parameters[k],2)) for k in parameters.keys()])
		ylabel("Latency")
		title(mouse+" | "+desc, fontsize = 11)
		grid()
		legend()
		savefig(file_name)

	def plotall(self, data, latency, filename):
		rcParams['xtick.labelsize'] = 8
		rcParams['ytick.labelsize'] = 8
		fig = figure(figsize=(12,10))

		for s,i in zip(data.keys(), xrange(len(data.keys()))):
			mean_time = np.mean(data[s], 0)
			std_time = np.std(data[s], 0)
			ax = fig.add_subplot(4,4,i+1)
			ax.plot(np.arange(len(mean_time)), mean_time, 'o-', label='Model')		
			ax.fill_between(np.arange(len(mean_time)), mean_time-(std_time/2.), mean_time+(std_time/2.), alpha=0.5)
			ax.plot(latency[s.split("_")[0]], 'o-', color = 'red', label='Mouse')
			ax.set_title(s)
		tight_layout()
		savefig(filename)

					
