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
		self.pos_guiding = np.array(['1a','1b','I','10','II','8','III'])
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
		self.label = {'late stage':['B28','B61','B137','B155','B163','B166','B150','B62','B76','B84','B139','B154','B74','B86','B152'],
                    'late stage approx':['B18','B20','B143','B78','B141','B144','B165','B82'],
                    'slow learner':['B22','B25','B58','B68','B70','B148','B161']}

	def move(self, action):
		# 0 : Forward, 1 : Left, 2 : Backward, 3 : Right
		self.world.moveALaLouche(self.actions.index(action))

	def guidage(self):		
		self.model.startTrial()
		for i in xrange(6):			
			self.model.computeValue(self.pos_guiding[i], self.state_guiding[i],self.action_guiding[i],self.possible_guiding[i])
			self.model.updateValue(self.reward_guiding[i], self.state_guiding[i+1])		

	def sferes(self, data):
		np.seterr(all='ignore')
		nb_point = data['info']['nb_point']
		nb_trial = data['info']['nb_trial']
		loglike = np.zeros(nb_point)
		self.model.startExp()
		self.world.startingPos()
		for i in xrange(nb_trial):
		# for i in xrange(1):
			self.model.startTrial()
			tmp = 0.0
			for j in xrange(len(data[i]['action'])):								
				pa = self.model.computeValue(data[i]['pos'][j], data[i]['state'][j], data[i]['action'][j], data[i]['possible'][j])								
				# if i==23:
				# 	if np.sum(data[i]['possible'][j])>1:											
				# 		print self.model.q_values				
				# if i == 3 and np.sum(data[i]['possible'][j])>1:
				# 	print i,j									
				# 	a = self.model.mask['I']['1b'][1]*np.atleast_3d(self.model.Pgoal)
				# 	# for k in xrange(100):
				# 	# 	# print k, np.sum(a[:,:,k])
				# 	# 	print k, self.model.Ppos.flatten()[k]					
				# 	if 
					# print self.model.q_values

				self.model.updateValue(data[i]['reward'][j], data[i]['state'][j+1])
				loglike[data[i]['ind'][j]] = np.log(pa)
				tmp += np.log(pa)				
			if data[i]['reward'][-1] == 0:
				self.guidage()
			
			# print i, tmp
			print tmp
		# for i in xrange(900):
		# 	print self.model.xy[i], self.model.Pgoal.flatten()[i]


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
			print n
			self.model.startExp()
			for i in xrange(nb_trials):
				# print i, self.model.varGoal
				self.model.startTrial()	
				self.reward_found = False		
				self.world.startingPos()
				state = self.pos_to_state[self.world.mousePos]				
				for j in xrange(self.nb_steps_max):					
					possible = self.world.readPathwaysALaLouche()
					action = self.model.chooseAction(self.world.mousePos, state, possible)
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

		data = {m:{s.split("_")[0]:data[m][s] for s in data[m].keys()} for m in data.keys()}
		colors = {'Graph':'green','VMWM':'blue','PI':'gray'}
		
		# for g in self.label.keys():		
		for g in ['late stage']:
			fig = figure(figsize=(14,10))
			for s,i in zip(self.label[g], xrange(len(self.label[g]))):
				ax = fig.add_subplot(3,5,i+1)				
				for m in data.keys():
					mean_time = np.mean(data[m][s], 0)
					std_time = np.std(data[m][s], 0)					
					ax.plot(np.arange(len(mean_time)), mean_time, '-', label=m, color = colors[m], linewidth=3)		
					ax.fill_between(np.arange(len(mean_time)), mean_time-(std_time/2.), mean_time+(std_time/2.), alpha=0.24, color=colors[m])					
					ax.set_title(s)
					ax.set_ylabel("Latency", fontsize=8)
					ax.set_xlabel("Trials", fontsize=8)
				ax.plot(latency[s.split("_")[0]], 'o-', color = 'red', label='Mouse')
				if i == 0:
					legend(bbox_to_anchor=(1.2, 1.05))	
			tight_layout()
			savefig(filename+"_group_"+g.replace(" ", "_")+"_test.pdf")


		# joining = [filename+"_group_"+g.replace(" ", "_")+"_test.pdf" for g in ['late stage', 'late stage approx', 'slow learner']]
		# os.system("pdftk "+" ".join(joining)+" cat output "+"../test/SFERES9_group_test_all_models.pdf")
			# os.system("evince "+filename+"_group_"+g+"_test.pdf")
		# # probleme car tout les essais ne font pas la meme taille
		# mean_all = np.zeros(84)
		# std_all = np.zeros(84)
		# mean_mice = np.zeros(84)
		# std_mice = np.zeros(84)
		# for i in xrange(84):
		# 	tmp = []
		# 	tmp2 = []
		# 	for s in data.iterkeys():
		# 		if data[s].shape[1] > i:
		# 			tmp.append(data[s][:,i])
		# 			tmp2.append(latency[s.split("_")[0]][i])
		# 	mean_all[i] = np.mean(tmp)
		# 	std_all[i] = np.std(tmp)
		# 	mean_mice[i] = np.mean(tmp2)
		# 	std_mice[i] = np.std(tmp2)

		# ax = fig.add_subplot(3,5,15)		
		# ax.plot(np.arange(len(mean_all)), mean_all, 'o-', label='Model')
		# ax.fill_between(np.arange(len(mean_all)), mean_all-(std_all/2.), mean_all+(std_all/2.), alpha=0.5)
		# ax.plot(np.arange(len(mean_mice)), mean_mice, '-', color = 'red', label='Mouse')
		# ax.fill_between(np.arange(len(mean_all)), mean_all-(std_all/2.), mean_all+(std_all/2.), alpha=0.5)




					
