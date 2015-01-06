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
		self.actions = ['F', 'L', 'R', 'U']
		# SPECIFIC TO THE MAZE
		self.world = World("TYM")
		self.current_pos = None
		self.reward_found = None
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

	def sferes(mouse_number, mouse_data):		
		pass

	def test(self, nb_exp):
		for n in xrange(nb_exp):
			self.model.startExp()
			self.reward_found = False		
			self.world.startingPos()
			while self.reward_found:
				state = self.pos_to_state[self.world.mousePos]				
				action = self.model.chooseAction(state)
				self.move(action)
				new_state = self.pos_to_state[self.world.mousePos]
				reward = self.world.readRew()
				self.model.updateValue(reward, next_state)





	