#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
import sys
sys.path.append("../src")
from Models import VMWM

states = np.array(['I', 'Y', 'U'])
sequences = states[np.random.randint(0, 3, 100)]

parameters = dict({'length':3,
					'beta':1.0,
					'eta':0.9,
					'gamma':0.3})

model = VMWM(parameters)
model.startExp()

for i in xrange(len(sequences)-1):
	model.chooseAction(sequences[i])
	r = np.random.randint(2)
	model.updateValue(r, sequences[i+1])	



