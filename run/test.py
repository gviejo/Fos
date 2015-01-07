#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
import sys
sys.path.append("../src")
from Wrap import TYMaze
from Models import VMWM

parameters = dict({'length':2,
					'beta':5.0,
					'eta':0.1,
					'gamma':0.9})

model = VMWM(parameters)

wrap = TYMaze(model)

wrap.test(1)



