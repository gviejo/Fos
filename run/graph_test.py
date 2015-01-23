#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
import sys
sys.path.append("../src")
from Wrap import TYMaze
from Models import *
from matplotlib import *
from pylab import *




parameters = {'beta': 10.405800000000001,
 'eta': 0.9,
 'gamma': 0.9}
 
model = Graph(parameters)

wrap = TYMaze(model)

data = wrap.test(parameters, 10, 40)

plot(np.mean(data,0))

show()