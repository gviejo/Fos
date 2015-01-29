import numpy as np
import sys
sys.path.append("../src")
from Wrap import TYMaze
from Models import *
from world import World
from matplotlib import *
from pylab import *
import cPickle as pickle



with open("mask.pickle") as handle:
	data = pickle.load(handle)

parameters = {'beta': 20.405800000000001,
 'eta': 0.9,
 'gamma': 0.96420399990357952}
 
model = PI(parameters)

wrap = TYMaze(model)
world = wrap.world
