#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys, os

all_mice = ['B163_52trials',
			'B76_68trials',
			'B144_84trials',
			'B165_84trials',
			'B78_84trials',
			'B139_68trials',
			'B143_68trials',
			'B58_68trials',
			'B25_68trials',
			'B22_68trials',
			'B18_68trials',
			'B82_84trials',
			'B154_68trials',
			'B166_52trials',
			'B161_84trials',
			'B20_68trials',
			'B155_52trials',
			'B68_68trials',
			'B152_84trials',
			'B70_68trials',
			'B141_84trials',
			'B84_68trials',
			'B86_84trials',
			'B62_68trials',
			'B137_52trials',
			'B148_84trials',
			'B61_52trials',
			'B28_52trials',
			'B74_84trials',
			'B150_52trials']

few_mice = ['B28','B61','B137','B155','B163','B166','B150','B62','B76','B84','B139','B154','B74','\
B86','B152']

list_mice = [i[0:-6] for i in all_mice if i.split("_")[0] in few_mice]

for i in list_mice:
	for j in xrange(10):
		command_line = "/home/viejo/sferes2/build/debug/exp/PI/PI --subject="+i+" --n_run="+str(j)+" --n_trials="+i[-2:]
		print command_line