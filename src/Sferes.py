#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Sferes.py

    class for multi-objective optimization
    to interface with sferes2 : see
    http://sferes2.isir.upmc.fr/    

Copyright (c) 2014 Guillaume VIEJO. All rights reserved.
"""

import sys
import os
import mmap
import numpy as np
sys.path.append("../../src")
from Models import *
if os.uname()[1] in ['atlantis', 'paradise']:
    from multiprocessing import Pool, Process
    from pylab import *
    import cPickle as pickle


def unwrap_self_load_data(arg, **kwarg):
    return pareto.loadPooled(*arg, **kwarg)

def unwrap_self_re_test(arg, **kwarg):
    return pareto.poolTest(*arg, **kwarg)

class maxlikelihood():
    """
    Explore Pareto Front from Sferes Optimization
    """
    def __init__(self, directory, speed = 'slowloading'):
        self.directory = directory
        self.data = dict()
        self.models = dict({"VMWM":VMWM(),
                            "Graph":Graph(),
                            "PI":PI()})
        self.p_order = dict({'VMWM':['beta', 'gamma', 'eta'],
                            'Graph':['beta', 'gamma', 'eta'],
                            'PI':['beta','gamma','eta']})
        # self.p_order = dict({'VMWM':['beta', 'gamma', 'eta', 'length']})        
        self.best = dict() # the best parameters set for each mouse        
        self.best_log = dict()
        self.upper_front = dict() # the fitness max at each generation for each mouse and each model 
        self.lower_front = dict() # the fitness min at each geneartion for each mouse and each model 
        self.label = {'late stage':['B28','B61','B137','B155','B163','B166','B150','B62','B76','B84','B139','B154','B74','B86','B152'],
                    'late stage approx':['B18','B20','B143','B78','B141','B144','B165','B82'],
                    'slow learner':['B22','B25','B58','B68','B70','B148','B161']}
        if speed == 'fastloading':
            print "Loading is fast, only last generation is loaded \n Change with speed='slowloading'\n"
            print("Loaded :")
            self.loadData() # fast loading of the last generation
        elif speed == 'slowloading':            
            print "Loading is slow, all the generation are loaded \n Change with speed='fastloading'\n"
            print("Loaded :")
            self.simpleLoadData() # slow loading of all the generation

    def simpleLoadData(self):
        model_in_folders = os.listdir(self.directory)
        if len(model_in_folders) == 0:
            sys.exit("No model found in directory "+self.directory)
        for m in model_in_folders:
            print("\t Modele :"+m)
            self.data[m] = dict()
            lrun = os.listdir(self.directory+"/"+m)
            order = self.p_order[m]
            scale = self.models[m].bounds
            for r in lrun:
                s = "_".join(r.split("_")[2:4])[:-6]
                n = int(r.split("_")[4].split(".")[0])
                print("\t\t Mouse : "+s)
                n = int(r.split("_")[4].split(".")[0])
                if s in self.data[m].keys():
                    self.data[m][s][n] = np.genfromtxt(self.directory+"/"+m+"/"+r)
                else :
                    self.data[m][s] = dict()
                    self.data[m][s][n] = np.genfromtxt(self.directory+"/"+m+"/"+r)                                
                for p in order:
                    self.data[m][s][n][:,order.index(p)+4] = scale[p][0]+self.data[m][s][n][:,order.index(p)+4]*(scale[p][1]-scale[p][0])
        for m in self.data.iterkeys():
            for s in self.data[m].iterkeys():
                tmp = []
                for n in self.data[m][s].iterkeys():
                    tmp.append(np.hstack((np.ones((len(self.data[m][s][n]),1))*n,self.data[m][s][n])))
                self.data[m][s] = np.vstack(tmp)

    def loadData(self):
        model_in_folders = os.listdir(self.directory)
        if len(model_in_folders) == 0:
            sys.exit("No model found in directory "+self.directory)

        pool = Pool(len(model_in_folders))
        tmp = pool.map(unwrap_self_load_data, zip([self]*len(model_in_folders), model_in_folders))
        
        for d in tmp:
            self.data[d.keys()[0]] = d[d.keys()[0]]

    def loadPooled(self, m):
        print("\t Modele :"+m)      
        data = {m:{}}
        list_file = os.listdir(self.directory+"/"+m)
        order = self.p_order[m]
        scale = self.models[m].bounds
        for r in list_file:
            s = r.split("_")[2]
            print("\t\t Mouse : "+s)
            n = int(r.split("_")[4].split(".")[0])
            filename = self.directory+"/"+m+"/"+r            
            nb_ind = int(self.tail(filename, 1)[0].split(" ")[1])
            last_gen = np.array(map(lambda x: x[0:-1].split(" "), self.tail(filename, nb_ind+1))).astype('float')
            if s in data[m].keys():
                data[m][s] = last_gen
            else:
                data[m][s] = {n:last_gen}
            for p in order:
                data[m][s][:,order.index(p)+4] = scale[p][0]+data[m][s][:,order.index(p)+4]*(scale[p][1]-scale[p][0])                    
        return data

    def tail(self, filename, n):
        size = os.path.getsize(filename)
        with open(filename, "rb") as f:
            fm = mmap.mmap(f.fileno(), 0, mmap.MAP_SHARED, mmap.PROT_READ)
            for i in xrange(size-1, -1, -1):
                if fm[i] == '\n':
                    n -= 1
                    if n == -1:
                        break
            return fm[i+1 if i else 0:].splitlines()

    def extractFrontLimits(self):
        for m in self.data.iterkeys():
            self.upper_front[m] = dict()
            self.lower_front[m] = dict()
            for s in self.data[m].iterkeys():                                
                ind = np.where(self.data[m][s][:,1] == 0)[0]
                self.upper_front[m][s] = self.data[m][s][ind,2]
                self.lower_front[m][s] = np.array([np.min(self.data[m][s][self.data[m][s][:,0]==g,2]) for g in np.unique(self.data[m][s][:,0])])
                # tmp = np.array([[np.max(self.data[m][s][self.data[m][s][:,0]==g,2]),np.min(self.data[m][s][self.data[m][s][:,0]==g,2])] for g in np.unique(self.data[m][s][:,0])])
                # self.upper_front[m][s] = tmp[:,0]
                # self.lower_front[m][s] = tmp[:,1]

    def plotEvolution(self, mouse=None):
        rcParams['ytick.labelsize'] = 8
        rcParams['xtick.labelsize'] = 8      
        n_model = len(self.data.keys())          
        fig = figure()         
        for m,i in zip(self.data.iterkeys(), xrange(n_model)):
            ax = fig.add_subplot(n_model,1,i+1)
            for s in self.data[m].iterkeys():                                
                # ax.fill_between(np.unique(self.data[m][s][:,0]), self.upper_front[m][s], self.lower_front[m][s], alpha = 0.2)
                ax.plot(self.upper_front[m][s], label = s, color= 'black')
            ax.set_title(m)            
            # ax.set_legend()

            
    def extractBestLog(self):
        for m in self.data.iterkeys():
            self.best[m] = dict()
            self.best_log[m] = dict()
            for s in self.data[m].iterkeys():
                # best_ind = np.where(self.data[m][s][:,0] == np.max(self.data[m][s][:,0]))[0][0] 
                best_ind = np.argmax(self.data[m][s][:,3])
                self.best[m][s] = {p:self.data[m][s][best_ind,5+self.p_order[m].index(p)] for p in self.p_order[m]}
                self.best_log[m][s] = self.data[m][s][best_ind,3]        
        
        # Winner par model
        subject = self.data[self.data.keys()[0]].keys()
        order = self.data.keys()
        self.winner = {m:[] for m in order}
        for s in subject:
            tmp = []
            for m in order:
                if s in self.best_log[m].keys():
                    tmp.append(self.best_log[m][s])
                else:
                    tmp.append(-10000.0)
            win = np.argmax(tmp)
            self.winner[order[win]].append(s)
            
        # group of best model
        winner = {m:[s.split("_")[0] for s in self.winner[m]] for m in self.winner.keys()}        
        self.group = dict()
        for g in self.label.iterkeys():
            self.group[g] = dict()
            for m in order:
                self.group[g][m] = []
                for s in winner[m]:                    
                    if s in self.label[g]:
                        self.group[g][m].append(s)                            

    def write(self, name):
        with open(name+"_parameters.pickle", 'wb') as f:
            pickle.dump(self.best, f)
        with open(name+"_parameters.txt", 'w') as f:
            f.write("BEST \n")
            # for g in self.group.iterkeys():
            for g in ['late stage', 'late stage approx', 'slow learner']:
                for m in self.group[g].iterkeys():
                    line = "group="+g+"\tmodel="+m+"\t"
                    line += " ".join([s for s in self.group[g][m]])
                    line += "\n"
                    f.write(line)
                f.write("\n")
            f.write("\n")
            for m in self.best.keys():
                f.write(m+"\n")
                tmp = {i.split("_")[0]:i for i in self.best[m].keys()}
                for g in ['late stage', 'late stage approx', 'slow learner']:
                    for ss in self.label[g]:                        
                        if ss in tmp.keys():
                            s = tmp[ss]
                            line=g+"\t mouse="+s.split("_")[0]+"\t"
                            line += " \t ".join([k+"="+str(self.best[m][s][k]) for k in self.p_order[m]])
                            line += "\t\tloglikelihood = "+str(self.best_log[m][s])+"\n"                            
                            f.write(line)
                f.write("\n")

        



    def plotHistBest(self, name):

        best_log = {m:{s.split("_")[0]:self.best_log[m][s] for s in self.best_log[m].iterkeys()} for m in self.best_log.iterkeys()}
        figure()         
        for l,i in zip(self.label.keys(),range(len(self.label.keys()))):
            subplot(len(self.label.keys()),1,i+1)            
            x_pos = np.arange(len(self.label[l]))
            for s,j in zip(self.label[l],x_pos):
                bar(j, best_log['VMWM'][s], 0.1)

        tight_layout()
        grid()
        legend()
        show()
        # savefig(file_name)

class pareto():
    """
    Explore Pareto Front from Sferes Optimization
    """
    def __init__(self, directory):
        self.directory = directory        
        self.data = dict()
        self.models = dict({"VMWM":VMWM()})

        self.p_order = dict({'VMWM':['beta', 'gamma', 'eta']})

        self.m_order = ['VMWM']
        self.colors_m = dict({'VMWM':'r'})
        self.opt = dict()
        self.pareto = dict()
        self.distance = dict()
        self.owa = dict()
        self.tche = dict()
        self.p_test = dict()
        self.mixed = dict()        
        self.indd = dict()
        self.zoom = dict()
        self.timing = dict()
        self.simpleLoadData()

    def simpleLoadData(self):
        model_in_folders = os.listdir(self.directory)
        if len(model_in_folders) == 0:
            sys.exit("No model found in directory "+self.directory)
        for m in model_in_folders:
            print("\t Modele :"+m)
            self.data[m] = dict()
            lrun = os.listdir(self.directory+"/"+m)
            order = self.p_order[m]
            scale = self.models[m].bounds
            for r in lrun:
                s = "_".join(r.split("_")[2:4])[:-6]
                n = int(r.split("_")[4].split(".")[0])
                print("\t\t Mouse : "+s)
                n = int(r.split("_")[4].split(".")[0])
                if s in self.data[m].keys():
                    self.data[m][s][n] = np.genfromtxt(self.directory+"/"+m+"/"+r)
                else :
                    self.data[m][s] = dict()
                    self.data[m][s][n] = np.genfromtxt(self.directory+"/"+m+"/"+r)                                
                for p in order:
                    self.data[m][s][n][:,order.index(p)+4] = scale[p][0]+self.data[m][s][n][:,order.index(p)+4]*(scale[p][1]-scale[p][0])

    def constructParetoFrontier(self):
        for m in self.data.iterkeys():
            self.pareto[m] = dict()
            for s in self.data[m].iterkeys():
                self.pareto[m][s] = dict()   
                tmp={n:self.data[m][s][n][self.data[m][s][n][:,0]==np.max(self.data[m][s][n][:,0])] for n in self.data[m][s].iterkeys()}
                tmp=np.vstack([np.hstack((np.ones((len(tmp[n]),1))*n,tmp[n])) for n in tmp.iterkeys()])
                ind = tmp[:,3] != 0
                tmp = tmp[ind]
                tmp = tmp[tmp[:,3].argsort()][::-1]
                pareto_frontier = [tmp[0]]
                for pair in tmp[1:]:
                    if pair[4] >= pareto_frontier[-1][4]:
                        pareto_frontier.append(pair)
                self.pareto[m][s] = np.array(pareto_frontier)

    def preview(self):
        rcParams['ytick.labelsize'] = 8
        rcParams['xtick.labelsize'] = 8        
        fig_model = figure(figsize = (10,10)) # for each model all subject            

        for m,i in zip(self.pareto.iterkeys(), xrange(len(self.pareto.keys()))):
            ax2 = fig_model.add_subplot(1,1,i+1)
            for s in self.pareto[m].iterkeys():
                # ax2.plot(self.pareto[m][s][:,3], self.pareto[m][s][:,4], "-o", alpha = 1.0, label = s)        
                ax2.plot(self.pareto[m][s][:,3], self.pareto[m][s][:,4], "-o", color = self.colors_m[m], alpha = 1.0)        
                ax2.plot(self.indd['distance'][m][s][3], self.indd['distance'][m][s][4], "*", markersize = 15)                
            ax2.set_title(m)
            ax2.set_xlim(-2000,0)

    def rankDistance(self):
        self.p_test['distance'] = dict()        
        self.indd['distance'] = dict()
        for m in self.pareto.iterkeys():
            self.distance[m] = dict()
            self.p_test['distance'][m] = dict()
            self.indd['distance'][m] = dict()
            for s in self.pareto[m].iterkeys():
                self.distance[m][s] = np.sqrt(np.sum(np.power(self.pareto[m][s][:,3:5]-np.array([-500.0,-8000.0]), 2),1))
                ind_best_point = np.argmin(self.distance[m][s])
                best_point = self.pareto[m][s][ind_best_point,3:5]
                # Saving best individual                        
                best_ind = self.pareto[m][s][ind_best_point]
                self.indd['distance'][m][s] = best_ind                
                self.p_test['distance'][m][s] = dict(zip(self.p_order[m],best_ind[5:]))

    def write(self, name, o):
        with open(name+"_parameters.pickle", 'wb') as f:
            pickle.dump(self.p_test[o], f)
        with open(name+"_parameters.txt", 'w') as f:
            for m in self.p_test[o].keys():
                for s in self.p_test[o][m].keys():
                    line="mouse="+s.split("_")[0]+"\t"
                    line += " \t ".join([k+"="+str(self.p_test[o][m][s][k]) for k in self.p_order[m]])
                    # line += "\t\tloglikelihood = "+str(self.best_log[m][s])+"\n"
                    line += "\n"
                    f.write(line)
        




