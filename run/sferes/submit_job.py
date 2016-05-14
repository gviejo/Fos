#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""                                                                                   
submit_job.py                                                                         
                                                                                      
to launch sferes job on cluster                                                       
                                 
version du modele de guillaume

Copyright (c) 2013 Guillaume Viejo. All rights reserved.                              
"""

import sys,os
from optparse import OptionParser

# -----------------------------------                                                 
# ARGUMENT MANAGER                                                                    
# -----------------------------------                                                 
if not sys.argv[1:]:
   sys.stdout.write("Sorry: you must specify at least 1 argument")
   sys.stdout.write("More help avalaible with -h or --help option")
   sys.exit(0)
parser = OptionParser()
parser.add_option("-m", "--model", action="store", help="The name of the model to optimize", default="VMWM")
parser.add_option("-t", "--time", action="store", help="Time of execution", default=False)
# parser.add_option("-d", "--data", action="store", help="The data to fit", default="52trials")
parser.add_option("-n", "--n_run", action="store", help="Number of run per subject", default=1)
(options, args) = parser.parse_args()
# -----------------------------------                                                 


# ------------------------------------
# CREATE DIRECTORY RESULTS
# ------------------------------------
os.system("rm -r /home/benedicte.babayan/Fos/run/sferes/"+options.model)
os.system("mkdir /home/benedicte.babayan/Fos/run/sferes/"+options.model)

# -----------------------------------                                                 
# GENERATE SCRIPTS                                                               
# -----------------------------------                                                 

all_mice = [i.split(".")[0] for i in os.listdir("../../data/")]
few_mice = ['B28','B61','B137','B155','B163','B166','B150','B62','B76','B84','B139','B154','B74','B86','B152']
list_mice = [i for i in all_mice if i.split("_")[0] in few_mice]
# list_mice = [i.split(".")[0] for i in os.listdir("datacut/")]

for s in xrange(len(list_mice)):	
    for i in xrange(int(options.n_run)):
    	# # ----------------------------------
    	# # Generate Bash Scripts
    	# # ----------------------------------
        filename = "submit_"+options.model+"_"+list_mice[s]+"_"+str(i)+".sh"
        f = open(filename, "w")
        f.writelines("#!/bin/sh\n")
        f.writelines("#PBS -N sferes_"+options.model+"_"+list_mice[s]+"_"+str(i)+"\n")
        f.writelines("#PBS -o /home/Fos/log/sferes_"+options.model+"_"+options.time+"_"+list_mice[s]+"_"+str(i)+".out\n")
        f.writelines("#PBS -b /home/Fos/log/sferes_"+options.model+"_"+options.time+"_"+list_mice[s]+"_"+str(i)+".err\n")
        f.writelines("#PBS -m abe\n")
        f.writelines("#PBS -M guillaume.viejo@gmail.com\n")
        f.writelines("#PBS -l walltime="+options.time+"\n")
        f.writelines("#PBS -l nodes=1:ppn=8"+"\n")                                                  
        f.writelines("#PBS -d /home/benedicte.babayan/Fos\n")
        # f.writelines("#PBS -v PYTHONPATH=/home/benedicte.babayan/lib/python/lib/python\n") 
        f.writelines("/home/benedicte.babayan/sferes2/build/debug/exp/"+options.model+"/"+options.model+" --subject="+list_mice[s]+" --n_run="+str(i)+"\n")        
        f.close()
        #-----------------------------------

        # ----------------------------------                                                  
        # GENERATE PYTHON SCRIPTS                                                             
        # ----------------------------------                                                  
        pythonfile = "/home/benedicte.babayan/sferes2/exp/"+options.model+"/sferes_"+options.model+"_"+list_mice[s]+".py"
        pf = open(pythonfile, 'w')                                                          
        pf.write("""#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from optparse import OptionParser                                                   
sys.path.append("/home/Fos/src/")
from Wrap import TYMaze
from Models import VMWM
import cPickle as pickle
parser = OptionParser()                                                             
parser.add_option("-b", "--beta", action="store", type = 'float')
parser.add_option("-g", "--gamma", action="store", type = 'float')
parser.add_option("-e", "--eta", action="store", type = 'float')
parser.add_option("-n", "--length", action="store", type='float')
(options, args) = parser.parse_args()
parameters = vars(options)
model = VMWM()
for p in parameters.iterkeys():                                                                                                              
    if parameters[p] is not None:
        parameters[p] = model.bounds[p][0]+parameters[p]*(model.bounds[p][1]-model.bounds[p][0])
parameters['length'] = 3
model = VMWM(parameters)
opt = TYMaze(model)\n""")
        pf.writelines("""with open('/home/benedicte.babayan/Fos/data/"""+list_mice[s]+""".pickle') as f:\n""")
        pf.writelines("""	data = pickle.load(f)\n""")
        pf.writelines("""llh = opt.sferes(data)\n""")
        pf.writelines("""print llh""")       
        pf.close()
        # ------------------------------------                                                                                            
        # SUBMIT                                                                                                                          
        # ------------------------------------                                                                                            
        # os.system("chmod +x "+filename)
        # os.system("qsub "+filename)
        # os.system("rm "+filename)



