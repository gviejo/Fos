# Fos

Project Sequence-based learning in a double T-Maze 

Packages
-----------------------
works with :
- python 2.7.3
- numpy 1.8.0
- scipy 0.13.0
- matplotlib 1.4
- cPickle 1.71

Folders/Files
-----------------------
- *run/* : All the main files to be called
  - *run/ccall* : The c++ version of Path Integration for called by Sferes2 for optimisation
    - *run/ccall/main.cpp* : The test call to compute log-likelihood
    - *run/ccall/normal.cpp* : A sub-routine for computing multivariate normal probabilities
    - *run/ccall/pi_call.cpp* : The main function for evaluating log-likelihood of path-integration
    - *run/ccall/submit_cjob.py* : To submit c++ job on the cluster for individual parameters optimisation
  - *run/sferes* : The main files for analysing the results of sferes2 optimisation
    - *run/sferes/analysis.py* : To load and and plot multi objective results from Sferes 2 optimisation.
    - *run/sferes/analysis2.py* : Construct and analyse a pareto frontier of two dimensions for special cases
    - *run/sferes/convert_to_pickle.py* : Convert original text data to pickle data for gain of speed processing
    - *run/sferes/logTest.py* : Test files for log-optimisation
    - *run/sferes/main.py* : Test files for parameters
    - *run/sferes/sferes_model_name-of-the-mouse.py* : Examples of files generated when lauching job on the cluster and called by sferes2 to compute log-likelihood
    - *run/sferes/submit_job* : To submit python job on the cluster for individual parameters optimisation
    - *run/sferes/testSferes.py* : Take pickle parameters set and test it to generate learning curves.
    
- *src/* : All the python classes 
  - *src/Models.py* : All the models that can be instanciated for testing parameters set
  - *src/Sferes.py* : Classes for analysing sferes2 output. Generates parameters set according to maximum log-likelihood
  - *src/Wrap.py* : Contains the information about the world for testing parameters set. 
  - *src/world.py* : Benoit's code. Allow the agent to move in the world. Called by Wrap.py
