#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 19:17:14 2019

@author: abhishek.umrawal
"""

"debugging functionality"    
# clear : to clear the console
# reset: to clear all variables
# %run filename.py to run from console

DEBUG_MODE = False
import pdb
if DEBUG_MODE:
    pdb.set_trace()
else:
    pdb.set_trace = lambda: 1

"important debugging tips"    
# if don't specify pdb.set_trace() here then do python -m pdb <filename.py> from terminal
# type c to continue running until next breakpoint()
# type a to see arguments with their values for the current function
# type print(<var_name>) to a print any variable created so far or
# type p <var_name> for regular printing
# type pp <var_name> for pretty printing
# type !<var_name> to create a new variable
# run
# type q to quit -- all variables made until then would be saved 
# type pdb.set_trace = lambda: 1 to neutralize pdb for the entire session

"importing required built-in modules"
import itertools
import networkx as nx
import timeit
from networkx.classes.function import density
import os
import pickle
import shutil
import logging

"importing required user-defined modules"
from im_functions.non_adaptive_im import non_adaptive_im

"start time"
start = timeit.default_timer()

#################################################################################################################################
"user inputs for non_adaptive_im"
weighting_schemes = ['wc']                          # 'tv' 'rn'
algorithms = ['celfpp']                             # ['optimal','celfpp','ccelfpp1','ris'] 
#algorithms = ['heuristic']                         # ['heuristic]
heuristics = ['']                                   # ['']  
#heuristics = ['degree','degdiscount']              # ['degree','ivgreedy','degdiscount']
max_budgets = [3]#[x for x in range(1,21)]                                 # just provide the max seed set size
diffusion_models = ['independent_cascade']        # 'linear_threshold'
community_methods = ['louvain']                     # infomap, greedy_modularity, label_propagation, girvan_newman
n_sim = 1000                                       # no. of diffusions
communities = []                                    # pre-defined communities
community_size_threshold = 0.01                     # merging communities with nodes less than threshold % of all nodes
is_graph_already_weighted = False                   # False if the graph is not already weighted
graph_type = 'directed'
#################################################################################################################################

#################################################################################################################################
"user input for this main file"
name_id = '_epinions'                             # network to work with and identifier to save the results
#################################################################################################################################

"reading the network"
if 0:#'florentine' in name_id:
    network = nx.florentine_families_graph()
    network.name = name_id
    network = network.to_directed()

elif 'florida' in name_id:
    network=nx.read_weighted_edgelist('network_data/florida_network.csv',delimiter=',',create_using= nx.DiGraph(),nodetype=int)
    nx.set_edge_attributes(network,values= nx.get_edge_attributes(network,'weight'), name='act_prob')
    network.name = name_id

elif graph_type == 'directed':
    network = nx.read_edgelist("network_data/"+name_id[1:]+"_network.txt",create_using=nx.DiGraph(), nodetype = int)
    network.name = name_id
    
elif graph_type == 'undirected':
    network = nx.read_edgelist("network_data/"+name_id[1:]+"_network.txt",create_using=nx.Graph(), nodetype = int)
    network.name = name_id 
    network = network.to_directed()

"relabeling the nodes as positive integers viz. 1,2,..."
network = nx.convert_node_labels_to_integers(network,first_label=1)

"looking up the n_sim -- unused block"
if 0:
    try: 
        with open('network_data/n_sim_data.pkl', 'rb') as f:
             n_sims = pickle.load(f)       
        
        n_sim = n_sims[name_id]
    except:
        pass

"create a list of all parameter lists, then use product"
tmp = [ [network], weighting_schemes, algorithms, heuristics,max_budgets, diffusion_models ]
tmp += [ [n_sim], [name_id] ]
tmp += [ community_methods, [communities], [community_size_threshold], [is_graph_already_weighted] ]
inputs = itertools.product( *tmp )
inputs = [tuple(i) for i in inputs]

"call non_adaptive_im for all input combinations -- no parallelization"
"python doesn't allow daemon processes to have children"
"there is already parallelization in greedy_im and hence in cgreedy_im"
for inpt in inputs:  
    
    "creating log files folder within the results folder" 
    results_folder_log_files = 'results'+os.sep+'results'+network.name+os.sep+'log_files'
    if not os.path.exists(results_folder_log_files):
        os.makedirs(results_folder_log_files)
    
    "remove the log file from previous runs"
    if os.path.exists(results_folder_log_files+os.sep+'log_'+str(inpt[2])+'_'+str(inpt[3])+'.log'):
        os.remove(results_folder_log_files+os.sep+'log_'+str(inpt[2])+'_'+str(inpt[3])+'.log')
        
    "removing exisiting log handlers"
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    
    "set up logging to file"
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(name)-6s %(levelname)-6s %(message)s',
                        datefmt='%m-%d %H:%M',
                        filename=results_folder_log_files+os.sep+str(inpt[2])+'_'+str(inpt[3])+'.log',
                        filemode='w')
    
    "define a Handler which writes INFO messages or higher to the sys.stderr"
    console = logging.StreamHandler()
    console.setLevel(logging.INFO) 
    
    "set a format which is simpler for console use"
    formatter = logging.Formatter('%(name)-6s: %(levelname)-6s %(message)s')
    
    "tell the handler to use this format"
    console.setFormatter(formatter)
    
    "add the handler to the root logger"
    logging.getLogger().addHandler(console)
    
    "Now, we can log to the root logger, or any other logger. First the root..."
    logging.info('I am running '+inpt[2]+'_im'+' upto budget '+str(inpt[4])+' for '+name_id[1:]+' network.')
    logging.info('The network is '+graph_type+'.')
    logging.info('The network has '+str(len(network.nodes))+' nodes and '+str(len(network.edges))+' edges.')
    logging.info('The network density is '+str(density(network))+'.')
    logging.info('I am using '+community_methods[0]+' community detection method.')
    logging.info('I am using '+str(n_sim)+' Monte-Carlo simulations.')
    logging.info('I am using '+weighting_schemes[0]+ ' weighting scheme.')
        
    "calling non_adaptive_im"
    non_adaptive_im(inpt)
    
    logging.info('I finished running '+inpt[2]+'_im'+' upto budget '+str(inpt[4])+' for '+name_id[1:]+' network.')
    logging.info('The network is '+graph_type+'.')
    logging.info('The network has '+str(len(network.nodes))+' nodes and '+str(len(network.edges))+' edges.')
    logging.info('I used '+community_methods[0]+' community detection method.')
    logging.info('I used '+str(n_sim)+' Monte-Carlo simulations.')
    logging.info('I used '+weighting_schemes[0]+ ' weighting scheme.')

    
if 'ccelfpp1' in algorithms:    
    "moving community wise results into a folder named all_community_results within results+network.name folder"
    if not os.path.exists('./results/results'+network.name+os.sep+'all_community_results'):
        os.makedirs('./results/results'+network.name+os.sep+'all_community_results')
    for folder_name in sorted(os.listdir('./results')):
        if network.name+'_community' in folder_name:
            #print(folder_name)
            shutil.move('./results/'+folder_name,'./results/results'+network.name+os.sep+'all_community_results')

"end time"
end = timeit.default_timer()

"time taken"
logging.info('Total time taken is ' + str(round(end - start,4)) + ' seconds.')    
