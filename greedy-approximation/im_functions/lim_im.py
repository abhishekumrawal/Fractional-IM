#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 21:30:04 2019

@author: abhishek.umrawal
"""

"importing required built-in modules"
import numpy as np
import random
import os
import pickle
import timeit
import json
import logging

def lim_im(network, budget, diffusion_model, n_sim=100, all_upto_budget=True):
    
    "set a random seed"
    np.random.seed(int(random.uniform(0, 1000000)))
    
    "lim code folder"
    lim_folder = './lim_code_release'
                
    "set budget as len(network.nodes) if the budget > len(network.nodes)"
    budget = min(budget,len(network.nodes))
        
    "creating pickle files folder within the results folder"
    results_folder_pickle_files = 'results'+os.sep+'results'+network.name+os.sep+'pickle_files'
    if not os.path.exists(results_folder_pickle_files):
        os.makedirs(results_folder_pickle_files) 
     
    "creating log files folder within the results folder"    
    results_folder_log_files = 'results'+os.sep+'results'+network.name+os.sep+'log_files'
    if not os.path.exists(results_folder_log_files):
        os.makedirs(results_folder_log_files) 
    
    "creating runtime files folder within the results folder"    
    results_folder_runtime_files = 'results'+os.sep+'results'+network.name+os.sep+'runtime_files'
    if not os.path.exists(results_folder_runtime_files):
        os.makedirs(results_folder_runtime_files) 
             
    """
    # Generating input for lim, i.e.
        - network.name[1:]: the actual diretced graph
    """
   
    "## Number of nodes and number of edges"
    num_nodes = str(network.number_of_nodes())
    num_edges = str(network.number_of_edges())
    
    "## Weighted edges"
    weighted_edges = []
    for edge in network.edges():
        weighted_edges.append([edge[0],edge[1],network[edge[0]][edge[1]]['act_prob']])

    "## Writing to file"
    fstr = lim_folder+"/data/"+network.name[1:]+".txt"
    with open(fstr, 'w') as f:
        f.write(num_nodes)
        f.write('\n')
        f.write(num_edges)
        f.write('\n')
        for weighted_edge in weighted_edges:
            f.write(' '.join([str(x) for x in weighted_edge]))
            f.write('\n')
    
    "### Running lim software which is written in C and saving the outputs"
    os.chdir(lim_folder+'/src')
    os.system("make")
    start = timeit.default_timer()
    os.system("./run main.cpp")    
    end = timeit.default_timer()
    runtime = end - start
    os.chdir("..") 
    os.chdir("..") 
    
    "##### Saving runtime info to a text file"
    runtime_info = {'lim':runtime}
    fstr = results_folder_runtime_files+os.sep+'runtime_info_lim.txt'
    with open(fstr, 'w') as f:
        f.write(json.dumps(runtime_info))
    
    "#### Output filename for best seed set"
    out_filename_allocation = lim_folder+"/allocation/"+network.name[1:]+".txt_cimm_eps=0.500000_group_0_new"
    
    "#### Output filename for exp influences"
    out_filename_exp = lim_folder+"/result/"+network.name[1:]+".txt_cimm_eps=0.500000_group_0_new"
    
    "#### Output filename for time"
    out_filename_time = lim_folder+"/time/"+network.name[1:]+".txt_cimm_eps=0.500000_group_0_new"
    
    "#### Getting the best seed sets (allocations) and exp influence"
    best_seed_sets = [[float(0) for x in range(network.number_of_nodes())]]
    
    for best_seed_set in open(out_filename_allocation).readlines():
        best_seed_sets.append([float(x) for x in best_seed_set.split(' ')[:-1]])
    
    "#### Getting the exp influences"
    exp_influence = [x.split(' ')[7] for x in open(out_filename_exp).readlines()]
    exp_influence = [float(x) for x in exp_influence]
    
    "#### Getting the runtimes (cumulative in seconds)"
    run_times = [x.split(' ')[7] for x in open(out_filename_time).readlines()]
    run_times = np.cumsum([0]+[float(x) for x in run_times])

    ##### Saving all runtimes to a text file"
    fstr = results_folder_runtime_files+os.sep+'runtime_info_lim_all.txt'
    with open(fstr, 'w') as f:
        for val in run_times:
            f.write(str(val))
            f.write('\n')
        
    #best_seed_set += random.sample(set(list(network.nodes)).difference(best_seed_set),budget-len(best_seed_set))
    #exp_influence += [exp_influence[-1] for x in range(budget-len(exp_influence))]
    
    "# Saving the results"
    if all_upto_budget == True:
        results = {'budget':budget, 'diffusion_model':diffusion_model, 'algorithm':'lim', 'n_sim':n_sim, \
                   'best_seed_set': best_seed_sets,\
                       'network_name':network.name, 'exp_influence':[0] + exp_influence}
            
        fstr = results_folder_pickle_files+os.sep+'output_lim__%i__.pkl'%(budget)
        with open(fstr,'wb') as f:
            pickle.dump(results, f)
                
        logging.info('The estimated exp influences are as follows.')
        logging.info(str([0] + exp_influence))
        logging.info('Total time taken by lim-IM is '+' '+str(round(runtime,2))+' seconds.') 
                        
        return [[None]] + [best_seed_set[:k+1] for k,_ in enumerate(best_seed_set)], [0] + exp_influence, runtime
    
    else:
        results = {'budget':budget, 'diffusion_model':diffusion_model, 'algorithm':'lim', 'n_sim':n_sim, \
                   'best_seed_set':best_seed_sets[-1], 'network_name':network.name, 'exp_influence':exp_influence[-1]}
            
        fstr = results_folder_pickle_files+os.sep+'output_lim__%i__.pkl'%(budget)
        with open(fstr,'wb') as f:
            pickle.dump(results, f)
        
        logging.info('The estimated exp influence is as follows.')
        logging.info(str(exp_influence[-1]))
        logging.info('Total time taken by lim-IM is '+' '+str(round(runtime,2))+' seconds.') 
                
        return best_seed_set, exp_influence[-1], runtime
