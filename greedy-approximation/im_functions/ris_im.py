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

def ris_im(network, budget, diffusion_model, n_sim=100, all_upto_budget=True):
    
    "set a random seed"
    np.random.seed(int(random.uniform(0, 1000000)))
    
    "ris code folder"
    ris_folder = './ris_code_release'
                
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
    # Generating input for RIS, i.e.
        - network.name[1:]: the actual diretced graph
    """
    
    "## Number of nodes and number of edges"
    num_nodes = str(network.number_of_nodes())
    num_edges = str(network.number_of_edges())
    first_row = [num_nodes,num_edges]

    "## Writing to file"
    fstr = ris_folder+"/graphInfo/"+network.name[1:]
    with open(fstr, 'w') as f:
        f.write(' '.join([str(x) for x in first_row]))
        f.write('\n')
        for edge in list(network.edges()):
            f.write(' '.join([str(x) for x in edge]))
            f.write('\n')
    
    "### Running RIS software which is written in C and saving the outputs"
    os.chdir(ris_folder)
    os.system("make")
    start = timeit.default_timer()
    os.system("./subsim -func=format -gname="+network.name[1:]+" -pdist=wc")  
    os.system("./subsim -func=im -gname="+network.name[1:]+" -seedsize="+str(budget)+" -eps=0.01")
    end = timeit.default_timer()
    runtime = end - start
    os.chdir("..")    
    
    "##### Saving runtime info to a text file"
    runtime_info = {'ris':runtime}
    fstr = results_folder_runtime_files+os.sep+'runtime_info_ris.txt'
    with open(fstr, 'w') as f:
        f.write(json.dumps(runtime_info))
    
    "#### Output filename for best seed set"
    out_filename = ris_folder+"/result/seed/seed_"+network.name[1:]+"_subsim_k"+str(budget)+"_wc"
    
    "#### Getting the best seed set and exp influence"
    best_seed_set = [x.split(' ')[0] for x in open(out_filename).readlines()]
    exp_influence = [0 for x in best_seed_set]
    #exp_influence = [x.split(' ')[1] for x in open(out_filename).readlines()]
    best_seed_set = [int(x) for x in best_seed_set]  
    #exp_influence = [float(x) for x in exp_influence]

    best_seed_set += random.sample(set(list(network.nodes)).difference(best_seed_set),budget-len(best_seed_set))
    #exp_influence += [exp_influence[-1] for x in range(budget-len(exp_influence))]
    
    "# Saving the results"
    if all_upto_budget == True:
        results = {'budget':budget, 'diffusion_model':diffusion_model, 'algorithm':'ris', 'n_sim':n_sim, \
                   'best_seed_set':[[None]] + [best_seed_set[:k+1] for k,_ in enumerate(best_seed_set)],\
                       'network_name':network.name, 'exp_influence':[0] + exp_influence}
            
        fstr = results_folder_pickle_files+os.sep+'output_ris__%i__.pkl'%(budget)
        with open(fstr,'wb') as f:
            pickle.dump(results, f)
                
        logging.info('The final solution is as follows.')
        logging.info(str([[None]] + [best_seed_set[:k+1] for k,_ in enumerate(best_seed_set)]))
        logging.info(str([0] + exp_influence))
        logging.info('Total time taken by RIS-IM is '+' '+str(round(runtime,2))+' seconds.') 
                        
        return [[None]] + [best_seed_set[:k+1] for k,_ in enumerate(best_seed_set)], [0] + exp_influence, runtime
    
    else:
        results = {'budget':budget, 'diffusion_model':diffusion_model, 'algorithm':'ris', 'n_sim':n_sim, \
                   'best_seed_set':best_seed_set, 'network_name':network.name, 'exp_influence':exp_influence[-1]}
            
        fstr = results_folder_pickle_files+os.sep+'output_ris__%i__.pkl'%(budget)
        with open(fstr,'wb') as f:
            pickle.dump(results, f)
        
        logging.info('The final solution is as follows.')
        logging.info(str(best_seed_set)) 
        logging.info(exp_influence[-1])
        logging.info('Total time taken by RIS-IM is '+' '+str(round(runtime,2))+' seconds.') 
                
        return best_seed_set, exp_influence[-1], runtime
