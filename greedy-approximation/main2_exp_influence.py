if __name__ == '__main__':

    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-
    """
    Created on Wed May 22 21:49:19 2019
    
    @author: abhishek.umrawal
    """
    import numpy as np
    import pickle as pickle
    import os as os
    import networkx as nx
    
    from multiprocessing import Pool
    import itertools
    import logging
    import timeit
    
    import pandas as pd
    import tikzplotlib
    import matplotlib.pyplot as plt
    plt.rcParams['figure.figsize'] = [8,6]
    plt.rcParams.update({'font.size': 12})
    
    from im_functions.weighted_network import weighted_network
    from im_functions.true_influence import true_influence
    
    "inputs"
    algorithms = ['lim']  #"one at a time"
    chosen_algorithm_indices = [0]
    budget =  20  
    interval = 1                   
    name_id = '_epinions'
    is_graph_already_weighted = False
    weighting_scheme = 'wc'
    num_procs = 33
    graph_type = 'directed'
    probs = [0.20, 0.40, 0.60, 0.80]
    
    "start timer"
    start = timeit.default_timer()
    
    "reading the network"
    if 'florentine' in name_id:
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
    
    "adding weights if the network is unweighted"
    if is_graph_already_weighted == False:
        network = weighted_network(network,method = weighting_scheme)
    
    "creating log files folder within the results folder" 
    results_folder_log_files = 'results'+os.sep+'results'+network.name+os.sep+'log_files'
    if not os.path.exists(results_folder_log_files):
        os.makedirs(results_folder_log_files)
    
    "remove the log file from previous runs"
    if os.path.exists(results_folder_log_files+os.sep+'exp_influence.log'):
        os.remove(results_folder_log_files+os.sep+'exp_influence.log')
        
    "removing exisiting log handlers"
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    
    "set up logging to file"
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(name)-6s %(levelname)-6s %(message)s',
                        datefmt='%m-%d %H:%M',
                        filename=results_folder_log_files+os.sep+'exp_influence.log',
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
    logging.info('I am running main2_exp_influence.py for '+network.name[1:]+' network.')
    logging.info('The network is '+graph_type+'.')
    logging.info('The network has '+str(len(network.nodes))+' nodes and '+str(len(network.edges))+' edges.')
    logging.info('I am using '+weighting_scheme+' weighting scheme.')
    
    "results folder path"
    results_folder = 'results'+os.sep+'results'+name_id
    
    "reading seed set and exp_influence pair from algo results"
    budget = min(budget,len(network.nodes))
    filename = 'output_'+algorithms[0]+'__'+str(budget)+'__.pkl'
    filename_with_path = results_folder+os.sep+'pickle_files'+os.sep+filename
    if os.path.exists(filename_with_path):
        with open(filename_with_path, 'rb') as f:
            results_algo = pickle.load(f) 
    
    #"creating a list of expected influences for integer budgets"
    #exp_influences_list_int_budgets = results_algo['exp_influence']
    
    "diffusion model and number of simulations"
    diffusion_model = results_algo['diffusion_model']
    n_sim = results_algo['n_sim']
    
    "calculating the expected influence and saving in the exp_influences_dict"
    logging.info('I am using '+str(n_sim)+' Monte-Carlo simulations.')    
    
    "for CELF++ and RIS"
    if algorithms[0] in ['celfpp','ris']:    
        "set to list for seed sets in unique_seed_sets"
        unique_seed_sets = [[]]+results_algo['best_seed_set'][1:]
    
        "creating input combinations for integer budgets"
        inputs = []
        for i,seed_set in enumerate(unique_seed_sets):
            tmp = [ [network], [seed_set], [diffusion_model], [n_sim], [[]], [name_id]]
            inputs.append(itertools.product( *tmp ))
        inputs_all = itertools.chain(*inputs)
        inputs_all = [tuple(i) for i in inputs_all]
        
        "parallelization"
        pool = Pool(processes=num_procs)
        exp_influences_list_int_budgets_raw = list(pool.map(true_influence, inputs_all))
        pool.close()
        pool.join()  
        
        "creating a list of expected influences for integer budgets"
        exp_influences_list_int_budgets = []
        for [seed_set,exp_influence] in exp_influences_list_int_budgets_raw:
             exp_influences_list_int_budgets.append(exp_influence)
    
        "creating combinations of spontaneous probabilities"
        spontaneous_probs = []
        for seed_set in unique_seed_sets[1:]:
            for prob in probs:
                spontaneous_prob = list(np.zeros(len(network)))
                spontaneous_prob[seed_set[-1]-1] = prob
                spontaneous_probs.append(spontaneous_prob)
                
        "creating all input combinations"
        inputs = []
        for i,seed_set in enumerate(unique_seed_sets):
            tmp = [ [network], [seed_set], [diffusion_model], [n_sim], spontaneous_probs[i*len(probs):(i+1)*len(probs)], [name_id]]
            inputs.append(itertools.product( *tmp ))
        inputs_all = itertools.chain(*inputs)
        inputs_all = [tuple(i) for i in inputs_all]
        
        "parallelization"
        pool = Pool(processes=num_procs)
        exp_influences_list_frac_budgets_raw = list(pool.map(true_influence, inputs_all))
        pool.close()
        pool.join()  
                
        "creating a list of expected influences for fractional budgets"
        exp_influences_list_frac_budgets = []
        for [seed_set,exp_influence] in exp_influences_list_frac_budgets_raw:
             exp_influences_list_frac_budgets.append(exp_influence)
        
        "combining all expected influences and sorting"
        exp_influences_list = sorted(exp_influences_list_int_budgets + exp_influences_list_frac_budgets)
    
    if algorithms[0] in ['lim']:  
        "for LIM"
        spontaneous_probs = results_algo['best_seed_set']
        "creating all input combinations"
        inputs = []
        for spontaneous_prob in spontaneous_probs:
            tmp = [ [network], [[]], [diffusion_model], [n_sim], [spontaneous_prob], [name_id]]
            inputs.append(itertools.product( *tmp ))
        inputs_all = itertools.chain(*inputs)
        inputs_all = [tuple(i) for i in inputs_all]
        
        "parallelization"
        pool = Pool(processes=num_procs)
        exp_influences_list_frac_budgets_raw = list(pool.map(true_influence, inputs_all))
        print(exp_influences_list_frac_budgets_raw)
        pool.close()
        pool.join()  
    
        "creating a list of expected influences for fractional budgets"
        exp_influences_list = []
        for [seed_set,exp_influence] in exp_influences_list_frac_budgets_raw:
            exp_influences_list.append(exp_influence)
    
    "creating a list of all budgets"
    budgets = [0]
    for int_budget in range(budget):
        for frac_budget in probs+[1]:
            budgets.append(int_budget+frac_budget)    
    
    "savings the results"
    filename = 'output_'+algorithms[0]+'__'+str(budget)+'__frac__.pkl'
    filename_with_path = results_folder+os.sep+'pickle_files'+os.sep+filename
    results = {'budget':budgets,'exp_influence':exp_influences_list}
    fstr = filename_with_path
    with open(fstr,'wb') as f:
        pickle.dump(results, f)
    
    logging.info('I finished running main2_exp_influence.py for '+network.name[1:]+' network.')
    logging.info('The network is '+graph_type+'.')
    logging.info('The network has '+str(len(network.nodes))+' nodes and '+str(len(network.edges))+' edges.')
    logging.info('I used '+weighting_scheme+ ' weighting scheme.')
    logging.info('I used '+str(n_sim)+' Monte-Carlo simulations.')
    
    "end timer"
    end = timeit.default_timer()
    logging.info('Total time taken by main_exp_influence.py is '+str(round(end - start,2))+' seconds.') 
    
    
    "plotting"
    if not os.path.exists(results_folder+os.sep+'plots'):
        os.makedirs(results_folder+os.sep+'plots')
    
    "creating a dataframe using budgets and exp_influences_list"
    output_df = pd.DataFrame()
    output_df['exp_influence'] = exp_influences_list
    output_df.index = budgets           
    ax = output_df.plot(lw=2,marker='o')
    ax.set_xlabel('Budget')
    ax.set_ylabel('Expected Influence')
    ax.grid(linestyle='-', linewidth=1)
    ax.set_title('Influence maximization for '+name_id[1:]+' network')
    legend_entries = algorithms
    ax.legend(legend_entries)
    tikzplotlib.clean_figure()
    tikzplotlib.save(results_folder+os.sep+'plots'+'/sim_'+algorithms[0]+'_'+name_id+'_plot_'+str(budget)+'.tex')
