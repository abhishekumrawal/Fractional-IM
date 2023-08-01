#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 18:29:52 2018

@author: abhishek.umrawal
"""
from independent_cascade import independent_cascade
from linear_threshold import linear_threshold
import networkx as nx
import numpy as np

def greedy_im(network, budget, diffusion_model, n_sim, influence_dict = {}, spontaneous_prob = []):
    
    """
    Greedy Algorithm for finding the best seed set of a user-specified budget
    
    Inputs:
        - network is a networkx object
        - budget is the user-specified marketing budget which represents the no.
          of individuals to be given the freebies
        - diffusion model is either "independent_cascade" or "linear_threshold"
        - n_sim is the no. of simulations to be perfomed to estimate the expected influence
        - influence_dict is a dictionary of seed set and expected influence pairs to save computation
        - spontaneous_prob is a vector of spontaneous adoption probabiities for
          each node
          
    Outputs:
        - best_seed_set is a subset of the set of nodes with the cardinality  
          same as the budget such that it maximizes the spread of marketing
        - max_influence is the value of maximum influence
  
    """
    nodes = list(nx.nodes(network))
    max_influence = []
    best_seed_set = []
    
    if budget == 0:
        
        influence = 0
        
        for j in range(n_sim):
            spontaneously_infected = []
            
            if len(spontaneous_prob) != 0:
                
                for m in range(len(nodes)):
                    
                    if np.random.rand() < spontaneous_prob[m]:
                        spontaneously_infected.append(nodes[m])
                                      
            if diffusion_model == "independent_cascade":
                layers = independent_cascade(network, spontaneously_infected)  
                    
            elif diffusion_model == "linear_threshold":
                layers = linear_threshold(network, spontaneously_infected)    
                    
            for k in range(len(layers)):
                influence = influence + len(layers[k])
                        
        influence = influence/n_sim
        max_influence.append(influence)
        best_seed_set.append(None)
    
    else:
        
        for l in range(budget):
            nodes_to_try = list(set(nodes)-set(best_seed_set))
            influence = np.zeros(len(nodes_to_try))
            
            for i in range(len(nodes_to_try)):
                   
                spontaneously_infected = []

                if len(spontaneous_prob) != 0:

                    for m in range(len(nodes_to_try)):
                        if np.random.rand() < spontaneous_prob[m]:
                            spontaneously_infected.append(nodes[m])

                best_seed_set_plus_ith_node = \
                    list(set(spontaneously_infected + best_seed_set + [nodes_to_try[i]]))
                
                if tuple(best_seed_set_plus_ith_node) in influence_dict.keys():
                    influence[i] = influence_dict[tuple(best_seed_set_plus_ith_node)]
                
                else: 
                    for j in range(n_sim):

                        if diffusion_model == "independent_cascade":
                            layers = independent_cascade(network, best_seed_set_plus_ith_node)  

                        elif diffusion_model == "linear_threshold":
                            layers = linear_threshold(network, best_seed_set_plus_ith_node)    

                        for k in range(len(layers)):
                            influence[i] = influence[i] + len(layers[k])

                    influence[i] = influence[i]/n_sim
                    influence_dict[tuple(best_seed_set_plus_ith_node)] = influence[i]
                    
            max_influence.append(np.max(influence))    
            best_seed_set.append(nodes_to_try[np.argmax(influence)])
    
    return best_seed_set, max_influence, influence_dict