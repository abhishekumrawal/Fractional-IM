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
from scipy.special import comb
from itertools import combinations

def optimal_im(network, budget, diffusion_model, n_sim, influence_dict = {}, spontaneous_prob = []):
    
    """
    Brute Force Algorithm for finding the best seed set of a user-specified budget
    
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
        
        combs_to_try = []
        for i in range(int(comb(len(nodes), budget))):
            combs_to_try.append(list(list(combinations(nodes, budget))[i]))
        
        influence = np.zeros(len(combs_to_try))
        
        for i in range(len(combs_to_try)):
            
            for j in range(n_sim):
                spontaneously_infected = []
                
                if len(spontaneous_prob) != 0:
                    
                    for m in range(len(nodes)):
                        if np.random.rand() < spontaneous_prob[m]:
                            spontaneously_infected.append(nodes[m])
                                  
                spontaneously_infected_plus_ith_comb = \
                    list(set(spontaneously_infected + combs_to_try[i]))
                
                if tuple(spontaneously_infected_plus_ith_comb) in influence_dict.keys():
                    influence[i] = influence_dict[tuple(spontaneously_infected_plus_ith_comb)]
                
                if diffusion_model == "independent_cascade":
                    layers = independent_cascade(network, spontaneously_infected_plus_ith_comb)  
                
                elif diffusion_model == "linear_threshold":
                    layers = linear_threshold(network, spontaneously_infected_plus_ith_comb)    
                
                for k in range(len(layers)):
                    influence[i] = influence[i] + len(layers[k])
                
                influence_dict[tuple(spontaneously_infected_plus_ith_comb)] = influence[i]
                    
            influence[i] = influence[i]/n_sim
        
        max_influence = np.max(influence)    
        best_seed_set = combs_to_try[np.argmax(influence)]

    
    return best_seed_set, max_influence, influence_dict