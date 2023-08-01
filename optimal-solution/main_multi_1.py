#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  1 22:45:43 2018

@author: abhishek.umrawal
"""

"Importing necessary modules"

import numpy as np
import pandas as pd
import time
from scipy.special import comb
import math
import os
os.system("ls -l")

import networkx as nx
import matplotlib.pyplot as plt
%matplotlib
plt.rcParams['figure.figsize'] = [15,15]

import warnings
warnings.filterwarnings('ignore')

import sys
sys.path.append('im_functions')

from influence import influence

from degree_im import degree_im
from wdegree_im import wdegree_im
from greedy_im import greedy_im
from optimal_im import optimal_im


"Defining an empty networkx digraph object"
network = nx.DiGraph()

"Reading edges from a text file"
edge_list = pd.read_table('network_files/small_network1.txt', sep=' ', names=('from','to'))
#edge_list = pd.read_table("facebook_network.txt", sep=' ', names=('from','to'))
#edge_list = edge_list[0:101]

"Defining random edge weights"
weight = np.random.rand(len(edge_list))/1.5

"Adding edges to the network"
for i in range(len(edge_list)):
    network.add_edge(edge_list['from'][i],edge_list['to'][i],act_prob=weight[i])

"Visualizing out-degree"
nodes = list(nx.nodes(network))
out_degree_list = []
for node in nodes:
    out_degree_list.append(network.out_degree(node))

out_degree_df=pd.DataFrame()
out_degree_df['nodes']=nodes
out_degree_df['out_degree']=out_degree_list
sorted_out_degree_df = out_degree_df.sort_values(by='out_degree',ascending=False)

"Visualizing the network"
pos=nx.spring_layout(network)
nx.draw(network,pos,node_color='#A0CBE2',edge_color='#BB0000',width=2,edge_cmap=plt.cm.Blues,with_labels=True)
plt.savefig('im_outputs/multilinear/budget_1/network.png', dpi=500, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format=None,transparent=False, bbox_inches=None, pad_inches=0.1) 
        
"Marketing budget"
budget = 1 #int(np.floor(len(network)/3))

"Defining vector of random spontaneous adoption probabilities"
spontaneous_prob = np.random.rand(len(network))/2

"Declaring the diffusion model"
diffusion_model = "independent_cascade"

"Number of simulations"
n_sim = 10000

"Runnig Greedy Algorithm with Spontaneous Adoption to find the best seed set of size 2"
best_seed_set, max_influence, influence_dict = greedy_im(network,budget,diffusion_model,n_sim=n_sim,influence_dict={},spontaneous_prob=[])

"Finding influence of 'one' best seed (no interventions) with varying spontaneous infection probabilities"
spontaneous_prob = np.zeros(len(network))
prob1 = np.arange(0,1.01,0.01)

influence1= []
for i in range(len(prob1)):
    spontaneous_prob[best_seed_set[0]-1] = prob1[i]
    influence1.append(influence(network,[],diffusion_model,n_sim=n_sim,spontaneous_prob=spontaneous_prob))

output1 = pd.DataFrame()
output1['prob1'] = prob1; output1['influence'] = influence1

"Approximating the influence with multi-linear extension"
influence_corners = output1[output1.prob1.isin([1,0])].influence

multilinear_influence1 = []
for i in range(len(prob1)):
    x = np.array([(1-prob1[i]),prob1[i]])
    multilinear_influence1.append(np.dot(influence_corners,x))

output1['multilinear_influence'] = multilinear_influence1
output1['absolute_error'] = np.abs(output1.influence - output1.multilinear_influence)
mape = np.mean(np.abs(output1.influence - output1.multilinear_influence))
sde = np.sqrt(np.var(np.abs(output1.influence - output1.multilinear_influence)))
output1.to_csv('im_outputs/multilinear/budget_1/results.csv')

"Plotting influence"
fig1 = plt.figure()
plt.plot(prob1,influence1)
plt.xlabel('initial probability')
plt.ylabel('expected influence')
fig1.savefig('im_outputs/multilinear/budget_1/actual_influence.png')

"Plotting multilinear influence"
fig1 = plt.figure()
plt.scatter(prob1,multilinear_influence1)
plt.xlabel('initial probability')
plt.ylabel('expected influence')
fig1.savefig('im_outputs/multilinear/budget_1/multilinear_influence.png')

"Plotting influence error"
fig1 = plt.figure()
plt.plot(prob1,output1.absolute_error)
plt.xlabel('initial probability')
plt.ylabel('influence error')
fig1.savefig('im_outputs/multilinear/budget_1/influence_error.png')