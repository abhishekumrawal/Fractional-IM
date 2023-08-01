#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 13 04:14:43 2018

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
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

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
budget = 2 #int(np.floor(len(network)/3))

"Defining vector of random spontaneous adoption probabilities"
spontaneous_prob = np.random.rand(len(network))/2

"Declaring the diffusion model"
diffusion_model = "independent_cascade"

"Number of simulations"
n_sim = 10000

"Runnig Greedy Algorithm with Spontaneous Adoption to find the best seed set of size 2"
best_seed_set, max_influence, influence_dict = greedy_im(network,budget,diffusion_model,n_sim=n_sim,influence_dict={},spontaneous_prob=[])

"Finding influence of the best seed set (no interventions) with varying spontaneous infection probabilities"
start_time=time.time()

"Approximating the influence with multi-linear extension"

"Calculating the influence corners"
spontaneous_prob = np.zeros(len(network))
prob = np.arange(0,1.01)
influence_corners= []
for i in range(len(prob)):
    for j in range(len(prob)):
        spontaneous_prob[best_seed_set[0]-1] = prob[i] 
        spontaneous_prob[best_seed_set[1]-1] = prob[j] 
        influence_corners.append(influence(network,[],diffusion_model,n_sim,spontaneous_prob))

end_time=time.time()
print("run time: " + str((end_time-start_time)) + " sec")

"Calculating the influence for the interior points"
prob = np.arange(0,1.001,.001)
prob1 = []; prob2 =[]; influence2 = []
for i in range(len(prob)):
    for j in range(len(prob)):
        prob1.append(prob[i])
        prob2.append(prob[j])
        x = np.array([(1-prob[i])*(1-prob[j]), (1-prob[i])*(prob[j]), (prob[i])*(1-prob[j]), (prob[i])*(prob[j])])
        influence2.append(np.dot(influence_corners,x))

output2 = pd.DataFrame()
output2['prob1'] = prob1; output2['prob2'] = prob2; 
output2['sum_prob'] = output2['prob1'] + output2['prob2']
output2['influence'] = influence2
output2.to_csv('im_outputs/solution/budget_2/output.csv')

budgets = [0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0]
hyperplanes2 = []
for i in range(len(budgets)):
    hyperplanes2.append(output2[output2.sum_prob == budgets[i]])

"Plotting and saving figure - Influence Hyperplane when x1 + x2 + x3 = k"
for i in range(len(budgets)-2):
    %matplotlib
    fig = plt.figure(i+1)
    plt.plot(hyperplanes2[i+1].prob1,hyperplanes2[i+1].influence)
    plt.xlabel('x1')
    plt.ylabel('influence')
    plt.title('Influence Hyperplane when x1 + x2 = '+str(budgets[i+1]))
    fig.savefig('im_outputs/solution/budget_2/influence_hyperplane2-'+str(budgets[i+1])+'.png')
    plt.pause(0.00001)
    plt.show()

"Optimal solutions"
optimal_solutions2 = []
for i in range(len(hyperplanes2)):
    optimal_solutions2.append(hyperplanes2[i][hyperplanes2[i]['influence'] == max(hyperplanes2[i]['influence'])])
optimal_solutions2 = pd.concat(optimal_solutions2)

%matplotlib
fig = plt.figure(i+1)
plt.plot(optimal_solutions2.prob1,optimal_solutions2.prob2)
plt.xlabel('x1')
plt.ylabel('x2')
plt.title('Optimal Solutions for different Budgets')
fig.savefig('im_outputs/solution/budget_2/optimal_solutions.png')
