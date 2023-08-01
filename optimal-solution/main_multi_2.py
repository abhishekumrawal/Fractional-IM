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
    
"Relabling the nodes so that {1,2} is the greedy best"
mapping = {1:10,2:9,9:2,10:1,3:3,4:4,5:5,6:6,7:7,8:8}
network = nx.relabel_nodes(network, mapping)

"Printing act_prob"
#nx.get_edge_attributes(network,'act_prob')

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
options = {
    #'node_color': 'blue',
    'node_size': 1000,
    #'width': 4,
    'arrowstyle': '-|>',
    'arrowsize': 30,
}

pos=nx.spring_layout(network)
nx.draw(network,pos, font_size = 15, font_family='times', node_color='w',edge_color='#BB0000',width=1,with_labels=True, **options)
nx.draw_networkx_nodes(network,pos, node_color='#A0CBE2', node_size=1000, node_shape='o', alpha=None)        
nx.draw_networkx_edges(network,pos, width = 2, **options)
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

spontaneous_prob = np.zeros(len(network))
prob = np.arange(0,1.05,0.05)

prob1 = []; prob2 =[]; influence2= []
for i in range(len(prob)):
    for j in range(len(prob)):
        prob1.append(prob[i])
        prob2.append(prob[j])
        spontaneous_prob[best_seed_set[0]-1] = prob[i] 
        spontaneous_prob[best_seed_set[1]-1] = prob[j] 
        influence2.append(influence(network,[],diffusion_model,n_sim=n_sim,spontaneous_prob=spontaneous_prob))

output2 = pd.DataFrame()
output2['prob1'] = prob1; output2['prob2'] = prob2; output2['influence'] = influence2

end_time=time.time()
print("run time: " + str((end_time-start_time)) + " sec")

"Approximating the influence with multi-linear extension"
influence_corners = output2[output2.prob1.isin([1,0]) & output2.prob2.isin([1,0])].influence

multilinear_influence2 = []
for i in range(len(prob)):
    for j in range(len(prob)):
        x = np.array([(1-prob[i])*(1-prob[j]), (1-prob[i])*(prob[j]), (prob[i])*(1-prob[j]), (prob[i])*(prob[j])])
        multilinear_influence2.append(np.dot(influence_corners,x))

output2['multilinear_influence'] = multilinear_influence2
output2['absolute_error'] = np.abs(output2.influence - output2.multilinear_influence)
mape = np.mean(np.abs(output2.influence - output2.multilinear_influence))
sde = np.sqrt(np.var(np.abs(output2.influence - output2.multilinear_influence)))
output2.to_csv('im_outputs/multilinear/budget_2/results.csv')

"Reading output2"
output2 = pd.read_csv('im_outputs/multilinear/budget_2/results.csv')

output2.error = output2.influence - output2.multilinear_influence

"Plotting the influence"
plt.rcParams['figure.figsize'] = [3,4]
plt.rcParams.update({'font.size': 11})
width = 3
#%matplotlib
fig = plt.figure()
ax = Axes3D(fig)
surface = ax.plot_trisurf(output2.prob1, output2.prob2, output2.influence, cmap=cm.jet, linewidth=0.1)
ax.set_xlabel('Discount to user 1 ($d_1$)')
ax.set_ylabel('Discount to user 2 ($d_2$)')
ax.set_zlabel('Influence')
#fig.colorbar(surface, shrink=0.5, aspect=5)
fig.savefig('im_outputs/multilinear/budget_2/actual_influence.png',dpi=300, bbox_inches = "tight")

"Plotting the influence"
#plt.rcParams['figure.figsize'] = [3,3]
plt.rcParams.update({'font.size': 11})
width = 3
#%matplotlib
fig = plt.figure()
ax = Axes3D(fig)
surface = ax.plot_trisurf(output2.prob1, output2.prob2, output2.multilinear_influence, cmap=cm.jet, linewidth=0.1)
ax.set_xlabel('Discount to user 1 ($d_1$)')
ax.set_ylabel('Discount to user 2 ($d_2$)')
ax.set_zlabel('Multilinear ext.')
#fig.colorbar(surface, shrink=0.5, aspect=5)
fig.savefig('im_outputs/multilinear/budget_2/multilinear_influence.png',dpi=300, bbox_inches = "tight")


"Plotting the influence eroors"
#plt.rcParams['figure.figsize'] = [3,3]
plt.rcParams.update({'font.size': 11})
width = 3
#%matplotlib
fig = plt.figure()
ax = Axes3D(fig)
surface = ax.plot_trisurf(output2.prob1, output2.prob2, output2.error, cmap=cm.jet, linewidth=0.1)
ax.set_xlabel('Discount to user 1 ($d_1$)')
ax.set_ylabel('Discount to user 2 ($d_2$)')
ax.set_zlabel('Error')
#fig.colorbar(surface, shrink=0.5, aspect=5)
fig.savefig('im_outputs/multilinear/budget_2/influence_error.png',dpi=300, bbox_inches = "tight")

# "Plotting the influence error"
# fig23 = plt.figure(3)
# ax = Axes3D(fig23)
# surface = ax.plot_trisurf(output2.prob1, output2.prob2, output2.absolute_error, cmap=cm.jet, linewidth=0.1)
# ax.set_title('Multilinear Influence Surface')
# fig23.colorbar(surface, shrink=0.5, aspect=5)
# fig23.savefig('im_outputs/multilinear/budget_2/influence_error.png')
