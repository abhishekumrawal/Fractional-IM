#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 19:17:14 2019

@author: abhishek.umrawal
"""

"Importing necessary modules"
import numpy as np
import pandas as pd
import time
#from scipy.special import comb
#import math
import os
os.system("ls -l")

import networkx as nx
import matplotlib.pyplot as plt
#%matplotlib
#%matplotlib inline
plt.rcParams['figure.figsize'] = [3,3]
plt.rcParams.update({'font.size': 22})
from mpl_toolkits.mplot3d import Axes3D
#from matplotlib import cm

import warnings
warnings.filterwarnings('ignore')

import sys
sys.path.append('im_functions')

from influence import influence

#from degree_im import degree_im
#from wdegree_im import wdegree_im
#from greedy_im import greedy_im
from optimal_im import optimal_im


"Defining an empty networkx digraph object"
network = nx.DiGraph()

"Reading edges from a text file"
edge_list = pd.read_table('network_files/small_network3.txt', sep=' ', names=('from','to'))
#edge_list = pd.read_table("facebook_network.txt", sep=' ', names=('from','to'))
#edge_list = edge_list[0:101]

"Defining random edge weights"
#weight = np.random.rand(len(edge_list))/1.5
#weight = np.array([0.5, 0.7, 0.7, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5])
weight = np.array([0.6, 0.5, 0.7, 0.4, 0.1, 0.7, 0.5, 0.5, 0.4, 0.5, 0.8, 0.5])

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
#nx.draw(network,pos,node_color='#A0CBE2',edge_color='#BB0000',width=1,edge_cmap=plt.cm.Blues)
        
"Marketing budget"
budget = 3 #int(np.floor(len(network)/3))

"Defining vector of random spontaneous adoption probabilities"
spontaneous_prob = np.random.rand(len(network))/2

"Declaring the diffusion model"
diffusion_model = "independent_cascade"

"Number of simulations"
n_sim = 1000

"Runnig Greedy Algorithm with Spontaneous Adoption to find the best seed set of size 3"
#best_seed_set, max_influence, influence_dict = greedy_im(network,budget,diffusion_model,n_sim=n_sim,influence_dict={},spontaneous_prob=[])

best_seed_set, max_influence, influence_dict = optimal_im(network,budget,diffusion_model,n_sim=n_sim,influence_dict={},spontaneous_prob=[])


#best_seed_set = [1,2,3];

"Finding influence of the best seed set (no interventions) with varying spontaneous infection probabilities"
start_time=time.time()

"Approximating the influence with multi-linear extension"

"Calculating the influence corners"
spontaneous_prob = np.zeros(len(network))
prob = np.arange(0,1.01)
prob1 = []; prob2 = []; prob3 = []; influence_corners= []
for i in range(len(prob)):
    for j in range(len(prob)):
        for k in range(len(prob)):
            prob1.append(prob[i])
            prob2.append(prob[j])
            prob3.append(prob[k])
            spontaneous_prob[best_seed_set[0]-1] = prob[i] 
            spontaneous_prob[best_seed_set[1]-1] = prob[j] 
            spontaneous_prob[best_seed_set[2]-1] = prob[k] 
            influence_corners.append(influence(network,[],diffusion_model,n_sim,spontaneous_prob))
            
influence_corners_df = pd.DataFrame()
influence_corners_df['prob1'] = prob1; influence_corners_df['prob2'] = prob2; influence_corners_df['prob3'] = prob3;  
influence_corners_df['sum_prob'] = influence_corners_df['prob1'] + influence_corners_df['prob2'] + influence_corners_df['prob3']
influence_corners_df['influence'] = influence_corners          
print(influence_corners_df)

"Calculating the influence for the interior points"
prob = np.arange(0,1.01,.01)
prob1 = []; prob2 = []; prob3 = []; influence3 = []
for i in range(len(prob)):
    for j in range(len(prob)):
        for k in range(len(prob)):
            prob1.append(prob[i])
            prob2.append(prob[j])
            prob3.append(prob[k])
            x = np.array([(1-prob[i])*(1-prob[j])*(1-prob[k]), (1-prob[i])*(1-prob[j])*(prob[k]), 
                          (1-prob[i])*(prob[j])*(1-prob[k]), (1-prob[i])*(prob[j])*(prob[k]),
                          (prob[i])*(1-prob[j])*(1-prob[k]), (prob[i])*(1-prob[j])*(prob[k]),
                          (prob[i])*(prob[j])*(1-prob[k]), (prob[i])*(prob[j])*(prob[k])])
            influence3.append(np.dot(influence_corners,x))

output3 = pd.DataFrame()
output3['prob1'] = prob1; output3['prob2'] = prob2; output3['prob3'] = prob3;  
output3['sum_prob'] = output3['prob1'] + output3['prob2'] + output3['prob3']
output3['influence'] = influence3
output3.to_csv('im_outputs/solution/budget_3/output.csv')

output3= pd.read_csv('im_outputs/solution/budget_3/output.csv')

#budgets = [0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.5, 2.6, 2.8, 3.0]
#budgets = [0, 0.1, 0.5, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.51, 1.52, 1.53, 1.54, 1.55, 1.56, 1.57, 1.58, 1.59, 1.6, 1.8, 1.85, 1.9, 1.91, 1.92, 1.95, 1.99, 2.0, 2.1, 2.5, 2.9, 3.0]
#budgets = np.round(np.linspace(0,3,50),2)
#budgets = [0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.5, 2.6, 2.8, 3.0]
#budgets = np.sort(np.append(budgets,[1,1.48,1.49,1.50,1.51,1.52,1.53,1.54,1.55,1.56,1.57,1.58,2]))
#budgets = [0, 0.1, 0.5, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.5, 2.9, 3.0]
budgets = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.76, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0]

hyperplanes3 = []
for i in range(len(budgets)):
    hyperplanes3.append(output3[output3.sum_prob == budgets[i]])

# =============================================================================
# "Plotting the influence hyperplanes when x1 + x2 + x3 = k"
# for i in range(len(budgets)-2):
#     %matplotlib
#     fig = plt.figure(i+1)
#     ax = Axes3D(fig)
#     surface = ax.plot_trisurf(hyperplanes3[i+1].prob1, hyperplanes3[i+1].prob2, hyperplanes3[i+1].influence, cmap=cm.jet, linewidth=0.1)
#     ax.set_title('Influence Hyperplane when x1 + x2 + x3 = '+str(budgets[i+1]))  
#     fig.colorbar(surface, shrink=0.5, aspect=5)
#     fig.savefig('fractional_incentives/results/optimal_solution/influence_hyperplane3-'+str(budgets[i+1])+'.png')
#     plt.pause(0.00001)
# =============================================================================

"influence along the line from best to best pair."
x1 = [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0]
x2 = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
x3 = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]

output3_line = pd.DataFrame()
for i in range(len(x1)):
    output3_line = pd.concat([output3_line, output3.loc[(output3.prob1 == x1[i]) & (output3.prob2 == x2[i]) & (output3.prob3 == x3[i])]]) 

"Optimal solutions"
optimal_solutions3 = []
for i in range(len(hyperplanes3)):
    optimal_solutions3.append(hyperplanes3[i][hyperplanes3[i]['influence'] == max(hyperplanes3[i]['influence'])])
optimal_solutions3 = pd.concat(optimal_solutions3)
print(optimal_solutions3)


"Corner to Corner"
corner_to_corner1 = pd.DataFrame()
corner_to_corner2 = pd.DataFrame()
corner_to_corner1 = pd.concat([corner_to_corner1,influence_corners_df[(influence_corners_df['prob1'] == 1) & (influence_corners_df['prob2'] == 0) & (influence_corners_df['prob3'] == 0)]])
corner_to_corner1 = pd.concat([corner_to_corner1,influence_corners_df[(influence_corners_df['prob1'] == 1) & (influence_corners_df['prob2'] == 0) & (influence_corners_df['prob3'] == 1)]])
corner_to_corner2 = pd.concat([corner_to_corner2,influence_corners_df[(influence_corners_df['prob1'] == 0) & (influence_corners_df['prob2'] == 0) & (influence_corners_df['prob3'] == 1)]])
corner_to_corner2 = pd.concat([corner_to_corner2,influence_corners_df[(influence_corners_df['prob1'] == 0) & (influence_corners_df['prob2'] == 1) & (influence_corners_df['prob3'] == 1)]])
  

"Plotting optimal influence as a function of the budget"   
plt.rcParams['figure.figsize'] = [5.2,5.2] 
plt.rcParams.update({'font.size': 11})
width = 3
#%matplotlib
fig1 = plt.figure()
plt.plot(optimal_solutions3.sum_prob,optimal_solutions3.influence,linewidth=width, label='$F(\mathbf{x}^*)$')
#plt.plot(output3_line.sum_prob,output3_line.influence)
plt.plot(corner_to_corner1.sum_prob,corner_to_corner1.influence, linewidth=width, label='line segment connecting $F(1,0,0)$ to $F(1,0,1)$')
plt.plot(corner_to_corner2.sum_prob,corner_to_corner2.influence, linewidth=width, label='line segment connecting $F(0,0,1)$ to $F(0,1,1)$')
plt.text(0.07,0.05,'$F(0,0,0)$')
plt.text(0.48,5.10,'$F(1,0,0)$')
plt.text(1.42,7.15,'$F(0,1,1)$')
plt.text(0.90,3.30,'$F(0,0,1)$')
plt.text(2.02,6.70,'$F(1,0,1)$')
plt.text(2.62,7.48,'$F(1,1,1)$')
plt.text(1.01,6.35,'$F(1,0,.76)$')
plt.text(1.80,6.00,'$F(0,.76,1)$')
#plt.legend()
plt.xlabel('$k$')
plt.ylabel('$F(\mathbf{x}^*)$')
#plt.xlabel('budget')
#plt.ylabel('expected optimal influence')
plt.grid()
fig1.savefig('im_outputs/solution/budget_3/opt_influence_vs_budget.png',dpi=300)

"Plotting the optimal solutions"
plt.rcParams['figure.figsize'] = [3,3]
plt.rcParams.update({'font.size': 11})
width = 3
#%matplotlib
fig = plt.figure()
ax = Axes3D(fig)
#optimal_solutions3 = optimal_solutions3[(optimal_solutions3.sum_prob>=1) & (optimal_solutions3.sum_prob<=2)]
scatter = ax.scatter(optimal_solutions3.prob1, optimal_solutions3.prob2, optimal_solutions3.prob3, c='b', marker='o')
#ax.set_title('Optimal Solutions for different $k$')
ax.set_xlabel('$x_1^*$')
ax.set_ylabel('$x_2^*$')
ax.set_zlabel('$x_3^*$')
ax.text(1.04,0,.77,'(1,0,.76)', color='b')
ax.text(-0.52,.78,.8,'(0,.76,1)', color='b')
#ax.set_title('Optimal Solutions for different Budgets')
#ax.set_xlabel('discount to person 1 (x1)')
#ax.set_ylabel('discount to person 2 (x2)')
#ax.set_zlabel('discount to person 3 (x3)')
fig.savefig('im_outputs/solution/budget_3/optimal_solutions3.png',dpi=300, bbox_inches = "tight")

######################
#### DEFENSE #########
######################
"Plotting optimal influence as a function of the budget -- defense presentation"   
plt.rcParams['figure.figsize'] = [5.2,5.2] 
plt.rcParams.update({'font.size': 11})
width = 3
#%matplotlib
fig1 = plt.figure()
plt.plot(optimal_solutions3.sum_prob,optimal_solutions3.influence,linewidth=width, label='$F(\mathbf{x}^*)$')
#plt.plot(output3_line.sum_prob,output3_line.influence)
#plt.plot(corner_to_corner1.sum_prob,corner_to_corner1.influence, linewidth=width, label='line segment connecting $F(1,0,0)$ to $F(1,0,1)$')
#plt.plot(corner_to_corner2.sum_prob,corner_to_corner2.influence, linewidth=width, label='line segment connecting $F(0,0,1)$ to $F(0,1,1)$')
plt.text(0.07,0.05,'$F(0,0,0)$')
plt.text(0.48,5.10,'$F(1,0,0)$')
plt.text(1.42,7.15,'$F(0,1,1)$')
#plt.text(0.90,3.30,'$F(0,0,1)$')
#plt.text(2.02,6.70,'$F(1,0,1)$')
plt.text(2.62,7.48,'$F(1,1,1)$')
plt.text(1.01,6.35,'$F(1,0,.76)$')
plt.text(1.80,6.00,'$F(0,.76,1)$')
#plt.legend()
plt.xlabel('$k$')
plt.ylabel('$F(\mathbf{d}^*)$')
#plt.xlabel('budget')
#plt.ylabel('expected optimal influence')
plt.grid()
fig1.savefig('im_outputs/solution/budget_3/opt_influence_vs_budget-defense1.png',dpi=300)

"Plotting optimal influence as a function of the budget -- defense presentation"   
plt.rcParams['figure.figsize'] = [5.2,5.2] 
plt.rcParams.update({'font.size': 11})
width = 3
#%matplotlib
fig1 = plt.figure()
plt.plot(optimal_solutions3.sum_prob,optimal_solutions3.influence,linewidth=width, label='$F(\mathbf{x}^*)$')
#plt.plot(output3_line.sum_prob,output3_line.influence)
plt.plot(corner_to_corner1.sum_prob,corner_to_corner1.influence, linewidth=width, label='line segment connecting $F(1,0,0)$ to $F(1,0,1)$')
plt.plot(corner_to_corner2.sum_prob,corner_to_corner2.influence, linewidth=width, label='line segment connecting $F(0,0,1)$ to $F(0,1,1)$')
plt.text(0.07,0.05,'$F(0,0,0)$')
plt.text(0.48,5.10,'$F(1,0,0)$')
plt.text(1.42,7.15,'$F(0,1,1)$')
plt.text(0.90,3.30,'$F(0,0,1)$')
plt.text(2.02,6.70,'$F(1,0,1)$')
plt.text(2.62,7.48,'$F(1,1,1)$')
plt.text(1.01,6.35,'$F(1,0,.76)$')
plt.text(1.80,6.00,'$F(0,.76,1)$')
#plt.legend()
plt.xlabel('$k$')
plt.ylabel('$F(\mathbf{d}^*)$')
#plt.xlabel('budget')
#plt.ylabel('expected optimal influence')
plt.grid()
fig1.savefig('im_outputs/solution/budget_3/opt_influence_vs_budget-defense2.png',dpi=300)


"Plotting greedy influence as a function of the budget -- defense presentation"   

"Corner to Corner"
corner_to_corner_greedy1 = pd.DataFrame()
corner_to_corner_greedy2 = pd.DataFrame()
corner_to_corner_greedy3 = pd.DataFrame()
corner_to_corner_greedy1 = pd.concat([corner_to_corner_greedy1,influence_corners_df[(influence_corners_df['prob1'] == 0) & (influence_corners_df['prob2'] == 0) & (influence_corners_df['prob3'] == 0)]])
corner_to_corner_greedy1 = pd.concat([corner_to_corner_greedy1,influence_corners_df[(influence_corners_df['prob1'] == 1) & (influence_corners_df['prob2'] == 0) & (influence_corners_df['prob3'] == 0)]])
corner_to_corner_greedy2 = pd.concat([corner_to_corner_greedy2,influence_corners_df[(influence_corners_df['prob1'] == 1) & (influence_corners_df['prob2'] == 0) & (influence_corners_df['prob3'] == 0)]])
corner_to_corner_greedy2 = pd.concat([corner_to_corner_greedy2,influence_corners_df[(influence_corners_df['prob1'] == 1) & (influence_corners_df['prob2'] == 0) & (influence_corners_df['prob3'] == 1)]])
corner_to_corner_greedy3 = pd.concat([corner_to_corner_greedy3,influence_corners_df[(influence_corners_df['prob1'] == 1) & (influence_corners_df['prob2'] == 0) & (influence_corners_df['prob3'] == 1)]])
corner_to_corner_greedy3 = pd.concat([corner_to_corner_greedy3,influence_corners_df[(influence_corners_df['prob1'] == 1) & (influence_corners_df['prob2'] == 1) & (influence_corners_df['prob3'] == 1)]])

plt.rcParams['figure.figsize'] = [5.2,5.2] 
plt.rcParams.update({'font.size': 11})
width = 3
#%matplotlib
fig1 = plt.figure()
plt.plot(corner_to_corner_greedy1.sum_prob,corner_to_corner_greedy1.influence, linewidth=width, label='line segment connecting $F(1,0,0)$ to $F(1,0,1)$',color='purple')
plt.plot(corner_to_corner_greedy2.sum_prob,corner_to_corner_greedy2.influence, linewidth=width, label='line segment connecting $F(1,0,0)$ to $F(1,0,1)$',color='purple')
plt.plot(corner_to_corner_greedy3.sum_prob,corner_to_corner_greedy3.influence, linewidth=width, label='line segment connecting $F(1,0,0)$ to $F(1,0,1)$',color='purple')
plt.text(0.07,0.05,'$F(0,0,0)$')
plt.text(0.48,5.10,'$F(1,0,0)$')
plt.text(1.40,6.70,'$F(1,0,1)$')
plt.text(2.40,8.15,'$F(1,1,1)$')
#plt.legend()
plt.xlabel('$k$')
plt.ylabel('$F(\mathbf{d}^*)$')
#plt.xlabel('budget')
#plt.ylabel('expected optimal influence')
plt.grid()
fig1.savefig('im_outputs/solution/budget_3/greedy_influence_vs_budget-defense.png',dpi=300)


"Plotting the optimal solutions -- defense presentation"
plt.rcParams['figure.figsize'] = [3,3]
plt.rcParams.update({'font.size': 11})
width = 3
#%matplotlib
fig = plt.figure()
ax = Axes3D(fig)
#optimal_solutions3 = optimal_solutions3[(optimal_solutions3.sum_prob>=1) & (optimal_solutions3.sum_prob<=2)]
scatter = ax.scatter(optimal_solutions3.prob1, optimal_solutions3.prob2, optimal_solutions3.prob3, c='b', marker='o')
#ax.set_title('Optimal Solutions for different $k$')
ax.set_xlabel('Discount to user 1 ($d_1^*$)')
ax.set_ylabel('Discount to user 2 ($d_2^*$)')
ax.set_zlabel('Discount to user 3 ($d_3^*$)')
ax.text(1.04,0,.77,'(1,0,.76)', color='b')
ax.text(-0.52,.78,.8,'(0,.76,1)', color='b')
#ax.set_title('Optimal Solutions for different Budgets')
#ax.set_xlabel('discount to person 1 (x1)')
#ax.set_ylabel('discount to person 2 (x2)')
#ax.set_zlabel('discount to person 3 (x3)')
fig.savefig('im_outputs/solution/budget_3/optimal_solutions3-defense.png',dpi=300, bbox_inches = "tight")

end_time=time.time()
print("run time: " + str((end_time-start_time)) + " sec")
