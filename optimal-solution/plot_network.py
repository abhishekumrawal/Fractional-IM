#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 20:22:51 2019

@author: abhishek.umrawal
"""

"Importing necessary modules"

import networkx as nx
import matplotlib.pyplot as plt
%matplotlib


"Declaring the Florentine families network"
network = nx.florentine_families_graph()
network = network.to_directed()
#network = nx.convert_node_labels_to_integers(network,first_label=1)

options = {
    #'node_color': 'blue',
    'node_size': 1000,
    #'width': 4,
    'arrowstyle': '-|>',
    'arrowsize': 30,
}

pos=nx.spring_layout(network)
nx.draw(network,pos, font_size = 15, font_family='times', node_color='w',edge_color='#BB0000',width=1,with_labels=True, **options)
nx.draw_networkx_nodes(network,pos, node_color='#A0CBE2', with_labels=True, node_size=1000, node_shape='o', alpha=None,)        
nx.draw_networkx_edges(network,pos, width = 3, **options)

        
plt.savefig("small_network.png", dpi=500, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None,transparent=False, bbox_inches=None, pad_inches=0.1)



#pos=nx.spring_layout(network)
#nx.draw(network,pos,node_color='#A0CBE2',edge_color='#BB0000',width=2,edge_cmap=plt.cm.Blues,with_labels=True)
#plt.savefig("small_network.png", dpi=500, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None,transparent=False, bbox_inches=None, pad_inches=0.1)
