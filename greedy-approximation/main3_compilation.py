#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 21:49:19 2019

@author: abhishek.umrawal
"""
import pickle as pickle
import os as os
import pandas as pd
import matplotlib.pyplot as plt
#matplotlib
import tikzplotlib
plt.rcParams['figure.figsize'] = [8,6]
plt.rcParams.update({'font.size': 12})

"list of algorithms"
all_algorithms = ['celfpp','ris','lim']

"inputs"
chosen_algorithm_indices = [0,1,2]
budget =  20 
interval = 1                      
name_id = '_facebook'
probs = [0.20, 0.40, 0.60, 0.80]

"looking up for nested non-adaptive methods"
chosen_algorithms = [item for i,item in enumerate(all_algorithms) if i in chosen_algorithm_indices]
results_folder = 'results'+os.sep+'results'+name_id
output_dict={}
results_dicts = {}
for algorithm in chosen_algorithms:  
    filename = 'output_'+algorithm+'__'+str(budget)+'__frac__.pkl'
    filename_with_path = results_folder+os.sep+'pickle_files'+os.sep+filename
    if os.path.exists(filename_with_path):
        with open(filename_with_path, 'rb') as f:
            results_dict = pickle.load(f) 
        results_dicts[algorithm] = results_dict
        #best_seed_sets = results_dict['best_seed_set']  
        exp_influences = results_dict['exp_influence']
        
        if len(exp_influences) > 1:              
            #output_dict[algorithm+'_best_seed_sets'] = best_seed_sets  
            output_dict[algorithm+'_exp_influences'] = exp_influences

if 0 in chosen_algorithm_indices:
    "celf++ floor calculations"
    celfpp_exp_influences = output_dict['celfpp_exp_influences']
    celfpp_floor_values = [celfpp_exp_influences[i] for i in range(budget*(len(probs)+1)) if (i+1)%(len(probs)+1) == 1]
    celfpp_floor_exp_influences = []
    for val in celfpp_floor_values:
        celfpp_floor_exp_influences+=[val for i in range(len(probs)+1)]
    celfpp_floor_exp_influences+=[celfpp_exp_influences[-1]]
    output_dict['celfpp_floor_exp_influences'] = celfpp_floor_exp_influences


if 1 in chosen_algorithm_indices:
    "ris floor calculations"
    ris_exp_influences = output_dict['ris_exp_influences']
    ris_floor_values = [ris_exp_influences[i] for i in range(budget*(len(probs)+1)) if (i+1)%(len(probs)+1) == 1]
    ris_floor_exp_influences = []
    for val in ris_floor_values:
        ris_floor_exp_influences+=[val for i in range(len(probs)+1)]
    ris_floor_exp_influences+=[celfpp_exp_influences[-1]]
    output_dict['ris_floor_exp_influences'] = ris_floor_exp_influences

"saving output as a dataframe"
if not os.path.exists(results_folder+os.sep+'plots'):
    os.makedirs(results_folder+os.sep+'plots')
columns = list(output_dict.keys())
output_df = pd.DataFrame()
for column in columns:
    output_df[column] = output_dict[column]
    
"creating a list of all budgets"
budgets = [0]
for int_budget in range(budget):
    for frac_budget in probs+[1]:
        budgets.append(int_budget+frac_budget)        
output_df.index = budgets

"saving output_df as a csv file"
output_df.to_csv(results_folder+os.sep+'plots'+'/sim_exp_influence_data_'+str(budget)+'.csv',index=False,header=True)

"plotting"
chosen = [colname for colname in output_df.columns if 'exp_' in colname]             
ax = output_df.loc[:,chosen].plot(lw=2,marker='o')
ax.set_xlabel('$k$')
ax.set_ylabel('Expected Influence')
legend_entries = ['GreedyMultiLinExt-CELF++','GreedyMultiLinExt-RIS','LIM','FloorMultiLinExt']
chosen_legend_entries = [item for i,item in enumerate(legend_entries) if i in chosen_algorithm_indices]
ax.legend(chosen_legend_entries)
ax.grid(linestyle='-', linewidth=1)
ax.set_title('Influence maximization for '+name_id[1:]+' network')
#ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))  
ax.get_figure().savefig(results_folder+os.sep+'plots'+'/sim'+name_id+'_plot_'+str(budget)+'.eps')
ax.get_figure().savefig(results_folder+os.sep+'plots'+'/sim'+name_id+'_plot_'+str(budget)+'.jpg')
tikzplotlib.clean_figure()
tikzplotlib.save(results_folder+os.sep+'plots'+'/sim'+name_id+'_plot_'+str(budget)+'.tex')