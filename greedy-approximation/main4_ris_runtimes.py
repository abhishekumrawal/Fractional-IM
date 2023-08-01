#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 23:10:19 2022

@author: abhishek.umrawal
"""

import numpy as np
import pandas as pd
import os
import tikzplotlib
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [8,6]
plt.rcParams.update({'font.size': 12})


"input"                  
name_id = '_deezer'
budget = 20
probs = [0.20, 0.40, 0.60, 0.80]

"ris code folder"
ris_folder = './ris_code_release'

run_times_ris=[0]
for budget in range(1,budget+1):
    "output filename for run time"
    out_filename = ris_folder+"/result/"+name_id[1:]+"_subsim_k"+str(budget)+"_wc"
    
    "reading run times from file and appending a list"
    run_times_ris.append(float(open(out_filename).readlines()[1].split(' ')[2]))

run_times_ris_with_frac = []
for val in run_times_ris[:-1]:
    run_times_ris_with_frac+=[val for i in range(len(probs)+1)]
run_times_ris_with_frac+=[run_times_ris[-1]]

"saving run times for ris a text file"
results_folder_runtime_files = './results'+os.sep+'results'+name_id+os.sep+'runtime_files'
fstr = results_folder_runtime_files+os.sep+'runtime_info_ris_all.txt'
with open(fstr, 'w') as f:
    for val in run_times_ris_with_frac:
        f.write(str(val))
        f.write('\n')

"lim code folder"
lim_folder = './lim_code_release'

"output filename for run time"
out_filename_time = lim_folder+"/time/"+name_id[1:]+".txt_cimm_eps=0.500000_group_0_new"

"#### Getting the runtimes (cumulative in seconds)"
run_times_lim = [x.split(' ')[7] for x in open(out_filename_time).readlines()]
run_times_lim = [0] + [float(x) for x in run_times_lim]
run_times_lim = np.cumsum(run_times_lim)

"creating a list of all budgets"
budgets = [0]
for int_budget in range(budget):
    for frac_budget in probs+[1]:
        budgets.append(int_budget+frac_budget)  

"plotting"
results_folder = 'results'+os.sep+'results'+name_id
if not os.path.exists(results_folder+os.sep+'plots'):
    os.makedirs(results_folder+os.sep+'plots')

"actual run times"
output_df = pd.DataFrame()
output_df['ris'] = run_times_ris_with_frac
output_df['lim'] = run_times_lim
output_df.index = budgets           
ax = output_df.plot(lw=2,marker='o')
ax.set_xlabel('Budget')
ax.set_ylabel('Run-time')
ax.grid(linestyle='-', linewidth=1)
ax.set_title('Influence maximization for '+name_id[1:]+' network')
legend_entries = ['MLE-Greedy(RIS)','LIM-Greedy(RIS)']
ax.legend(legend_entries)
tikzplotlib.clean_figure()
tikzplotlib.save(results_folder+os.sep+'plots'+'/sim'+name_id+'_plot_'+str(budget)+'_runtimes.tex')

"log run times"
output_df = pd.DataFrame()
output_df['ris'] = np.log(run_times_ris_with_frac)
output_df['lim'] = np.log(run_times_lim)
output_df.index = budgets           
ax = output_df.plot(lw=2,marker='o')
ax.set_xlabel('Budget')
ax.set_ylabel('ln(Run-time)')
ax.grid(linestyle='-', linewidth=1)
ax.set_title('Influence maximization for '+name_id[1:]+' network')
legend_entries = ['MLE-Greedy(RIS)','LIM-Greedy(RIS)']
ax.legend(legend_entries)
tikzplotlib.clean_figure()
tikzplotlib.save(results_folder+os.sep+'plots'+'/sim'+name_id+'_plot_'+str(budget)+'_ln_runtimes.tex')