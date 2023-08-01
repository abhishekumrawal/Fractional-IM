#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 21:24:44 2022

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
name_id = '_wikipedia'
budget = 20
probs = [0.20, 0.40, 0.60, 0.80]
num_runs = 10

"run time files folder"
results_folder_runtime_files = './results'+os.sep+'results'+name_id+os.sep+'runtime_files'

"reading relevant file names"
os.chdir(results_folder_runtime_files)
filenames = []
for filename in os.listdir():
    if 'runtime_info_ris_all' in filename:
        filenames.append(filename)
os.chdir('..')
os.chdir("..")
os.chdir("..")

"reading run times for ris from a text file"
run_times_ris_all_runs = []
for filename in filenames:
    out_filename_time = results_folder_runtime_files+os.sep+filename
    run_times_ris_all_runs.append([float(x) for x in open(out_filename_time).readlines()])

run_times_ris = [float(sum(col))/len(col) for col in zip(*run_times_ris_all_runs)]

"reading run times for lim from a text file"
out_filename_time = results_folder_runtime_files+os.sep+'runtime_info_lim_all.txt'
run_times_lim = [float(x) for x in open(out_filename_time).readlines()]
run_times_lim = [float(0)] + [float(x) for x in run_times_lim]
#run_times_lim = [float(x) for x in run_times_lim]
run_times_lim = np.cumsum(run_times_lim)

"creating a list of all budgets"
budgets = [float(0)]
for int_budget in range(budget):
    for frac_budget in probs+[1]:
        budgets.append(int_budget+frac_budget)  

"plotting"
results_folder = 'results'+os.sep+'results'+name_id
if not os.path.exists(results_folder+os.sep+'plots'):
    os.makedirs(results_folder+os.sep+'plots')

"actual run times"
output_df = pd.DataFrame()
output_df['ris'] = run_times_ris
output_df['lim'] = run_times_lim
output_df.index = budgets           
ax = output_df.plot(lw=0,marker='o')
ax.set_xlabel('Budget')
ax.set_ylabel('Run-time')
ax.grid(linestyle='-', linewidth=1)
ax.set_title('Influence maximization for '+name_id[1:]+' network')
legend_entries = ['MLE-Greedy','LIM-Greedy']
ax.legend(legend_entries)
tikzplotlib.clean_figure()
tikzplotlib.save(results_folder+os.sep+'plots'+'/sim'+name_id+'_plot_'+str(budget)+'_runtimes_avg.tex')

"log scale plot"           
ax = output_df.plot(lw=0,marker='o',logy=True)
ax.set_xlabel('Budget')
ax.set_ylabel('Run-time')
ax.grid(linestyle='-', linewidth=1)
ax.set_title('Influence maximization for '+name_id[1:]+' network')
legend_entries = ['MLE-Greedy','LIM-Greedy']
ax.legend(legend_entries)
tikzplotlib.clean_figure()
tikzplotlib.save(results_folder+os.sep+'plots'+'/sim'+name_id+'_plot_'+str(budget)+'_ln_runtimes_avg.tex')