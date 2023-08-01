#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 21:01:00 2019

@author: abhishek.umrawal
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 19:17:14 2019

@author: abhishek.umrawal
"""

"importing required user-defined modules"
from im_functions.weighted_network import weighted_network
from im_functions.celfpp_im import celfpp_im
from im_functions.ris_im import ris_im
from im_functions.lim_im import lim_im
from im_functions.ccelfpp1_im import ccelfpp1_im
from im_functions.heuristic_im import heuristic_im

def non_adaptive_im(inpt):
    
    network, weighting_scheme, algorithm, heuristic, max_budget, diffusion_model, n_sim, name_id, community_method, communities, community_size_threshold, is_graph_already_weighted = inpt
    
    if is_graph_already_weighted == False:
        "## adding weights to the unweighted network"
        network = weighted_network(network,method = weighting_scheme)
    
    if (algorithm == "celfpp"):
        best_seed_sets, exp_influences, runtime = celfpp_im(network, max_budget, diffusion_model, n_sim, all_upto_budget=True)
    elif (algorithm == "ris"):
        best_seed_sets, exp_influences, runtime = ris_im(network, max_budget, diffusion_model, n_sim, all_upto_budget=True)
    elif (algorithm == "lim"):
        best_seed_sets, exp_influences, runtime = lim_im(network, max_budget, diffusion_model, n_sim, all_upto_budget=True)
    elif (algorithm == "ccelfpp1"):
        best_seed_sets, exp_influences,runtime = ccelfpp1_im(network, max_budget, diffusion_model, n_sim, community_method, communities, community_size_threshold, all_upto_budget=True)

    elif (algorithm == "heuristic"):
        best_seed_sets, exp_influences,runtime = heuristic_im(network, heuristic, max_budget, diffusion_model, n_sim, all_upto_budget=True)
    return 