#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 22:33:27 2022

@author: abhishek.umrawal
"""

import os

for i in range(5):
    print('I AM AT ITERATION NO. '+str(i)+'.')
    os.system('python main1_offline.py')
    os.system('python main4_runtimes.py')

#os.system('python main5_avg_runtimes.py')