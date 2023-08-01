#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 21:30:04 2019

@author: abhishek.umrawal
"""

"importing required built-in modules"
import os
import timeit


    
"#### Changing directory to the Goyal's software folder, compiling the C codes by doing make and then coming back"
os.chdir(celfpp_folder_new)
os.system("make")
os.chdir("..")
os.chdir("..")

"##### Calling CELF++"    
start = timeit.default_timer()
os.system(celfpp_folder_new+"/InfluenceModels -c "+celfpp_folder_new+"/config_test.txt")
end = timeit.default_timer()
runtime = end - start
