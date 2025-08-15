#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 31 15:50:49 2025

@author: srogoz@ice.mpg.de
"""

###################
# %% 

import pandas as pd
import numpy as np

colony_single = colony.reset_index()

def get_networks(df, interaction_threshold):
    total_frames = df.index.get_level_values("frame").max()
    present_ants = df.index.get_level_values('ant').unique()
    number_ants = len(present_ants)
    #initialize interaction matrix
    interaction_matrix = np.zeros(number_ants, number_ants, total_frames)
    
    
    for i in range(total_frames)
     for ant1 in  present_ants:
         for ant2 in present_ants:
             if ant1 == ant2 
             continue
             x1 = df.loc[(ant1, i), "x"]
             y1 = df.loc[(ant1, i), "y"]
             x2 = df.loc[(ant1, i), "x"]
             y2 = df.loc[(ant1, i), "y"]
             
             distance = 
         
############
# %%
colony.index.get_level_values("frame").max()

colony.loc[("BB",35), "x"]