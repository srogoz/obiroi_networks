#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 31 15:50:49 2025

@author: srogoz@ice.mpg.de
"""
import pandas as pd
import numpy as np
import os
#os.path.isdir('EXP1_NetworkAnalysis/EXP1_colonies/multindex')

#INPUT: multiindex ONLY over ant and frame
# OUTPUT: temporal interaction matrix with dim: present ants x present ants x number of total frames
# if ant is NA but present not registered as having any interactions
###################
# %% load xy-data
# multindex
colony_dir = 'EXP1_NetworkAnalysis/EXP1_colonies/multindex/b16-1_interpo_20250731.pkl'
#colony_dir = 'EXP1_NetworkAnalysis/EXP1_colonies/singleindex/b16-1_interpo_20250731.pkl'
colony = pd.read_pickle(colony_dir)

###
# %% format colonies 
colony = colony.reset_index('colony', drop = True)

#############
# %%


import pandas as pd
import numpy as np

interaction_threshold = 0.002  # bodylength ants

def get_network(df, interaction_threshold):
    total_frames = df.index.get_level_values("frame").max() + 1
    present_ants = df.index.get_level_values('ant').unique()
    number_ants = len(present_ants)
    
    # Create mapping from ant name to matrix index
    ant_to_idx = {ant: idx for idx, ant in enumerate(present_ants)}
    
    # Initialize 3D interaction matrix
    interaction_matrix = np.zeros((number_ants, number_ants, total_frames))

    for i in range(total_frames):
        for ant1 in present_ants:
            for ant2 in present_ants:
                if ant1 == ant2:
                    continue

                try:
                    x1 = df.loc[(ant1, i), "x"]
                    y1 = df.loc[(ant1, i), "y"]
                    x2 = df.loc[(ant2, i), "x"]
                    y2 = df.loc[(ant2, i), "y"]
                except KeyError:
                    continue  # missing data

                if pd.isna(x1) or pd.isna(y1) or pd.isna(x2) or pd.isna(y2):
                    continue

                distance = np.sqrt((x1 - x2)**2 + (y1 - y2)**2)

                if distance <= interaction_threshold:
                    idx1 = ant_to_idx[ant1]
                    idx2 = ant_to_idx[ant2]
                    interaction_matrix[idx1, idx2, i] = 1

    return interaction_matrix
                 
        
         
############
# %%
interaction_matrix = get_network(colony, interaction_threshold)

# %%
#interaction_matrix[:,:,3000]
#present_ants = colony.index.get_level_values('ant').unique()
#os.path.isdir('processed/networks')
############
# %% save networks

from datetime import datetime

today= datetime.today().strftime('%Y%m%d')

filename = f"b16-1_proxMatrInterpo_{today}.npy"

filepath = f"processed/networks/{filename}"

np.save(filepath, interaction_matrix)

###############
# %% SAVE PRESENT ANTS

from datetime import datetime

today= datetime.today().strftime('%Y%m%d')

filename = f"b16-1_present_ants_{today}.npy"

filepath = f"processed/present_ants/{filename}"

np.save(filepath, present_ants)
