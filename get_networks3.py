#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 31 15:50:49 2025

@author: srogoz@ice.mpg.de
"""
import pandas as pd
import numpy as np
import os

os.chdir('/home/srogoz@ice.mpg.de/ulr-lab/Users/Sarah_saro7615/EXP1_GROUPSIZEandCOMPOSITION/')

###check directory
#os.path.isdir('EXP1_NetworkAnalysis/EXP1_colonies/multindex')

#INPUT: multiindex ONLY over ant and frame
# OUTPUT: temporal interaction matrix with dim: present ants x present ants x number of total frames
# if ant is NA but present not registered as having any interactions
###################
# %% load xy-data
# multindex
colony_dir = 'processed/EXP1_colonies/a16-4_interpo_20250822.parquet'
#os.path.isdir('processed/EXP1_colonies/a16-4_interpo_20250822.parquet')

colony = pd.read_parquet(colony_dir)

###
# %% format colony: drop colony_name and multiindex ant and frame 

colony = (colony.drop(columns =["colony"]).set_index(["frame","ant"]).sort_index())

####################
#%% vectorized version for faster running
#step 1: prep array
def prepare_array (df):
    ants = df.index.get_level_values("ant").unique()
    frames = df.index.get_level_values("frame").unique()
    
    #into integers
    ant_to_idx = {ant: i for i, ant in enumerate(ants)}
    frame_to_idx ={frame: i for i, frame in enumerate(frames)}
    
    #structure of coordinate matrix
    coordinates = np.full((len(frames), len(ants),2), np.nan)
    
    #initialize by iterating through rows
    for (ant, frame), row in df[["x","y"]].iterrows():
        coordinates[frame_to_idx[frame],ant_to_idx[ant], :] = row.values
        
    return coordinates, ants, frames
#################
# %% distance measure

def pairwise_distance(points):
    diff = points[:, np.newaxis, :] - points[np.newaxis, :, :]  # shape (n_ants, n_ants, 2)
    dist = np.sqrt(np.sum(diff**2, axis=-1))                    # shape (n_ants, n_ants)
    return dist
    

###################
# %%
def get_network_vec(df, interaction_threshold=0.002):
    coordinates, ants, frames = prepare_array(df)
    n_frames, n_ants, _= coordinates.shape
    
    interaction_matrix = np.zeros((n_ants, n_ants, n_frames) , dtype = np.uint8 )
    
    for i in range(n_frames):
        points = coordinates[i]
        
        if np.isnan(points).all():
            continue #skip empty frames
        dist = pairwise_distance(points)
        mask = (dist <= interaction_threshold) & ~np.eye(n_ants, dtype=bool)
        
        #keep only upper tri i<j
        mask = np.triu(mask, k=1)
        
        interaction_matrix[:, :, i] = mask.astype(np.uint8)
    
    return interaction_matrix, ants, frames
        
#####################
#%%run 
interaction_matrix, ants, frames = get_network_vec(colony, interaction_threshold = 0.002)
###############
#%%
interaction_matrix[::60000]
colony.loc[("BP",525), "x"]
colony.index.get_level_values("frame").max()
############
#%% assignmentrate per frame
na_per_frame = colony[['x','y']].isna().any(axis=1)  # True if either x or y is NaN
na_count = na_per_frame.groupby(level='frame').sum() 
mean_assignmentrate = np.mean(na_count)
sd_assignment = np.sqrt(np.var(na_count))
###############
#%%
import matplotlib.pyplot as plt

plt.figure(figsize=(14,6))
plt.plot(na_count.index, na_count.values, marker='o')
plt.title("Assignmentrate per fr ")
plt.xlabel("Frame")
plt.ylabel("Number of ants with NaN (x or y)")
plt.grid(True)
plt.show()
######################
#%%
collapsed = interaction_matrix.sum(axis=2)

###############
#%%

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

filename = f"a16-4_proxMatrInterpo_{today}.npy"

filepath = f"processed/networks/{filename}"

np.save(filepath, interaction_matrix)

###############
# %% SAVE PRESENT ANTS

from datetime import datetime

today= datetime.today().strftime('%Y%m%d')

filename = f"b16-1_present_ants_{today}.npy"

filepath = f"processed/present_ants/{filename}"

np.save(filepath, present_ants)
