#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 17 12:06:01 2025

@author: srogoz@ice.mpg.de
"""
import numpy as np

####################
#%% vectorized version for faster running
#step 1: prep array
def prepare_array (df):
    #reset multindex
    df = df.reset_index()
    #multiindex and drop colony column 
    df = (df.drop(columns =["colony"]).set_index(["ant","frame"]).sort_index())

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
# %% step 2: definition distance measure

def pairwise_distance(points):
    diff = points[:, np.newaxis, :] - points[np.newaxis, :, :]  # shape (n_ants, n_ants, 2)
    dist = np.sqrt(np.sum(diff**2, axis=-1))                    # shape (n_ants, n_ants)
    return dist
    

###################today= datetime.today().strftime('%Y%m%d')
# %% runtime: max 5 min WORKING ONE!!!
#INPUT: multindex(ant, frame)
#OUTPUT: interactionmatrix:16x16xframes no col and row names
# ants : col and row names
#frames: number of frames

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
        