#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 17 13:34:39 2025

@author: srogoz@ice.mpg.de
"""

##############
##%%finding colonies
import h5py
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import datetime
import pathlib

#vectorized network extraction
from get_vecnetworks import prepare_array, pairwise_distance, get_network_vec
#find continous nan runs function
from find_continousNANruns import find_continuous_NAN_runs
#reassign nan runs function
from reassign_nan_runs import reassign_nan_runs
# load all files for one colony import h5py
from merge_files_per_colony import merge_files_for_colony
##########
#%%
#initialisting root base directory For EXP1
os.chdir('/home/srogoz@ice.mpg.de/ulr-lab/Users/Sarah_saro7615/EXP1_GROUPSIZEandCOMPOSITION/')
#saving raw and interpolated xy-data to:
folder = "processed/EXP1_colonies_inf"
#saving networks to:
folder_networks = "processed/networks_inf"
folder_presentants = "processed/present_ants_inf"


##############
#%% parameters and DATA
#infected data
colony_directories = { #cam1
                        "b16-1" : 'sarah_3216ab_1_inf_20240626_161641/b16-1/antdata',
                      "ba16-1" : 'sarah_3216ab_1_inf_20240626_161641/ba16-1/antdata',
                      "a16-1" : 'sarah_3216ab_1_inf_20240626_161641/a16-1/antdata',
                      #cam2
                      "b16-2" : 'sarah_3216ab_2_inf_20240626_161710/1/antdata/C1',
                      "ba16-2" : 'sarah_3216ab_2_inf_20240626_161710/1/antdata/C2',
                      "a16-2" : 'sarah_3216ab_2_inf_20240626_161710/2/antdata/C1',
                      "b16-3" : 'sarah_3216ab_2_inf_20240626_161710/2/antdata/C2',
                      "ba16-3" : 'sarah_3216ab_2_inf_20240626_161710/3/antdata/C1',
                      "a16-3" : 'sarah_3216ab_2_inf_20240626_161710/3/antdata/C2',
                      #cam3 
                      "b16-4" : 'sarah_3216ab_3_inf_20240626_161759/b16-4/antdata',
                      "ba16-4" : 'sarah_3216ab_3_inf_20240626_161759/ba16-5/antdata',
                      "a16-4" : 'sarah_3216ab_3_inf_20240626_161759/a16-5/antdata',
                      #cam4 
                      "b16-5" : 'sarah_3216ab_4_inf_20240626_161846/b16-5/antdata',
                      "ba16-5" : 'sarah_3216ab_4_inf_20240626_161846/ba16-5/antdata',
                      "a16-5" : 'sarah_3216ab_4_inf_20240626_161846/a16-5/antdata',
                      #cam5 
                     "b16-6" : 'sarah_3216ab_5_inf_20240626_161922/b16-6/antdata',
                     "ba16-6" : 'sarah_3216ab_5_inf_20240626_161922/ba16-6/antdata',
                     "a16-6" : 'sarah_3216ab_5_inf_20240626_161922/a16-6/antdata',
                     #cam6 
                     "b16-7" : 'sarah_3216ab_6_inf_20240626_161959/b16-7/antdata',
                     "ba16-7" : 'sarah_3216ab_6_inf_20240626_161959/ba16-7/antdata',
                     "a16-7" : 'sarah_3216ab_6_inf_20240626_161959/a16-7/antdata',
                     "b16-8" : 'sarah_3216ab_6_inf_20240626_161959/b16-8/antdata',
                     "ba16-8" : 'sarah_3216ab_6_inf_20240626_161959/ba16-8/antdata',
                     "a16-8" : 'sarah_3216ab_6_inf_20240626_161959/a16-8/antdata'
                                            
                      }


expected_ants = ["BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", 
                   "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"]

selected_colonies_16ants_10frames = ["a16-1", "a16-2", "a16-3", "a16-4",  "a16-6", "b16-1", "b16-3", "b16-4", "b16-5",  "b16-7", "ba16-1", "ba16-4", "ba16-5"]
selected_colonies_16ants_5frames = [ "a16-7", "a16-8","b16-7", "b16-8","ba16-8"]
selected_colonies_notfinished = ["b16-2", "ba16-2", "ba16-3","a16-5","b16-6", "ba16-6","ba16-7"]
colony_names = ["a16-1", "a16-2", "a16-3", "a16-4", "a16-5", "a16-6", "a16-7", "a16-8", "b16-1", "b16-2", "b16-3", "b16-4", "b16-5", "b16-6", "b16-7", "b16-8",  "ba16-1", "ba16-2", "ba16-3", "ba16-4", "ba16-5", "ba16-6", "ba16-7", "ba16-8"]
#distance for interpolation based on distance is one body length of an ant and 2 seconds
nanThreshold_d = 0.002
nanThreshold_t = 5 #change to 1s measured in frames

today= datetime.today().strftime('%Y%m%d')
##############
#%% find and sort all .mat files per colony
sorted_files_per_colony = {}

for colony_name in colony_names :
    ##list and sort all files in the directory, f.stem: file without extension
     ##split(x)by sign and look at first one
     
    colony_files = list(pathlib.Path(colony_directories[colony_name]).glob('*.mat'))
    colony_files_sorted = sorted(colony_files, key=lambda f: int(f.stem.split('_')[1]))
    sorted_files_per_colony [colony_name] = colony_files_sorted

##########
#%% loop 

for selected_colony in selected_colonies_16ants_5frames:
    selected_colonies = [selected_colony]


    combined_colony_data = merge_files_for_colony(selected_colony, colony_directories, expected_ants, colony_names)
    #assignmentrate per colony per ant
    
    assignmentrate_percolony_perant = pd.DataFrame(index=expected_ants, columns=colony_names)
    
    for (col, ant), df in combined_colony_data.groupby(level=["colony", "ant"]):
        if "x" not in df.columns:
            continue  # skip if 'x' column is missing
    
        x = df["x"]  # Get the x-column
        if len(x) == 0:
            continue  # avoid division by zero
    
        nan_percent = x.isna().sum() / len(x) * 100
        assignmentrate = 100 - nan_percent  # percentage of frames that are NOT NaN
    
        assignmentrate_percolony_perant.loc[ant, col] = assignmentrate
    
    
    
    filename = f"assignmentrates_{selected_colony}_{today}.npy"
    #filename = f"assignmentrates_{selected_colonies}_{today}.npy"
    filepath = f"processed/assignment_rates_before_interpo_inf/{filename}"
    
    np.save(filepath, assignmentrate_percolony_perant) 
    
    # apply run function
    all_nan_runs = []
    
    for (col, ant), df in combined_colony_data.groupby(level=["colony", "ant"]):
        result = find_continuous_NAN_runs(df, col, ant)
        all_nan_runs.append(result)
    
    # Concatenate just once
    result_continuous_nan_runs = pd.concat(all_nan_runs, ignore_index=True)
    #assign multiindex
    result_continuous_nan_runs.set_index(["colony", "ant"], inplace=True)   
    
    
    # interpolate na-runs
    combined_colony_data_interpolated, interpolation_summary, boundary_cases = reassign_nan_runs(
        combined_colony_data,
        result_continuous_nan_runs,
        selected_colonies,
        expected_ants,
        nanThreshold_d,
        nanThreshold_t
    )
    
    
    # adjust for multindexing leading to 1-16 being NA
    interpolation_summary_2 = interpolation_summary.iloc[0, 16:32]
    # save interpolation_summary
    filename = f"interpolation_summary_{selected_colony}_{today}.npy"
    filepath = f"processed/interpolation_summary_inf/{filename}"
    np.save(filepath, interpolation_summary_2) 
    
    
    # save as parquet to be used in R in folder for networks
    #reset multiindex, interpolated
    filename1 = f"{selected_colony}_interpo_{today}.parquet"
    file_path = os.path.join(folder, filename1)
    df_reset3 = combined_colony_data_interpolated.reset_index()
    df_reset3.to_parquet(file_path, index = False)
    
    #reset multiindex, ORIGINAL
    filename = f"{selected_colony}_{today}.parquet"
    file_path2 = os.path.join(folder, filename)
    df_reset2 = combined_colony_data.reset_index()
    df_reset2.to_parquet(file_path2, index = False)
    
    #run interaction_matrix in form of (ant, frame) multiindex dropped colony 
    
    interaction_matrix, ants, frames = get_network_vec(combined_colony_data_interpolated, interaction_threshold = 0.002)
    
    #save networks and present ants
    filename1 = f"{selected_colony}_proxMatrInterpo_{today}.npy"
    filepath = f"{folder_networks}/{filename1}"
    np.save(filepath, interaction_matrix)
    
    
    filename = f"{selected_colony}_present_ants_{today}.npy"
    filepath = f"{folder_presentants}/{filename}"
    np.save(filepath, ants)
