#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 16 08:45:36 2025

@author: srogoz@ice.mpg.de
"""
###########
#%%
#function to deal with nans
## INPUT: combined_colony_data , result_nan_runs , threshold_d, threshold_t, selected_colonies, expected_ants
### OUTPUT: combined_colony_data_interpolated, interpolation_summary(percentage of frames interpolated of total frames per ant), boundary_cases(collects nan at end and nan in beginning cases)
###
#check for potential jumps in coordinates afterwards
import pandas as pd
import numpy as np


###########
#%%
def reassign_nan_runs(df_data, df_nan_runs, selected_colonies, expected_ants, threshold_d, threshold_t):
    combined_colony_data_interpolated = df_data.copy()
    combined_colony_data_interpolated = combined_colony_data_interpolated.sort_index()
    print(combined_colony_data_interpolated.head())
    #sum interpolated frames in percent
    interpolation_summary = pd.DataFrame(index = selected_colonies, columns = expected_ants)
    boundary_cases = pd.DataFrame(0, index = selected_colonies, columns = expected_ants)
    df_nan_runs = df_nan_runs.sort_index()
    total_runs = sum(len(df_nan_runs.loc[(col, ant)]) for col in selected_colonies for ant in expected_ants)
    
    for col in selected_colonies:
     
        for ant in expected_ants:
            #if not all ants are present
            if (col, ant) not in df_nan_runs.index:
                continue
            total_num_frames = len(df_data.loc[(col, ant, slice(None))]["x"])
            df_runs = df_nan_runs.loc[(col, ant)]
            df_runs = df_runs.sort_index()
            interpolation_sum_ant =0
            
            counter= 0
            for i in range(len(df_runs["run_length"])):
                 counter += 1
                 print(f"Processing run {counter}/{total_runs} (Colony: {col}, Ant: {ant}, Run index: {i+1})")
                #get parameters
                 
                 dist = df_runs["distance"].iloc[i]
                 #get parameters
                 x_before = df_runs["x_before_NA"].iloc[i]
                 x_after = df_runs["x_after_NA"].iloc[i]
                 y_before = df_runs["y_before_NA"].iloc[i]
                 y_after = df_runs["y_after_NA"].iloc[i]
                 
                 run_start = df_runs["start_frame"].iloc[i]
                 run_stop = df_runs["stop_frame"].iloc[i]
                 run_len = df_runs["run_length"].iloc[i]
                    
                    
                    
                 
                    
                #skip cases with missing boundary values
                 if  x_before is None or x_after is None or y_before is None or y_after is None:
                     boundary_cases.loc[col, ant] +=1
                     continue 
                 
                 if (run_len <= threshold_t) or (dist <=threshold_d) :
                     
                     
                       
                    #evenly spaced replacement intervall (omitting beginning and end value)
                    x_replace = np.linspace(x_before, x_after, run_len+2)[1:-1] 
                    y_replace = np.linspace(y_before, y_after, run_len+2)[1:-1]
                    assignment_flag = np.full(run_len, 5)
                    
                    # Safety check
                    assert len(x_replace) == run_len, f"Length mismatch for {col}-{ant}, run {i}"
                    idx = pd.IndexSlice
                    print(f"{col}-{ant}: run_start={run_start}, run_stop={run_stop}, run_len={run_len}")
                    print("  x_replace len:", len(x_replace))
                    print("  num rows to assign:", len(combined_colony_data_interpolated.loc[idx[col, ant, run_start:run_stop]]))
                    print(type(run_start))
                    print(type(run_start))
                    print((col, ant, run_start) in combined_colony_data_interpolated.index)
                    print((col, ant, run_start) in combined_colony_data_interpolated.index)
                    #add flag that this was interpolated 
                   
                    combined_colony_data_interpolated.loc[idx[col, ant, run_start:run_stop], "x"] = x_replace
                    combined_colony_data_interpolated.loc[idx[col, ant, run_start:run_stop], "y"] = y_replace
                    combined_colony_data_interpolated.loc[idx[col, ant, run_start:run_stop], "assignment_type"] = assignment_flag
                    
                    #combined_colony_data_interpolated.loc[(col, ant)]["x"].iloc[run_start: run_stop+1] = x_replace
                    #combined_colony_data_interpolated.loc[(col, ant)]["y"].iloc[run_start: run_stop+1] = y_replace
                    #combined_colony_data_interpolated.loc[(col, ant)]["assignment_type"].iloc[run_start: run_stop+1] = assignment_flag
                    interpolation_sum_ant = interpolation_sum_ant + run_len
                    
            
            interpolation_summary[col, ant] = interpolation_sum_ant/(total_num_frames/100)
            
    return combined_colony_data_interpolated, interpolation_summary, boundary_cases


