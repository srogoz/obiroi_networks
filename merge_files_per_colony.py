#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 17 15:00:45 2025

@author: srogoz@ice.mpg.de
"""
import pandas as pd 
import h5py
import pathlib
# read and merge all .mat files in combined_colony_data--WORKING ONEE!!!!

def merge_files_for_colony(col, colony_directories, expected_ants, colony_names):
    # find and sort all .mat files per colony
   

    
         
    colony_files = list(pathlib.Path(colony_directories[col]).glob('*.mat'))
    colony_files_sorted = sorted(colony_files, key=lambda f: int(f.stem.split('_')[1]))
   

    combined_colony_data = pd.DataFrame()
    #for col in selected_colonies:
    print(f"Working on {col}")
          
    
    current_max_frame = 0 
     

    for i, file_path in enumerate(colony_files_sorted) :
        max_frames_in_this_file = 0 #frame number in this file
          
       
        
        for ant in expected_ants :
           
            with h5py.File(file_path, "r") as file:
                    #print(list(file.keys()))
                    if ant not in file:
                        continue
                    
                    
                  
                    #load the ant specific data
                    dataset = file[ant]
                    antrax_data = pd.DataFrame(dataset[()])
                    antrax_data = antrax_data.T
                    #add ant metadata
                    print(f"Shape of data for ant {ant} in colony {col} in {i}: {antrax_data.shape}")
                    if antrax_data.shape[1] == 4:
                       antrax_data.columns = ["x", "y", "orientation", "assignment_type"]
                    else:
                          print(f"Unexpected shape for {col}, {ant}, skipping...")
                          continue  # skip to the next one
                    antrax_data.columns = ["x", "y", "orientation", "assignment_type"]
                    antrax_data.index.name ="frames"
                    antrax_data ["colony"] = col
                    antrax_data["ant"] = ant
                    
                     #overall frame counter
                    num_frames = antrax_data.shape[0]
                    antrax_data.index = range(current_max_frame, current_max_frame + num_frames)
                    antrax_data.index.name = "frame"
         
                    # Track max number of frames seen
                    max_frames_in_this_file = max(max_frames_in_this_file, antrax_data.shape[0])


         
                    # Add a column to track which file it came from (optional)
                    antrax_data["source_file"] = i
                    
                    # Append to the combined dataframe
                    combined_colony_data = pd.concat([combined_colony_data, antrax_data])
        # Update for next run
        current_max_frame += max_frames_in_this_file
        
    # set multi-index "colony", "ant", keeps frame index as a multi          
    combined_colony_data.reset_index(inplace=True) 
    combined_colony_data.set_index(["colony", "ant", "frame"], inplace=True)
    combined_colony_data.sort_index(inplace= True)
    
    return combined_colony_data
