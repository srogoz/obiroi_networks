#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  8 14:55:53 2025

@author: srogoz@ice.mpg.de
"""

##############
##%%finding colonies
import h5py
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#initialisting root base directory
os.chdir('/home/srogoz@ice.mpg.de/ulr-lab/Users/Sarah_saro7615/EXP1_GROUPSIZEandCOMPOSITION/')

colony_directories = { "b16-1" : 'sarah_3216ab_1_0_lasthalf/b16-1/antdata',
                      "ba16-1" : 'sarah_3216ab_1_0_lasthalf/ba16-1/antdata',
                      "a16-1" : 'sarah_3216ab_1_0_lasthalf/a16-1/antdata',
                      "b16-2" : 'sarah_3216ab_2_0_lasthalf/b16-2/antdata/C1',
                      "ba16-2" : 'sarah_3216ab_2_0_lasthalf/b16-2/antdata/C2',
                      "a16-2" : 'sarah_3216ab_2_0_lasthalf/a16-2/antdata',
                      "b16-3" : 'sarah_3216ab_2_0_lasthalf/b16-3/antdata',
                      "ba16-3" : 'sarah_3216ab_2_0_lasthalf/ba16-3/antdata',
                      "a16-3" : 'sarah_3216ab_2_0_lasthalf/a16-3/antdata',
                      "b16-4" : 'sarah_3216ab_3_0_lasthalf/b16-4/antdata',
                      "ba16-4" : 'sarah_3216ab_3_0_lasthalf/ba16-4/antdata',
                      "a16-4" : 'sarah_3216ab_3_0_lasthalf/a16-4/antdata',
                      "b16-5" : 'sarah_3216ab_4_0_lasthalf/b16-5/antdata',
                      "ba16-5" : 'sarah_3216ab_4_0_lasthalf/ba16-5/antdata',
                      "a16-5" : 'sarah_3216ab_4_0_lasthalf/a16-5/antdata',
                     "b16-6" : 'sarah_3216ab_5_0_lasthalf/b16-6/antdata',
                     "ba16-6" : 'sarah_3216ab_5_0_lasthalf/ba16-6/antdata',
                     "a16-6" : 'sarah_3216ab_5_0_lasthalf/a16-6/antdata',
                     "b16-7" : 'sarah_3216ab_6_0_lasthalf/b16-7/antdata/C1',
                     "ba16-7" : 'sarah_3216ab_6_0_lasthalf/b16-7/antdata/C2',
                     "a16-7" : 'sarah_3216ab_6_0_lasthalf/2/antdata/C1',
                     "b16-8" : 'sarah_3216ab_6_0_lasthalf/2/antdata/C2',
                     "ba16-8" : 'sarah_3216ab_6_0_lasthalf/3/antdata/C1',
                     "a16-8" : 'sarah_3216ab_6_0_lasthalf/3/antdata/C2'
                                            
                      }

expected_ants = ["BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", 
                   "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"]

colony_names = ["a16-1", "a16-2", "a16-3", "a16-4", "a16-5", "a16-6", "a16-7", "a16-8", "b16-1", "b16-2", "b16-3", "b16-4", "b16-5", "b16-6", "b16-7", "b16-8",  "ba16-1", "ba16-2", "ba16-3", "ba16-4", "ba16-5", "ba16-6", "ba16-7", "ba16-8"]
#distance for interpolation based on distance is one body length of an ant and 2 seconds
nanThreshold_d = 0.002
nanThreshold_t = 20


#def read_antrax_xy_data (directory_list):
##########
#%% find and sort all .mat files
import pathlib

##for all colonies in colony directory
#for path in colony_directories.values():
 #   mat_files = list(pathlib.Path(pa.glob("*.mat"))
  #  print(mat_files)
  
  #to check if path is correct
#print(os.path.exists(colony_directories['b16-1']))


sorted_files_per_colony = {}

for colony_name in colony_names :
    ##list and sort all files in the directory, f.stem: file without extension
     ##split(x)by sign and look at first one
     
    colony_files = list(pathlib.Path(colony_directories[colony_name]).glob('*.mat'))
    colony_files_sorted = sorted(colony_files, key=lambda f: int(f.stem.split('_')[1]))
    sorted_files_per_colony [colony_name] = colony_files_sorted


####################
#%% read and merge all .mat files in combined_colony_data--WORKING ONEE!!!!
selected_colonies = ["b16-5"]

combined_colony_data = pd.DataFrame()


for col in selected_colonies:
    print(f"Working on {col}")
   
       
    colony_files_sorted = sorted_files_per_colony[col]
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
############
#%% set multi-index "colony", "ant", keeps frame index as a multi  
combined_colony_data.reset_index(inplace=True) 
combined_colony_data.set_index(["colony", "ant", "frame"], inplace=True)
combined_colony_data.sort_index(inplace= True)

#################
# %% safety check dimensions of combined_colony_data
print(combined_colony_data.reset_index().groupby(["colony", "ant"])["frame"].min())
print(combined_colony_data.reset_index().groupby(["colony", "ant"])["frame"].max())

#################
# %% safety check for indexing access
print(combined_colony_data.index.names)
idx = pd.IndexSlice
combined_colony_data.loc[idx["b16-1", "BG", slice(730,740)], "x" ]

###################
#%% save combined antrax data from all colonies 
#reset multindex data interpolated
from datetime import datetime

today= datetime.today().strftime('%Y%m%d')

filename = f"b16-2_interpo_{today}.pkl"

filepath = f"processed/{filename}"

df_reset = combined_colony_data_interpolated.reset_index()

df_reset.to_pickle(filepath)

###################
#%% save combined antrax data from all colonies 
#reset multindex data ORIGINAL
from datetime import datetime

today= datetime.today().strftime('%Y%m%d')

filename = f"b16-2_{today}.pkl"

filepath = f"processed/{filename}"

df_reset = combined_colony_data.reset_index()

df_reset.to_pickle(filepath)
########################
#%% save as parquet to be used in R
#reset multiindex, interpolated
from datetime import datetime 

today= datetime.today().strftime('%Y%m%d')

filename = f"b16-4_interpo_{today}.parquet"

folder = "processed"

file_path = os.path.join(folder, filename)

df_reset2 = combined_colony_data_interpolated.reset_index()

df_reset2.to_parquet(file_path, index = False)
#############################
#%% save as parquet to be used in R
#reset multiindex, ORIGINAL
from datetime import datetime 

today= datetime.today().strftime('%Y%m%d')

filename = f"b16-4_{today}.parquet"

folder = "processed"

file_path = os.path.join(folder, filename)

df_reset2 = combined_colony_data.reset_index()

df_reset2.to_parquet(file_path, index = False)