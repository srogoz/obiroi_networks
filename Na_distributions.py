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
#%% read and merge all .mat files in combined_colony_data

combined_colony_data = pd.DataFrame()


for col in colony_names:
    print(f"Working on {col}")
   
       
    colony_files_sorted = sorted_files_per_colony[col]
    current_max_frame = 0 
     

    for i in range(len(colony_files_sorted)) :
     
       
        
        for ant in expected_ants :
            file_path = colony_files_sorted[i]
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
         
                    # Update for next run
                    current_max_frame += num_frames
         
                    # Add a column to track which file it came from (optional)
                    antrax_data["source_file"] = i
                    
                    # Append to the combined dataframe
                    combined_colony_data = pd.concat([combined_colony_data, antrax_data])
    
               
############
#%%
      
#set multi-index "colony", "ant"     
combined_colony_data.set_index(["colony", "ant"], inplace=True)   

###################
#%% save combined antrax data from all colonies 
from datetime import datetime

##get date
today= datetime.today().strftime('%Y%m%d')

filename = f"combined_colony_data_UNPROCESSED_singleindex{today}.pkl"

filepath = f"processed/{filename}"

##save as .pkl in folder processed

combined_colony_data.to_pickle(filepath)



##########
#%%
combined_colony_data.loc[("a16-1", "GG"), "x"].head(100)
#################
#%%

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
    
#############################
#%% heat plot for assignment rate
import seaborn as sns
import matplotlib.pyplot as plt

assignmentrate_percolony_perant = pd.DataFrame(index=selected_colonies, columns=expected_ants)

# Fill in the assignment rate values
for (col, ant), df in combined_colony_data.groupby(level=["colony", "ant"]):
    if "x" not in df.columns:
        continue  # skip if 'x' column is missing
    x = df["x"]
    if len(x) == 0:
        continue
    nan_percent = x.isna().sum() / len(x) * 100
    assignmentrate = 100 - nan_percent  # percent of non-NaN
    assignmentrate_percolony_perant.loc[col, ant] = assignmentrate

# Convert values to float (they may be stored as object due to NaNs)
assignmentrate_percolony_perant = assignmentrate_percolony_perant.astype(float)

# Plot heatmap
plt.figure(figsize=(14, 8))
sns.heatmap(assignmentrate_percolony_perant, annot=True, fmt=".1f", cmap="viridis", cbar_kws={'label': 'Assignment Rate (%)'})
plt.title("Assignment Rate per Colony per Ant")
plt.xlabel("Ant")
plt.ylabel("Colony")
plt.tight_layout()
plt.show()
##############
# %% save assignment rates
# convert data
# %% save interpolation_summary
today= datetime.today().strftime('%Y%m%d')

filename = f"assignmentrates_b16-4_{today}.npy"

filepath = f"processed/assignment_rates_before_interpo/{filename}"

np.save(filepath, assignmentrate_percolony_perant) 
###########
################
#%% plot assignment rate old and interpolation values
plt.figure(figsize=(14, 8))
sns.heatmap(assignmentrate_percolony_perant, annot=True, fmt=".1f", cmap="viridis", cbar_kws={'label': 'Assignment Rate (%)'})
plt.title("Assignment Rate per Colony per Ant")
plt.xlabel("Ant")
plt.ylabel("Colony")
plt.tight_layout()
plt.show()

#############
#%% check average for assignment rate per colony
print(assignmentrate_percolony_perant.loc["a16-1"].mean())
#########
#%% function for all continous NAN runs--WORKING ONE
def find_continuous_NAN_runs(df, colony, ant):
    runs = []
    x_starts = []
    x_stops = []
    y_starts = []
    y_stops = []
    start_frames = []
    stop_frames = []
    distance = []
    count = 0 
    start = None

    for i, val in enumerate(df["x"].isna()):
        if val:
            if count == 0:
                start = i 
            count += 1
        elif count > 0:
            
            #if there is an End of Nan in frames
            runs.append(count)
            start_frames.append(start)
            stop = start + count - 1
            stop_frames.append(stop)
            
            #check if NA value is part of a run
            pre_index = start - 1
            post_index = start + count

            x_before = (df.iloc[pre_index]['x'] if pre_index >= 0 else None)
            y_before = (df.iloc[pre_index]['y'] if pre_index >= 0 else None)
            x_after = (df.iloc[post_index]['x'] if post_index < len(df) else None)
            y_after = (df.iloc[post_index]['y'] if post_index < len(df) else None)
            
            x_starts.append(x_before)
            y_starts.append(y_before)
            x_stops.append(x_after)
            y_stops.append(y_after)
            
            if x_before is not None and x_after is not None and y_before is not None and y_after is not None:
                dist = ((abs(x_before-x_after)**2 + abs(y_before-y_after)**2)**(1/2))
               
            else:
                dist = None
            
            distance.append(dist)
            
            count = 0
    #if nan run til end 
    if count > 0: 
        
        stop = start + count - 1
        runs.append(count)
        start_frames.append(start)
        stop_frames.append(stop)

        pre_index = start - 1
        post_index = stop + 1

        x_before = df.iloc[pre_index]['x'] if pre_index >= 0 else None
        y_before = df.iloc[pre_index]['y'] if pre_index >= 0 else None
        x_after = df.iloc[post_index]['x'] if post_index < len(df) else None
        y_after = df.iloc[post_index]['y'] if post_index < len(df) else None

        x_starts.append(x_before)
        y_starts.append(y_before)
        x_stops.append(x_after)
        y_stops.append(y_after)

        if x_before is not None and x_after is not None and y_before is not None and y_after is not None:
            dist = ((x_after - x_before)**2 + (y_after - y_before)**2)**0.5
        else:
            dist = None
            
        distance.append(dist)
        
        

    result_df = pd.DataFrame({
        "run_length" : runs, 
        "start_frame": start_frames,
        "stop_frame": stop_frames,
        "x_before_NA": x_starts,
        "x_after_NA" : x_stops,
        "y_before_NA": y_starts,
        "y_after_NA": y_stops,
        "distance" : distance
    })

    # 🔹 Add identifiers
    result_df["colony"] = colony
    result_df["ant"] = ant

    return result_df

###########
#%%
all_nan_runs = []

for (col, ant), df in combined_colony_data.groupby(level=["colony", "ant"]):
    result = find_continuous_NAN_runs(df, col, ant)
    all_nan_runs.append(result)

# Concatenate just once
result_continuous_nan_runs = pd.concat(all_nan_runs, ignore_index=True)
########
#%%assign multiindex
result_continuous_nan_runs.set_index(["colony", "ant"], inplace=True)   


##################
#%%

from datetime import datetime 

today= datetime.today().strftime('%Y%m%d')

filename = f"result_continuous_na_runs_b16-1_{today}.pkl"

filepath = f"processed/{filename}"

##save as .pkl in folder processed


result_continuous_nan_runs.to_pickle(filepath)
#############
#%% save as parquet to be used in R
from datetime import datetime 
#reset multiindex


today= datetime.today().strftime('%Y%m%d')
filename = f"b16-1_{today}.parquet"
#reset multiindex

df_reset = combined_colony_data.reset_index()

df_reset.to_parquet("processed/b16", index = False)
#############
#%%

print(combined_colony_data.index)

############
#%%



import reassign_nan_runs
selected_colonies = ["b16-4"]

combined_colony_data_interpolated, interpolation_summary, boundary_cases = reassign_nan_runs.reassign_nan_runs(
    combined_colony_data,
    result_continuous_nan_runs,
    selected_colonies,
    expected_ants,
    nanThreshold_d,
    nanThreshold_t
)

###############
# %%
runlength_difference = []

for i in range(len(result_continuous_nan_runs["run_length"])):
    run = result_continuous_nan_runs["run_length"].iloc[i]
    stop = result_continuous_nan_runs[ "stop_frame"].iloc[i]
    start = result_continuous_nan_runs[ "start_frame"].iloc[i]
    dif = stop-start+1
    
    if run == dif :
        runlength_difference.append(True)
    else :
        runlength_difference.append(False)
 
############
#%% they are all the same length
runlength_difference.count('False')

#############
# %%

#%% check if frames are missing 
runlength_difference = []

for runs in range(enumerate(result_continuous_nan_runs.loc["run_length"])):
    run = result_continuous_nan_runs["run_length"].iloc[runs]
    stop = result_continuous_nan_runs[ "stop_frame"].iloc[runs]
    start = result_continuous_nan_runs[ "start_frame"].iloc[runs]
    dif = stop-start
    if run == dif 
        runlength_difference.append(True)
        else runlength_difference.append(False)

#######################
#%%run check
combined_colony_data= combined_colony_data.sort_index
combined_colony_data.loc[pd.IndexSlice["b16-1", "GB", 773:776]]
##############
#%%
print(combined_colony_data.index)
print(combined_colony_data.head())
print("Index sample:", combined_colony_data.index.names)
print("Unique colonies:", combined_colony_data.index.get_level_values('colony').unique())
print("Unique ants:", combined_colony_data.index.get_level_values('ant').unique())

#############
# %% adjust for multindexing leading to 1-16 being NA

interpolation_summary_2 = interpolation_summary.iloc[0, 16:32]
#######
# %% save interpolation_summary
today= datetime.today().strftime('%Y%m%d')

filename = f"interpolation_summary_b16-4_{today}.npy"

filepath = f"processed/interpolation_summary/{filename}"

np.save(filepath, interpolation_summary_2) 



#%% plot histograms and plot 
#Plot histogram of NA run lengths with colorcode


distance = ((abs(df['x_after_NA']-df['x_before_NA']))**2 + abs(df['y_after_NA']-df['y_before_NA'])**2)**(1/2)
conditionlist = [ (distance<=0.0005) & (df['run_length']<10) ,
                  (distance <=0.0005) & (df['run_length'] >10)
                  ]
df ['distance'] = ((abs(df['x_after_NA']-df['x_before_NA']))**2 + abs(df['y_after_NA']-df['y_before_NA'])**2)**(1/2)
#color_list = ['green', 'orange']

#result_continuous_nan_runs.loc["a16-1"]['position_measure'] = np.select(conditionlist, color_list, default="red")
##########
#%% distance distribution
#print(df.max())
#df = result_continuous_nan_runs[result_continuous_nan_runs["colony"] == "a16-1"]
df["distance"]= df["distance"].astype(float).round(4)
df1= df[df["ant"]=="PG"]
# Define bin width (e.g., 0.001)
bin_width = 0.001
bins = np.arange(0, df1['distance'].max() + bin_width, bin_width)

plt.figure(figsize=(8, 4))
plt.hist(df1['distance'],
         bins=bins,
         color = "springgreen"
        # edgecolor = "pink"
       #  edgecolor= a['position_measure']
         )
plt.xlabel("Distance [m]")
plt.ylabel("Frequency")
plt.title("Distribution of distance at beginning and end of NaN runs for PG")
plt.xlim(0, 0.01)
plt.ylim(0,6700)
plt.tight_layout()
plt.show()
#################################
#%% run length distribution-used

#df = result_continuous_nan_runs[result_continuous_nan_runs["colony"] == "a16-1"]
plt.figure(figsize=(8, 4))
plt.hist(df['run_length'],
         bins=range(1, df['run_length'].max() + 2),
         edgecolor = "purple"
       #  edgecolor= a['position_measure']
         )
plt.xlabel("Length of continuous NaN sequences")
plt.ylabel("Frequency")
plt.title("Distribution of Continuous NaN Runs")
plt.xlim(0, 500)
plt.ylim(0,7000)
plt.tight_layout()
plt.show()

########################
#%%
interpolated_assignment_rates = pd.DataFrame(index=colony_names, columns=expected_ants)

colony = "a16-1"
total_frames = 466973
for ant in expected_ants:
    condition_0 = (
        (df["colony"] == colony) &
        (df["ant"] == ant) &
        (df["distance"] < 0.002) &
        (df["run_length"] < 20)
    )
    condition_d = (
        (df["colony"] == colony) &
        (df["ant"] == ant) &
        (df["distance"] < 0.002)
    )
    condition_t = (
        (df["colony"] == colony) &
        (df["ant"] == ant) &
        (df["run_length"] < 20)
    )

    ass_0 = df.loc[condition_0, "run_length"].sum()
    ass_d = df.loc[condition_d, "run_length"].sum()
    ass_t = df.loc[condition_t, "run_length"].sum()

    result = (((ass_t - ass_0) + (ass_d - ass_0))) / (total_frames / 100)
    interpolated_assignment_rates.loc[colony, ant] = result
###############
#%%  does not work yettt
assign_reduced = assignmentrate_percolony_perant.loc[["a16-1"]].astype(float)

interpolated = interpolated_assignment_rates.loc[["a16-1"]].astype(float)
interpolated.index = ["a16-1_interpolated"]

sum_df = assign_reduced + interpolated
sum_df.index = ["a16-1_total"]

combined_df = pd.concat([assign_reduced, interpolated, sum_df])

plt.figure(figsize=(14, 6))
sns.heatmap(combined_df, annot=True, fmt=".1f", cmap="viridis", cbar_kws={'label': 'Assignment Rate (%)'})
plt.title("Assignment Rates for Colony a16-1 (Original, Interpolated, Total)")
plt.xlabel("Ant")
plt.ylabel("Rate Type")
plt.tight_layout()
plt.show()
###############
#%% check how high assignment rate would be if remove all runs smaller then 2s and below 2mm
condition_0 = ((df["colony"] == "a16-1")&
               (df["ant"] == ant)&
             (df["distance"] < 0.002) &
             (df["run_length"] < 20)
             )
condition_d = ((df["colony"] == "a16-1")&
               (df["ant"] == ant)&
             (df["distance"] < 0.002) 
              )

condition_t = ((df["colony"] == "a16-1")&
               (df["ant"] == ant)&
               (df["run_length"] < 20)
               )

print(condition.sum())  # number of True matches
print(df.loc[condition, "run_length"].sum())  # the rows that match
combined_colony_data.sort_index(inplace=True)
print(len(combined_colony_data.loc[("a16-1", "GB"),"x"]))
print((combined_colony_data.loc[("a16-1", "GB"),"x"]).isna().sum())
#print(result_continuous_nan_runs.loc[condition, "run_length"].sum())
nas_per_ant = combined_colony_data["x"].isna().groupby(level=["colony", "ant"]).sum()
print(nas_per_ant)
############################
#%% max of nan runlength

print(result_continuous_nan_runs[result_continuous_nan_runs["colony"]=="a16-1"].max())
 

