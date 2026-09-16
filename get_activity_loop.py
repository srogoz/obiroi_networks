#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 09:29:00 2025

@author: srogoz@ice.mpg.de
"""

#################
#%%import libraries

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import re
from shapely.geometry import Polygon
import pickle 
import statistics
import math
from datetime import datetime
#######
#%% import DATA EXP1 
selected_colonies = ["a16-1", "a16-2", "a16-3", "a16-4", "a16-5", "a16-6", "a16-7", "a16-8", "b16-1", "b16-2", "b16-3", "b16-4", "b16-6", "b16-7", "b16-8",  "ba16-1", "ba16-2", "ba16-3", "ba16-4", "ba16-5", "ba16-6", "ba16-7", "ba16-8"]

r_scale_list = {"a16-1" : 0.0000692181774923531, 
                "a16-2" : 6.5588677500214875E-5, 
                "a16-3" : 6.6262348839314609E-5,
                "a16-4" : 6.76039164635912E-5,
                "a16-5" : 6.9075080950827393E-5,
                "a16-6" : 6.9196848745745241E-5,
                "a16-7": 6.4840496848792E-5,
                "a16-8": 6.4997579704194046E-5,
                "b16-1": 6.8428691920128E-5,
                "b16-2": 6.5813747884814654E-5,
                "b16-3": 6.6461533088672529E-5,
                "b16-4": 6.7381187955876053E-5,
                "b16-6" : 7.014591415803833E-5,
                "b16-7" : 6.50545992694844E-5,
                "b16-8": 6.4840496848792E-5,
                "ba16-1" : 6.8392404639652343E-5,
                "ba16-2" : 6.5813747884814654E-5,
                "ba16-3" : 6.5900149289449417E-5,
                "ba16-4" : 6.79059203962662E-5,
                "ba16-5" : 6.8976056111590375E-5,
                "ba16-6" : 6.8936356959651283E-5,
                "ba16-7" : 6.50545992694844E-5,
                "ba16-8" : 6.4997579704194046E-5}

expected_ants = ["BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", 
                   "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"]

tag_genotype_mixedcolony = {"BB": "b",
                      "BG": "b",
                      "BO": "b",
                      "BP": "b",
                      "GB": "b",
                      "GG": "b",
                      "GO": "b",
                      "GP": "b",
                      "OB": "a",
                      "OG": "a",
                      "OO": "a",
                      "OP": "a",
                      "PB": "a",
                      "PG": "a",
                      "PO": "a",
                      "PP": "a"}

###########
#%% EXP1 beforefungus, 1-6_2_last

selected_colonies = ["a16-1", "a16-2", "a16-3", "a16-4", "a16-5", "a16-6", "a16-7", "a16-8", "b16-1", "b16-2", "b16-3", "b16-4", "b16-5", "b16-6", "b16-7", "b16-8",  "ba16-1", "ba16-2", "ba16-3", "ba16-4", "ba16-5", "ba16-6", "ba16-7", "ba16-8"]

r_scale_list = {"a16-1" : 6.9230769230769237E-5, 
                "a16-2" : 6.61764705882353E-5, 
                "a16-3" : 6.61764705882353E-5,
                "a16-4" : 6.797583081570997E-5,
                "a16-5" : 6.6824466892463082E-5,
                "a16-6" : 6.9230769230769237E-5,
                "a16-7": 6.5789473684210525E-5,
                "a16-8": 6.61764705882353E-5,
                "b16-1": 6.81818181818182E-5,
                "b16-2": 6.61764705882353E-5,
                "b16-3": 6.61764705882353E-5,
                "b16-4": 6.8389057750759875E-5,
                "b16-5": 6.8962813944831221E-5,
                "b16-6" : 6.9230769230769237E-5,
                "b16-7" : 6.5789473684210525E-5,
                "b16-8": 6.5789473684210525E-5,
                "ba16-1" : 6.8807339449541291E-5,
                "ba16-2" : 6.61764705882353E-5,
                "ba16-3" : 6.61764705882353E-5,
                "ba16-4" : 6.797583081570997E-5,
                "ba16-5" : 6.812398025488937E-5,
                "ba16-6" : 6.8807339449541291E-5,
                "ba16-7" : 6.5789473684210525E-5,
                "ba16-8" : 6.61764705882353E-5}

expected_ants = ["BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", 
                   "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"]

tag_genotype_mixedcolony = {"BB": "b",
                      "BG": "b",
                      "BO": "b",
                      "BP": "b",
                      "GB": "b",
                      "GG": "b",
                      "GO": "b",
                      "GP": "b",
                      "OB": "a",
                      "OG": "a",
                      "OO": "a",
                      "OP": "a",
                      "PB": "a",
                      "PG": "a",
                      "PO": "a",
                      "PP": "a"}
#############
#%% after exposure 
selected_colonies = ["a16-1", "a16-2", "a16-3", "a16-4", "a16-5", "a16-6", "a16-7", "a16-8", "b16-1", "b16-2", "b16-3", "b16-4", "b16-5" , "b16-6", "b16-7", "b16-8",  "ba16-1", "ba16-2", "ba16-3", "ba16-4", "ba16-5", "ba16-6", "ba16-7", "ba16-8"]

r_scale_list = {"a16-1" : 6.9230769230769237E-5, 
                "a16-2" : 6.61764705882353E-5, 
                "a16-3" : 6.61764705882353E-5,
                "a16-4" : 6.797583081570997E-5,
                "a16-5" : 6.6824466892463082E-5,
                "a16-6" : 6.9230769230769237E-5,
                "a16-7": 6.5789473684210525E-5,
                "a16-8": 6.61764705882353E-5,
                "b16-1": 6.81818181818182E-5,
                "b16-2": 6.61764705882353E-5,
                "b16-3": 6.61764705882353E-5,
                "b16-4": 6.8389057750759875E-5,
                "b16-5": 6.8962813944831221E-5,
                "b16-6" : 6.9230769230769237E-5,
                "b16-7" : 6.5789473684210525E-5,
                "b16-8": 6.5789473684210525E-5,
                "ba16-1" : 6.8807339449541291E-5,
                "ba16-2" : 6.61764705882353E-5,
                "ba16-3" : 6.61764705882353E-5,
                "ba16-4" : 6.797583081570997E-5,
                "ba16-5" : 6.812398025488937E-5,
                "ba16-6" : 6.8807339449541291E-5,
                "ba16-7" : 6.5789473684210525E-5,
                "ba16-8" : 6.61764705882353E-5}

expected_ants = ["BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", 
                   "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"]

tag_genotype_mixedcolony = {"BB": "b",
                      "BG": "b",
                      "BO": "b",
                      "BP": "b",
                      "GB": "b",
                      "GG": "b",
                      "GO": "b",
                      "GP": "b",
                      "OB": "a",
                      "OG": "a",
                      "OO": "a",
                      "OP": "a",
                      "PB": "a",
                      "PG": "a",
                      "PO": "a",
                      "PP": "a"}
###########
#%% FUNCTIONS

os.chdir('/home/srogoz@ice.mpg.de/ulr-lab/Users/Sarah_saro7615/EXP1_GROUPSIZEandCOMPOSITION/')

folder= "processed/EXP1_colonies_inf/"
nest_folder = "nests/afterfungus/"

def get_matched_files(folder, selected_colonies):
   #
    file_list = {colony : [] for colony in selected_colonies}
    
    for file in os.listdir(folder):
    
        
        for colony in selected_colonies:
            pattern = re.compile(rf"^{re.escape(colony)}_interpo.*")
            if pattern.match(file):
                file_path = os.path.join(folder, file)
                file_list[colony].append(file_path)
                            
                        
    return file_list


def get_matched_polygons(folder, selected_colonies):
   #folder up to nests
    file_list = {colony : [] for colony in selected_colonies}
    for root, dirs, files in os.walk(folder):
        base=os.path.basename(root)
        parent = os.path.basename(os.path.dirname(root))
        
        if base== "masks" and parent in selected_colonies:
                for f in files:
                    if f in ["all_points.npy"]:
                        file_path = os.path.join(root, f)
                        file_list[parent].append(file_path)
                                    
                        
    return file_list

#######################
#%% get interpolated xy_file and nest polygons paths

xy_datainterpo = get_matched_files(folder, selected_colonies)

polygons = get_matched_polygons(nest_folder, selected_colonies)

##
#%%WORKING calculation of activity!!!


singleant_activity =pd.DataFrame(index=expected_ants, columns= selected_colonies )
framelimit=12*60*60*10


for colony in selected_colonies:
  
    #colony = "b16-2"
    points = np.load(polygons[colony][0], allow_pickle=True)
    
    #get center of mass for polygons for colony
    polygon1= Polygon(points[0])
    polygon2= Polygon(points[1])
    polygon3= Polygon(points[2])
    polygon4= Polygon(points[3])


    centroid1 = polygon1.centroid
    centroid2 = polygon2.centroid
    centroid3 = polygon3.centroid
    centroid4 = polygon4.centroid
    r_scale = r_scale_list[colony]
    centroidx_mean= round(statistics.mean([centroid1.x, centroid2.x, centroid3.x,centroid4.x])*r_scale,7)
    centroidy_mean= round(statistics.mean([centroid1.y, centroid2.y, centroid3.y,centroid4.y])*r_scale,7)

    #get colony data
    colony_dir = xy_datainterpo[colony]
    colony_xy = pd.read_parquet(colony_dir)
    
    #x_frame= x_i -x_centroid and y_frame= y_i-y_centroid and sum x_frame + y_frame over single frame
    colony_xy["centroid_diff_xsq"] = ((centroidx_mean-colony_xy["x"])**2)
    colony_xy["centroid_diff_ysq"] = ((centroidy_mean-colony_xy["y"])**2)
    colony_xy["sum_singleframe"]=colony_xy["centroid_diff_xsq"]+colony_xy["centroid_diff_ysq"]          

    
    colony_xy = (colony_xy.drop(columns =["colony"]).set_index(["ant","frame"]).sort_index())
    ants = colony_xy.index.get_level_values("ant").unique()
    frames = colony_xy.index.get_level_values("frame").unique()                                
    all_frames = colony_xy.index.get_level_values("frame").max()   
    start = all_frames-framelimit

    #crop to 12hs 
    
    colony_xy_crop= colony_xy[colony_xy.index.get_level_values("frame")>=start]



    for ant in ants:
        
        #divided by number of frames in which not na 
        ant_data=colony_xy_crop.loc[ ant, 'sum_singleframe']
        
        valid_data= ant_data.dropna()
        all_frames_sum=valid_data.sum()
        all_frames=len(valid_data)
       
        
        rmsd_ant = math.sqrt(all_frames_sum/all_frames)
        
        singleant_activity.loc[ant, colony]=rmsd_ant
        
    print(f"{colony}")

#convert to unit: mm of rmsd

singleant_activity_mm= singleant_activity*1000


##
#%%WORKING calculation of activity for infection timeframe, meaning only one nest per timeintervall!!!


singleant_activity =pd.DataFrame(index=expected_ants, columns= selected_colonies )
framelimit=12*60*60*10


for colony in selected_colonies:
  
    #colony = "b16-2"
    points = np.load(polygons[colony][0], allow_pickle=True)
    
    #get center of mass for polygons for colony
    polygon1= Polygon(points[0])
    


    centroid1 = polygon1.centroid
   
    r_scale = r_scale_list[colony]
    centroidx_mean= round(statistics.mean([centroid1.x])*r_scale,7)
    centroidy_mean= round(statistics.mean([centroid1.y])*r_scale,7)

    #get colony data
    colony_dir = xy_datainterpo[colony]
    colony_xy = pd.read_parquet(colony_dir)
    
    #x_frame= x_i -x_centroid and y_frame= y_i-y_centroid and sum x_frame + y_frame over single frame
    colony_xy["centroid_diff_xsq"] = ((centroidx_mean-colony_xy["x"])**2)
    colony_xy["centroid_diff_ysq"] = ((centroidy_mean-colony_xy["y"])**2)
    colony_xy["sum_singleframe"]=colony_xy["centroid_diff_xsq"]+colony_xy["centroid_diff_ysq"]          

    
    colony_xy = (colony_xy.drop(columns =["colony"]).set_index(["ant","frame"]).sort_index())
    ants = colony_xy.index.get_level_values("ant").unique()
    frames = colony_xy.index.get_level_values("frame").unique()                                
    all_frames = colony_xy.index.get_level_values("frame").max()   
    start = all_frames-framelimit

    #crop to 12hs 
    
    colony_xy_crop= colony_xy[colony_xy.index.get_level_values("frame")>=start]



    for ant in ants:
        
        #divided by number of frames in which not na 
        ant_data=colony_xy_crop.loc[ ant, 'sum_singleframe']
        
        valid_data= ant_data.dropna()
        all_frames_sum=valid_data.sum()
        all_frames=len(valid_data)
       
        
        rmsd_ant = math.sqrt(all_frames_sum/all_frames)
        
        singleant_activity.loc[ant, colony]=rmsd_ant
        
    print(f"{colony}")

#convert to unit: mm of rmsd

singleant_activity_mm= singleant_activity*1000

########
########
#%% save as csv

today= datetime.today().strftime('%Y%m%d')

np.save(f"processed/colony_activity/activity_colonies_mm_{today}", singleant_activity_mm)

singleant_activity_mm.to_csv(f"processed/colony_activity/activity_colonies_mm_{today}.csv", index= True)



############
#%% colony_activity_mean and sd as eaquivalent to division of labour

colony_activity_mean = pd.DataFrame(
    [ {colony: singleant_activity_mm[colony].mean() for colony in selected_colonies} ]
)

colony_activity_sd = pd.DataFrame(
    [ {colony: singleant_activity_mm[colony].std() for colony in selected_colonies} ]
)

#########
#%% save mean and sd 

#save mean in .npy 
today= datetime.today().strftime('%Y%m%d')

filename = f"activity_colonymean_mm_{today}.npy"

filepath = f"processed/colony_activity/{filename}"

np.save(filepath,colony_activity_mean) 

today= datetime.today().strftime('%Y%m%d')

filename1 = f"activity_colonymean_mm_{today}"

filepath1 = f"processed/colony_activity/{filename1}"

colony_activity_mean.to_csv(f"{filepath1}.csv", index= True)

#save sd in .npy and .csv

today= datetime.today().strftime('%Y%m%d')

filename2 = f"activity_colonysd_mm_{today}.npy"

filepath2 = f"processed/colony_activity/{filename2}"

np.save(filepath2,colony_activity_sd) 

today= datetime.today().strftime('%Y%m%d')

filename3 = f"activity_colonysd_mm_{today}"

filepath3 = f"processed/colony_activity/{filename3}"

colony_activity_sd.to_csv(f"{filepath3}.csv", index= True)


##########
#%% woorking with single individual values
mixed_colonies = ["ba16-1", "ba16-2", "ba16-3", "ba16-4", "ba16-5", "ba16-6", "ba16-7", "ba16-8"]

mean_a = pd.DataFrame(columns=mixed_colonies, index=["mean_a"])
mean_b = pd.DataFrame(columns=mixed_colonies, index=["mean_b"])

a_values = pd.DataFrame(index=expected_ants, columns=mixed_colonies)
b_values = pd.DataFrame(index=expected_ants, columns=mixed_colonies)

for colony in mixed_colonies:

    a_values_loop = []
    b_values_loop = []

    for ant in expected_ants:
        genotype = tag_genotype_mixedcolony[ant]

        value = singleant_activity_mm[colony].loc[ant]

        if genotype == "a":
            a_values_loop.append(value)
            a_values.loc[ant, colony] = value

        elif genotype == "b":
            b_values_loop.append(value)
            b_values.loc[ant, colony] = value

    # store mean values
    if len(a_values_loop) > 0:
        mean_a.loc["mean_a", colony] = np.nanmean(a_values_loop)

    if len(b_values_loop) > 0:
        mean_b.loc["mean_b", colony] = np.nanmean(b_values_loop)
    

##########
#%%   
split_means=pd.concat([mean_a,mean_b])


np.save(f"processed/colony_activity/split_means_{today}.npy", split_means)
split_means.to_csv(f"processed/colony_activity/split_means_{today}.csv", index= True)
