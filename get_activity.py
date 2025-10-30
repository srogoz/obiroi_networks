#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 27 11:10:33 2025

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
#######
#%% import xy_data interpolated
selected_colonies = ["a16-1", "a16-2", "a16-3", "a16-4", "a16-5", "a16-6", "a16-7", "a16-8", "b16-1", "b16-2", "b16-3", "b16-4", "b16-5", "b16-6", "b16-7", "b16-8",  "ba16-1", "ba16-2", "ba16-3", "ba16-4", "ba16-5", "ba16-6", "ba16-7", "ba16-8"]

###########
#%% FUNCTIONS

os.chdir('/home/srogoz@ice.mpg.de/ulr-lab/Users/Sarah_saro7615/EXP1_GROUPSIZEandCOMPOSITION/')

folder= "processed/EXP1_colonies/"

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

nest_folder = "nests/"
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
################
#%%

#def get_centerofmass(polygons_filepath):
    
 #   with open (ploygons)
#######################
#%% get interpolated xy_files

xy_datainterpo = get_matched_files(folder, selected_colonies)

polygons = get_matched_polygons(nest_folder, selected_colonies)

##
#%%
points = np.load(polygons["a16-1"][0], allow_pickle=True)

polygon1= Polygon(points[0])
polygon2= Polygon(points[1])
polygon3= Polygon(points[2])
polygon4= Polygon(points[3])

#############
#%%
centroid1 = polygon1.centroid
centroid2 = polygon2.centroid
centroid3 = polygon3.centroid
centroid4 = polygon4.centroid
r_scale = 0.0000692181774923531
centroidx_mean= round(statistics.mean([centroid1.x, centroid2.x, centroid3.x,centroid4.x])*r_scale,7)
centroidy_mean= round(statistics.mean([centroid1.y, centroid2.y, centroid3.y,centroid4.y])*r_scale,7)

##############
#%% get rmsd
colony_dir = xy_datainterpo["a16-1"]
colony_xy = pd.read_parquet(colony_dir)

colony_xy["centroid_diff_xsq"] = ((centroidx_mean-colony_xy["x"])**2)
colony_xy["centroid_diff_ysq"] = ((centroidy_mean-colony_xy["y"])**2)
colony_xy["sum_singleframe"]=colony_xy["centroid_diff_xsq"]+colony_xy["centroid_diff_ysq"]          
####
#%%
framelimit=12*60*60*10
#colony_xy = (colony_xy.drop(columns =["colony"]).set_index(["ant","frame"]).sort_index())
ants = colony_xy.index.get_level_values("ant").unique()
frames = colony_xy.index.get_level_values("frame").unique()                                
all_frames = colony_xy.index.get_level_values("frame").max()   
start = all_frames-framelimit

#crop to 12hs 

colony_xy_crop= colony_xy[colony_xy.index.get_level_values("frame")>=start]

singleant_activity = ()

for ant in ants:
    singleactivity[ant] = colony_xy.loc[ ant, "sum_singleframe"].sum()
    #divided by number of frames in which not na 