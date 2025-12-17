#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 17 13:09:58 2025

@author: srogoz@ice.mpg.de
"""
#%% function for all continous NAN runs--WORKING ONE
import pandas as pd

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
