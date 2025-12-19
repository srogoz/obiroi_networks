#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 13:19:06 2025

@author: srogoz@ice.mpg.de
"""

import cv2
import numpy as np
import os
import json
################
#%%collect  EXP1 inf background paths:
os.chdir('/home/srogoz@ice.mpg.de/ulr-lab/Users/Sarah_saro7615/EXP1_GROUPSIZEandCOMPOSITION/')
#list of bg_path_list of missing colony areas in EXP1
bg_path_list= { "a16-7":'sarah_3216ab_6_inf_20240626_161959/a16-7/parameters/backgrounds/background.png',
                "a16-8":'sarah_3216ab_6_inf_20240626_161959/a16-8/parameters/backgrounds/background.png',
                "b16-7":'sarah_3216ab_6_inf_20240626_161959/b16-7/parameters/backgrounds/background.png',
                "b16-8": 'sarah_3216ab_6_inf_20240626_161959/b16-8/parameters/backgrounds/background.png',
                "ba16-7": 'sarah_3216ab_6_inf_20240626_161959/ba16-7/parameters/backgrounds/background.png',
                "ba16-8": 'sarah_3216ab_6_inf_20240626_161959/ba16-8/parameters/backgrounds/background.png',
                "a16-6": 'sarah_3216ab_5_inf_20240626_161922/a16-6/parameters/backgrounds/background.png',
                "ba16-6": 'sarah_3216ab_5_inf_20240626_161922/ba16-6/parameters/backgrounds/background.png'
                    
   }

    
mask_path_list= { "a16-7":'sarah_3216ab_6_inf_20240626_161959/a16-7/parameters/masks/',
                "a16-8": 'sarah_3216ab_6_inf_20240626_161959/a16-8/parameters/masks/',
                "b16-7": 'sarah_3216ab_6_inf_20240626_161959/b16-7/parameters/masks/',
                "b16-8": 'sarah_3216ab_6_inf_20240626_161959/b16-8/parameters/masks/',
                "ba16-7": 'sarah_3216ab_6_inf_20240626_161959/ba16-7/parameters/masks/',
                "ba16-8": 'sarah_3216ab_6_inf_20240626_161959/ba16-8/parameters/masks/',
                "a16-6": 'sarah_3216ab_5_inf_20240626_161922/a16-6/parameters/masks/',
                "ba16-6": 'sarah_3216ab_5_inf_20240626_161922/ba16-6/parameters/masks/'
                    
   }

json_path_list= { "a16-7":'sarah_3216ab_6_inf_20240626_161959/a16-7/parameters/prmtrs.json',
                "a16-8": 'sarah_3216ab_6_inf_20240626_161959/a16-8/parameters/prmtrs.json',
                "b16-7": 'sarah_3216ab_6_inf_20240626_161959/b16-7/parameters/prmtrs.json',
                "b16-8": 'sarah_3216ab_6_inf_20240626_161959/b16-8/parameters/prmtrs.json',
                "ba16-7": 'sarah_3216ab_6_inf_20240626_161959/ba16-7/parameters/prmtrs.json',
                "ba16-8": 'sarah_3216ab_6_inf_20240626_161959/ba16-8/parameters/prmtrs.json',
                "a16-6": 'sarah_3216ab_5_inf_20240626_161922/a16-6/parameters/prmtrs.json',
                "ba16-6": 'sarah_3216ab_5_inf_20240626_161922/ba16-6/parameters/prmtrs.json'
                    
   }

########
#%% EXP2 timeframe A paths!!!!

os.chdir('/home/srogoz@ice.mpg.de/ulr-lab/Users/Sarah_saro7615/EXP2/')
#list of bg_path_list of missing colony areas in EXP1
bg_path_list= { #cam3:
                "bA16-3":'sarah_16bba_3_A_data/bA16-3/parameters/backgrounds/background.png',
                "b16-4":'sarah_16bba_3_A_data/b16-4/parameters/backgrounds/background.png',
                "B16-4":'sarah_16bba_3_A_data/bigb16-4/parameters/backgrounds/background.png',
                "bB16-4": 'sarah_16bba_3_A_data/bB16-4/parameters/backgrounds/background.png',
                "bA16-4": 'sarah_16bba_3_A_data/bA16-4/parameters/backgrounds/background.png',
                "b16-5": 'sarah_16bba_3_A_data/b16-5/parameters/backgrounds/background.png',
                #cam4:
                "B16-6":'sarah_16bba_4_A_data/bigb16-6/parameters/backgrounds/background.png',
                "b16-6":'sarah_16bba_4_A_data/b16-6/parameters/backgrounds/background.png',
                "bA16-5":'sarah_16bba_4_A_data/bA16-5/parameters/backgrounds/background.png',
                "bB16-5": 'sarah_16bba_4_A_data/bB16-5/parameters/backgrounds/background.png',
                "B16-5": 'sarah_16bba_4_A_data/bigb-5/parameters/backgrounds/background.png',
                 #cam5:
                "bA16-8":'sarah_16bba_5_A_data/bA16-8/parameters/backgrounds/background.png',
                "Bb16-8":'sarah_16bba_5_A_data/Bb16-8/parameters/backgrounds/background.png',
                "B16-8":'sarah_16bba_5_A_data/B16-8/parameters/backgrounds/background.png',
                "b16-8": 'sarah_16bba_5_A_data/smallb16-8/parameters/backgrounds/background.png',
                #cam6 
                 "bA16-7":'sarah_16bba_6_A_data/bA16-7/parameters/backgrounds/background.png',
                 "B16-7":'sarah_16bba_6_A_data/bigB16-7/parameters/backgrounds/background.png',
                 "Bb16-7":'sarah_16bba_6_A_data/Bb16-7/parameters/backgrounds/background.png',
                 "b16-7": 'sarah_16bba_6_A_data/b16-7/parameters/backgrounds/background.png',
                 "bA16-6": 'sarah_16bba_6_A_data/bA16-6/parameters/backgrounds/background.png',
                 "bB16-6": 'sarah_16bba_6_A_data/bB16-6/parameters/backgrounds/background.png',
                 
                  #cam2:
                  "Bb16-2":'sarah_16bba_2_A_data/Bb16-2/parameters/backgrounds/background.png',
                  "B16-2":'sarah_16bba_2_A_data/bigb16-2/parameters/backgrounds/background.png',
                  "b16-3":'sarah_16bba_2_A_data/b16-3/parameters/backgrounds/background.png',
                  "bA16-2": 'sarah_16bba_2_A_data/bA16-2/parameters/backgrounds/background.png',
                  "bB16-3": 'sarah_16bba_2_A_data/bB16-3/parameters/backgrounds/background.png',
                  "B16-3": 'sarah_16bba_2_A_data/bigb16-3/parameters/backgrounds/background.png'
                   
   }
    
mask_path_list= { #cam3:
                "bA16-3":'sarah_16bba_3_A_data/bA16-3/parameters/masks/',
                "b16-4":'sarah_16bba_3_A_data/b16-4/parameters/masks/',
                "B16-4":'sarah_16bba_3_A_data/bigb16-4/parameters/masks/',
                "bB16-4": 'sarah_16bba_3_A_data/bB16-4/parameters/masks/',
                "bA16-4": 'sarah_16bba_3_A_data/bA16-4/parameters/masks/',
                "b16-5": 'sarah_16bba_3_A_data/b16-5/parameters/masks/',
                 #cam4:
                 "B16-6":'sarah_16bba_4_A_data/bigb16-6/parameters/masks/',
                 "b16-6":'sarah_16bba_4_A_data/b16-6/parameters/masks/',
                 "bA16-5":'sarah_16bba_4_A_data/bA16-5/parameters/masks/',
                 "bB16-5": 'sarah_16bba_4_A_data/bB16-5/parameters/masks/',
                 "B16-5": 'sarah_16bba_4_A_data/bigb-5/parameters/masks/',
                  #cam5:
                 "bA16-8":'sarah_16bba_5_A_data/bA16-8/parameters/masks/',
                 "Bb16-8":'sarah_16bba_5_A_data/Bb16-8/parameters/masks/',
                 "B16-8":'sarah_16bba_5_A_data/B16-8/parameters/masks/',
                 "b16-8": 'sarah_16bba_5_A_data/smallb16-8/parameters/masks/',
                 #cam6 
                  "bA16-7":'sarah_16bba_6_A_data/bA16-7/parameters/masks/',
                  "B16-7":'sarah_16bba_6_A_data/bigB16-7/parameters/masks/',
                  "Bb16-7":'sarah_16bba_6_A_data/Bb16-7/parameters/masks/',
                  "b16-7": 'sarah_16bba_6_A_data/b16-7/parameters/masks/',
                  "bA16-6": 'sarah_16bba_6_A_data/bA16-6/parameters/masks/',
                  "bB16-6": 'sarah_16bba_6_A_data/bB16-6/parameters/masks/',
                  
                   #cam2:
                   "Bb16-2":'sarah_16bba_2_A_data/Bb16-2/parameters/masks/',
                   "B16-2":'sarah_16bba_2_A_data/bigb16-2/parameters/masks/',
                   "b16-3":'sarah_16bba_2_A_data/b16-3/parameters/masks/',
                   "bA16-2": 'sarah_16bba_2_A_data/bA16-2/parameters/masks/',
                   "bB16-3": 'sarah_16bba_2_A_data/bB16-3/parameters/masks/',
                   "B16-3": 'sarah_16bba_2_A_data/bigb16-3/parameters/masks/'
                     
   }

json_path_list= { #cam3:
                "bA16-3":'sarah_16bba_3_A_data/bA16-3/parameters/prmtrs.json',
                "b16-4":'sarah_16bba_3_A_data/b16-4/parameters/prmtrs.json',
                "B16-4":'sarah_16bba_3_A_data/bigb16-4/parameters/prmtrs.json',
                "bB16-4": 'sarah_16bba_3_A_data/bB16-4/parameters/prmtrs.json',
                "bA16-4": 'sarah_16bba_3_A_data/bA16-4/parameters/prmtrs.json',
                "b16-5": 'sarah_16bba_3_A_data/b16-5/parameters/prmtrs.json',
                #cam4:
                "B16-6":'sarah_16bba_4_A_data/bigb16-6/parameters/prmtrs.json',
                "b16-6":'sarah_16bba_4_A_data/b16-6/parameters/prmtrs.json',
                "bA16-5":'sarah_16bba_4_A_data/bA16-5/parameters/prmtrs.json',
                "bB16-5": 'sarah_16bba_4_A_data/bB16-5/parameters/prmtrs.json',
                "B16-5": 'sarah_16bba_4_A_data/bigb-5/parameters/prmtrs.json',
                 #cam5:
                "bA16-8":'sarah_16bba_5_A_data/bA16-8/parameters/prmtrs.json',
                "Bb16-8":'sarah_16bba_5_A_data/Bb16-8/parameters/prmtrs.json',
                "B16-8":'sarah_16bba_5_A_data/B16-8/parameters/prmtrs.json',
                "b16-8": 'sarah_16bba_5_A_data/smallb16-8/parameters/prmtrs.json',
                #cam6 
                 "bA16-7":'sarah_16bba_6_A_data/bA16-7/parameters/prmtrs.json',
                 "B16-7":'sarah_16bba_6_A_data/bigB16-7/parameters/prmtrs.json',
                 "Bb16-7":'sarah_16bba_6_A_data/Bb16-7/parameters/prmtrs.json',
                 "b16-7": 'sarah_16bba_6_A_data/b16-7/parameters/prmtrs.json',
                 "bA16-6": 'sarah_16bba_6_A_data/bA16-6/parameters/prmtrs.json',
                 "bB16-6": 'sarah_16bba_6_A_data/bB16-6/parameters/prmtrs.json',
                 
                  #cam2:
                  "Bb16-2":'sarah_16bba_2_A_data/Bb16-2/parameters/prmtrs.json',
                  "B16-2":'sarah_16bba_2_A_data/bigb16-2/parameters/prmtrs.json',
                  "b16-3":'sarah_16bba_2_A_data/b16-3/parameters/prmtrs.json',
                  "bA16-2": 'sarah_16bba_2_A_data/bA16-2/parameters/prmtrs.json',
                  "bB16-3": 'sarah_16bba_2_A_data/bB16-3/parameters/prmtrs.json',
                  "B16-3": 'sarah_16bba_2_A_data/bigb16-3/parameters/prmtrs.json'
                    
                  
   }


########
#%% pick colony 
colony_name="B16-7"

bg_path=bg_path_list[colony_name] 
mask_path = mask_path_list[colony_name]
json_path = json_path_list[colony_name]
#bg_path = 'sarah_3216ab_5_inf_20240626_161922/ba16-6/parameters/backgrounds/background.png'
#mask_path = 'sarah_3216ab_5_inf_20240626_161922/ba16-6/parameters/masks/'
#json_path = 'sarah_3216ab_5_inf_20240626_161922/ba16-6/parameters/prmtrs.json'


#########
#%%parameters
display_scale = 0.5  # e.g., 0.3–0.7 depending on your screen

point = None
dragging = False
radius = 6        # how close the mouse must be to "grab" the point

#visible radius for exp2: cam2: 350:
#exp1 cam 5 inf a16-6: 325
#exp1 cam 5 inf ba16-6: 328
     
scaled_radius=345
visible_radius=int(round(scaled_radius*display_scale,0))

########################
#FUNCTIONS
#%%rescale background image to see all colonies: updated!


def mouse_callback(event, x, y, flags, param):
    
    global point, dragging

    if event == cv2.EVENT_LBUTTONDOWN:
        # If point exists, check if click is near it
        if point is not None:
            px, py = point
            if (x - px)**2 + (y - py)**2 <= radius**2:
                dragging = True
        else:
            # Create new point
            point = (x, y)

    elif event == cv2.EVENT_MOUSEMOVE and dragging:
        # Move point with mouse
        point = (x, y)

    elif event == cv2.EVENT_LBUTTONUP:
        dragging = False


def get_colonyarea(bg_path, colony_name, display_scale):
    background = cv2.imread(bg_path)
    global point, dragging
    point = None
    dragging = False

    
    img_unscaled = background.copy()
    img = cv2.resize(img_unscaled, None, fx=display_scale, fy=display_scale)

    
    cv2.namedWindow(colony_name)
    cv2.setMouseCallback(colony_name, mouse_callback)

    while True:
        frame = img.copy()

        # Draw draggable point
        if point is not None:
            cv2.circle(frame, point, visible_radius, (0, 255, 0), -1)

        cv2.imshow(colony_name, frame)

        key = cv2.waitKey(20)
        if key == 13:  # Enter to confirm
            break
        if key == 27:  # ESC to cancel
            point = None
            break

    cv2.destroyWindow(colony_name)
    # --- Scale points back to original size ---
    point_scaled = [int(point[0] / display_scale), int(point[1] / display_scale)]

    return point_scaled
################
#%%
point_scaled=get_colonyarea(bg_path, colony_name, display_scale)

################
#############
#%% get vertices and get rscale and a and b 

def get_vertices(center, visible_radius, n_vertices):
    cx, cy = center
    
    #angles to equalitistantly cover area of circle
    theta = np.linspace(0, 2*np.pi, n_vertices, endpoint=False)

    x = cx + visible_radius * np.cos(theta)
    y = cy + visible_radius * np.sin(theta)

    # px/mm
    rscale=(2*visible_radius)/45
    #rscale for antrax
    r_scale_antrax= (1/rscale)*10**(-3)
    
    return (np.vstack((x, y)).T, r_scale_antrax)
############
#%%
vertices, r_scale = get_vertices(point_scaled, scaled_radius, 136)
########
#%%create roimask.png and replace in path
# openboundrymask.png automatically generated
# openboundryperimmask.png 
def get_masks(visible_radius, center, mask_path, bg_path):
    
    background = cv2.imread(bg_path)
    height, width = background.shape[:2]

    #black mask
    mask = np.zeros((height, width), dtype=np.uint8)
    
    cv2.circle(mask, center, visible_radius,255,-1)
   
    #save in antrax directory /parameters/mask
    
    combined_path = os.path.join(mask_path,
                                 "roimask.png")
    cv2.imwrite(combined_path, mask)
    
    


###########
#%%
get_masks(scaled_radius, point_scaled, mask_path, bg_path)

###########
#%% save parameters in prmtrs.json
def update_prmtrs_jsn(json_path, center, visible_radius, vertices, r_scale, out_path=None):
    
    
    with open(json_path, "r") as f:
        data = json.load(f)
    #replace rscale
    data["geometry_rscale"] = r_scale
    data["geometry_scale_tool_meas"] = 45
    
    #add colony area
    geom = data["geometry_scale_tool_params"]
    geom["type"]= "ellipse"
    #create new substructures
    geom["SemiAxes"]=[visible_radius , visible_radius]
    geom["Center"]= [center[0], center[1]]
    geom["RotationAngle"] = 0
    geom["Vertices"]= vertices.tolist()
    
    #overwrite original
    # write back
    if out_path is None:
        out_path = json_path   # overwrite original

    with open(out_path, "w") as f:
        json.dump(data, f, indent=4)

    print("Updated JSON saved")
    
##############
#%%
point_scaled=get_colonyarea(bg_path, colony_name, display_scale)
vertices, r_scale = get_vertices(point_scaled, scaled_radius, 136)
get_masks(scaled_radius, point_scaled, mask_path, bg_path)
update_prmtrs_jsn(json_path, point_scaled, scaled_radius, vertices, r_scale)