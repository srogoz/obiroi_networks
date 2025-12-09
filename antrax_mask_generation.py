#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 13:19:06 2025

@author: srogoz@ice.mpg.de
"""

import cv2
import matplotlib.pyplot as plt
import numpy as np
import os
############
## get files with already done backgrounds
os.chdir('/home/srogoz@ice.mpg.de/ulr-lab/Users/Sarah_saro7615/EXP2/')
#force to convert to one channel
bg_path = 'sarah_16bba_2_A_data/b16-3/parameters/backgrounds/background.png'
mask_path = 'sarah_16bba_2_A_data/b16-3/parameters/masks/'
background = cv2.imread(bg_path)

############paint elipse
#%%
n_vertices= 136

def get_colonyarea(background, colony_name):
    drawing = False  # True if mouse is pressed
    center = []  # Store center point
    
    
    if event == cv2.EVENT_LBUTTONDOWN:
        center.append((x, y))
        drawing = True
        # Draw a small circle at clicked point
        cv2.circle(background, (x, y), 3, (0, 255, 0), -1)
        
    ##store mask 
    #mask_name = f"{colony_name}/masks/roimask.png"
    return center)
    
    
####################
#%%
    
import cv2
import numpy as np
colony_name="b16-3"
# GLOBALS for simplicity
point = None
dragging = False
radius = 6        # how close the mouse must be to "grab" the point
##
#visible radius for exp2: cam2: 350:
    #rscale: 
visible_radius=350

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


def get_colonyarea(background, colony_name):
    global point, dragging
    point = None
    dragging = False

    img = background.copy()
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
    return point

################
#%%
get_colonyarea(background, colony_name)
################
#%%rescale image to see all colonies
#display_scale = 0.5  # e.g., 0.3–0.7 depending on your screen
#img_display = cv2.resize(img, None, fx=display_scale, fy=display_scale)


point = None
dragging = False
radius = 6        # how close the mouse must be to "grab" the point
##
#visible radius for exp2: cam2: 350:
    #rscale: 
visible_radius=350

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


def get_colonyarea(background, colony_name):
    global point, dragging
    point = None
    dragging = False


    img = background.copy()
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
    return point
########################
#%%rescale background image to see all colonies: updated!



point = None
dragging = False
radius = 6        # how close the mouse must be to "grab" the point
##
#visible radius for exp2: cam2: 350:
    #rscale: 
visible_radius=int(round(350*0.5,0))

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


def get_colonyarea(background, colony_name):
    global point, dragging
    point = None
    dragging = False

    display_scale = 0.5  # e.g., 0.3–0.7 depending on your screen
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
vertices, r_scale = get_vertices(point, visible_radius, 136)
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
get_masks(visible_radius,point, mask_path, bg_path)

###########
#%% save parameters in prmtrs.json

import json

json_path = 'sarah_16bba_2_A_data/b16-3/parameters/prmtrs.json'

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
update_prmtrs_jsn(json_path, point, visible_radius, vertices, r_scale)