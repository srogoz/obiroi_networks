#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 18 12:57:48 2025

@author: srogoz@ice.mpg.de
"""
##############
#%% nest location using background to larvae heap threahold on greyscale
import cv2
import matplotlib.pyplot as plt
import numpy as np
import os

os.chdir('/home/srogoz@ice.mpg.de/ulr-lab/Users/Sarah_saro7615/EXP1_GROUPSIZEandCOMPOSITION/')
#force to convert to one channel
image_path = 'sarah_3216ab_2_0_lasthalf/b16-2/parameters/backgrounds/background.png'
roi = cv2.imread( 'sarah_3216ab_2_0_lasthalf/b16-2/parameters/masks/colony_C1.png' , cv2.IMREAD_GRAYSCALE)
img = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)


_, roi_binary = cv2.threshold(roi, 127, 255, cv2.THRESH_BINARY)

masked_img = cv2.bitwise_and(img, img, mask=roi_binary)

# Apply threshold to the grayscale image
_, thresholded = cv2.threshold(img, 205, 255, cv2.THRESH_BINARY)

# Mask out everything outside ROI
thresholded_roi = cv2.bitwise_and(thresholded, thresholded, mask=roi_binary)

# Show results
cv2.imshow("Thresholded Inside ROI", thresholded_roi)
cv2.waitKey(0)
cv2.destroyAllWindows()


############
#%%
#cv2.imwrite('nests/cam_2/masks/mask-a16-2.png', thresholded_roi)´
img5 = cv2.imread('nests/cam_2/masksmask-a16-2.png')
############
#%%
import numpy as np
import os
import cv2
import matplotlib.pyplot as plt
import pickle
os.chdir('/home/srogoz@ice.mpg.de/ulr-lab/Users/Sarah_saro7615/EXP1_GROUPSIZEandCOMPOSITION/')

#######
#%% DATA SOURCING
cam1 = [ "nests/cam_1/cam1_tp1.png",
         "nests/cam_1/cam1_tp2.png",
         "nests/cam_1/cam1_tp3.png",
        "nests/cam_1/cam1_tp4.png"]

cam2 = [ "nests/cam_2/cam2_tp1.png",
         "nests/cam_2/cam2_tp2.png",
         "nests/cam_2/cam2_tp3.png",
        "nests/cam_2/cam2_tp4.png"]

cam3 = [ "nests/cam_3/cam3_tp1.png",
         "nests/cam_3/cam3_tp2.png",
         "nests/cam_3/cam3_tp3.png",
        "nests/cam_3/cam3_tp4.png"]

cam4 = [ "nests/cam_4/cam4_tp1.png",
         "nests/cam_4/cam4_tp2.png",
         "nests/cam_4/cam4_tp3.png",
        "nests/cam_4/cam4_tp4.png"]

cam5 = [ "nests/cam_5/cam5_tp1.png",
         "nests/cam_5/cam5_tp2.png",
         "nests/cam_5/cam5_tp3.png",
        "nests/cam_5/cam5_tp4.png"]

cam6 = [ "nests/cam_6/cam6_tp1.png",
         "nests/cam_6/cam6_tp2.png",
         "nests/cam_6/cam6_tp3.png",
        "nests/cam_6/cam6_tp4.png"]
#########
#%%

def get_nestoutlines(cam, colony_name):
    colony_path = f"nests/{cam}/{colony_name}"
    os.mkdir(colony_path)
    all_points = []
    for i, img_path in enumerate(cam):
    


        img = cv2.imread(img_path)
        
        ########################
        drawing = False  # True if mouse is pressed
        points = []  # Store polygon points
        mask_name = f"mask_tp_{i}"
        
        def draw_polygon(event, x, y, flags, param):
            global drawing, points, img_copy
        
            if event == cv2.EVENT_LBUTTONDOWN:
                points.append((x, y))
                drawing = True
                # Draw a small circle at clicked point
                cv2.circle(img_copy, (x, y), 3, (0, 255, 0), -1)
        
            elif event == cv2.EVENT_MOUSEMOVE and drawing:
                # Optional: show line from last point to current mouse position
                img_copy = img.copy()
                if len(points) > 0:
                    cv2.polylines(img_copy, [np.array(points + [(x, y)])], False, (0, 255, 0), 2)
                for p in points:
                    cv2.circle(img_copy, p, 3, (0, 255, 0), -1)
        
            elif event == cv2.EVENT_LBUTTONUP:
                drawing = False
        
        # Copy to draw on
        img_copy = img.copy()
        cv2.namedWindow('Draw ROI')
        cv2.setMouseCallback('Draw ROI', draw_polygon)
        
        while True:
            cv2.imshow('Draw ROI', img_copy)
            key = cv2.waitKey(1) & 0xFF
        
            # Press 'c' to clear points
            if key == ord('c'):
                points = []
                img_copy = img.copy()
        
            # Press 'q' to quit drawing and create mask
            elif key == ord('q'):
                break
        
        cv2.destroyAllWindows()
        
        # After quitting, create a mask from the polygon
        mask = np.zeros(img.shape[:2], dtype=np.uint8)
        if len(points) > 2:
            cv2.fillPoly(mask, [np.array(points)], 255)
        
        #save Mask
         
        cv2.imwrite(f"nests/{cam}/{colony_name}/masks/mask_{i}.png", mask)
        
        # Show mask
        cv2.imshow("ROI Mask", mask)
        cv2.waitKey(0)
        cv2.destroyAllWindows()

##q after all points are added
#c if want to remove all points in GUI
#ENTER to close mask display, otherwise GUI interferes with further use!


#cv2.imwrite("/nests/cam_2/masks/a16-2_1.png", mask)
# collect outlines in all_points, with option for saving
# .pkl and .npy both work, not displayed in folder, but seen on terminal


        all_points.append(points)

#with open('nests/cam_2/all_points.pkl', 'wb') as f:
    #pickle.dump(all_points, f)
    
        np.save('nests/{cam}/{colony_name}/all_points.npy', np.array(all_points, dtype=object))  # Save
##################
#%% display all outlines on img 1 
colors = [(162, 0, 255), (204, 129, 249), (188, 129, 247), (129, 247, 247)]

for idx, outline in enumerate(all_points):
    if len(outline) < 2:  # Skip if too few points
        continue
    pts = np.array(outline, np.int32).reshape((-1, 1, 2))
    color = colors[idx % len(colors)]  # Cycle through colors if more outlines
    cv2.polylines(img, [pts], isClosed=True, color=color, thickness=2)

################
#%%
# save outlines 
#
cv2.imwrite(f'nests/{cam}/{colony_name}/combined_outlines_{colony_name}.png', img)

####
#%% shape sanity check

##############################################
#%%not working


def get_nestoutlines(cam_list, colony_name):
    """
    Interactive ROI drawing tool for a set of images.
    Args:
        cam_list (list[str]): list of image paths
        colony_name (str): name for saving folder and files
    """
    # --- Setup ---
    cam_dir = os.path.dirname(cam_list[0])  # e.g. nests/cam_3
    save_dir = os.path.join(cam_dir, colony_name, "masks")
    os.makedirs(save_dir, exist_ok=True)
    all_points = []

    # --- Iterate over each timepoint image ---
    for i, img_path in enumerate(cam_list):
        img = cv2.imread(img_path)
        display_scale = 0.5  # e.g., 0.3–0.7 depending on your screen
        img_display = cv2.resize(img, None, fx=display_scale, fy=display_scale)
        
                
        drawing = False
        points = []
        

        def draw_polygon(event, x, y, flags, param):
            nonlocal drawing, points, img_copy
            if event == cv2.EVENT_LBUTTONDOWN:
                points.append((x, y))
                drawing = True
                cv2.circle(img_copy, (x, y), 3, (0, 255, 0), -1)

            elif event == cv2.EVENT_MOUSEMOVE and drawing:
                img_copy = img.copy()
                if len(points) > 0:
                    cv2.polylines(img_copy, [np.array(points + [(x, y)])],
                                  False, (0, 255, 0), 2)
                for p in points:
                    cv2.circle(img_copy, p, 3, (0, 255, 0), -1)

            elif event == cv2.EVENT_LBUTTONUP:
                drawing = False

        # --- Draw ROI interactively ---
        img_copy = img.copy()
        cv2.namedWindow('Draw ROI')
        cv2.setMouseCallback('Draw ROI', draw_polygon)

        
        print("   Left-click = add points | c = clear | q = confirm and save\n")

        while True:
            cv2.imshow('Draw ROI', img_copy)
            key = cv2.waitKey(1) & 0xFF
            if key == ord('c'):
                points = []
                img_copy = img.copy()
            elif key == ord('q'):
                break

        cv2.destroyAllWindows()
        
        points_scaled = [(int(x / display_scale), int(y / display_scale)) for (x, y) in points]
        # --- Create and save mask ---
        mask = np.zeros(img.shape[:2], dtype=np.uint8)
        if len(points_scaled) > 2:
            cv2.fillPoly(mask, [np.array(points_scaled)], 255)
        


        mask_path = os.path.join(save_dir, f"mask_tp_{i+1}.png")
        cv2.imwrite(mask_path, mask)
      

        # --- Show mask briefly ---
        cv2.imshow("ROI Mask", mask)
        cv2.waitKey(0)
        cv2.destroyAllWindows()

        # --- Save polygon points ---
        all_points.append(points)

    # --- Save all outlines ---
    np.save(os.path.join(save_dir, "all_points.npy"), np.array(all_points, dtype=object))
    with open(os.path.join(save_dir, "all_points.pkl"), 'wb') as f:
        pickle.dump(all_points, f)
   

    # --- Draw all outlines on first image ---
    colors = [(162, 0, 255), (204, 129, 249), (188, 129, 247), (129, 247, 247)]
    ###########
    for outline in all_points:
    pts = np.array(
        [(int(x / display_scale), int(y / display_scale)) for (x, y) in outline],
        np.int32
    ).reshape((-1, 1, 2))
    cv2.polylines(img, [pts], isClosed=True, color=color, thickness=2)
    
    #########
    overlay = cv2.imread(cam_list[0]).copy()
    for idx, outline in enumerate(all_points):
        if len(outline) < 2:
            continue
        pts = np.array(outline, np.int32).reshape((-1, 1, 2))
        color = colors[idx % len(colors)]
        cv2.polylines(overlay, [pts], isClosed=True, color=color, thickness=2)

    combined_path = os.path.join(cam_dir, colony_name,
                                 f"combined_outlines_{colony_name}.png")
    cv2.imwrite(combined_path, overlay)
    print(f" Saved combined outlines: {combined_path}")
    cv2.imshow("Combined Outlines", overlay)
    cv2.waitKey(0)
    cv2.destroyAllWindows()

    print("done")


################
#%%WORKING FUNCTION TO EXTRACT 

def get_nestoutlines(cam_list, colony_name):
    """
    Interactive ROI drawing tool for a set of images.
    Args:
        cam_list (list[str]): list of image paths
        colony_name (str): name for saving folder and files
    """
    # --- Setup ---
    cam_dir = os.path.dirname(cam_list[0])  # e.g. nests/cam_3
    save_dir = os.path.join(cam_dir, colony_name, "masks")
    os.makedirs(save_dir, exist_ok=True)
    all_points = []

    # --- Iterate over each timepoint image ---
    for i, img_path in enumerate(cam_list):
        img = cv2.imread(img_path)
        if img is None:
            print(f" no image at {img_path}")
            continue

        display_scale = 0.5  # e.g., 0.3–0.7 depending on your screen
        img_display = cv2.resize(img, None, fx=display_scale, fy=display_scale)

        drawing = False
        points = []

        def draw_polygon(event, x, y, flags, param):
            nonlocal drawing, points, img_display_copy
            if event == cv2.EVENT_LBUTTONDOWN:
                points.append((x, y))
                drawing = True
                cv2.circle(img_display_copy, (x, y), 3, (0, 255, 0), -1)

            elif event == cv2.EVENT_MOUSEMOVE and drawing:
                img_display_copy = img_display.copy()
                if len(points) > 0:
                    cv2.polylines(img_display_copy, [np.array(points + [(x, y)])],
                                  False, (0, 255, 0), 2)
                for p in points:
                    cv2.circle(img_display_copy, p, 3, (0, 255, 0), -1)

            elif event == cv2.EVENT_LBUTTONUP:
                drawing = False

        # --- Draw ROI interactively ---
        img_display_copy = img_display.copy()
        cv2.namedWindow('Draw ROI')
        cv2.setMouseCallback('Draw ROI', draw_polygon)

      
        print("   Left-click = add points | c = clear | q = confirm and save\n")

        while True:
            cv2.imshow('Draw ROI', img_display_copy)
            key = cv2.waitKey(1) & 0xFF
            if key == ord('c'):
                points = []
                img_display_copy = img_display.copy()
            elif key == ord('q'):
                break

        cv2.destroyAllWindows()

        # --- Scale points back to original size ---
        points_scaled = [(int(x / display_scale), int(y / display_scale)) for (x, y) in points]

        # --- Create and save mask ---
        mask = np.zeros(img.shape[:2], dtype=np.uint8)
        if len(points_scaled) > 2:
            cv2.fillPoly(mask, [np.array(points_scaled)], 255)

        mask_path = os.path.join(save_dir, f"mask_tp_{i+1}.png")
        cv2.imwrite(mask_path, mask)
       

        # --- Show mask briefly ---
        cv2.imshow("ROI Mask", mask)
        cv2.waitKey(0)
        cv2.destroyAllWindows()

        # --- Save polygon points ---
        all_points.append(points_scaled)

    # --- Save all outlines ---
    np.save(os.path.join(save_dir, "all_points.npy"), np.array(all_points, dtype=object))
    with open(os.path.join(save_dir, "all_points.pkl"), 'wb') as f:
        pickle.dump(all_points, f)
   

    # --- Draw all outlines on first image ---
    overlay = cv2.imread(cam_list[0]).copy()
    colors = [(162, 0, 255), (204, 129, 249), (188, 129, 247), (129, 247, 247)]

    for idx, outline in enumerate(all_points):
        if len(outline) < 2:
            continue
        pts = np.array(outline, np.int32).reshape((-1, 1, 2))
        color = colors[idx % len(colors)]
        cv2.polylines(overlay, [pts], isClosed=True, color=color, thickness=2)

    combined_path = os.path.join(cam_dir, colony_name,
                                 f"combined_outlines_{colony_name}.png")
    cv2.imwrite(combined_path, overlay)
   

    cv2.imshow("Combined Outlines", overlay)
    cv2.waitKey(0)
    cv2.destroyAllWindows()

    print("done")




################
#%%
get_nestoutlines(cam6, "a16-8")