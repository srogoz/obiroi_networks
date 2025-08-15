#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 31 10:41:58 2025

@author: srogoz@ice.mpg.de
"""

#########################
# %% create heatmaps of ants movement

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os


########drop colony index, if single colony pipeline
colony = combined_colony_data_interpolated.reset_index(level = "colony", drop = True)
###############
# %%
import os
import numpy as np
import matplotlib.pyplot as plt

export_dir = 'processed/plots/b16-1/heatmaps'
os.makedirs(export_dir, exist_ok=True)


def plot_heatmaps(df):
    """
    Generate and save heatmaps of x/y positions for each ant.
    Assumes MultiIndex: ['ant', 'frame']
    """
    all_ants = df.index.get_level_values('ant').unique()
    export_dir = 'processed/plots/b16-1/heatmaps'
    os.makedirs(export_dir, exist_ok=True)

    for ant in all_ants:
        
        ant_df = df.loc[ant]
        ant_df = ant_df.dropna(subset=['x', 'y'])

        if ant_df.empty:
            continue

        x = ant_df['x'].values
        y = ant_df['y'].values

        x_min, x_max = x.min(), x.max()
        y_min, y_max = y.min(), y.max()
        bins = 90

        heatmap, xedges, yedges = np.histogram2d(
            x, y, bins=bins, range=[[x_min, x_max], [y_min, y_max]]
        )

        plt.figure(figsize=(8, 6))
        plt.imshow(
            heatmap.T, origin='lower', extent=[x_min, x_max, y_min, y_max],
            cmap='viridis', aspect='auto'
        )
        plt.colorbar(label='Density [frames]')
        plt.xlabel('x [mm]')
        plt.ylabel('y [mm]')
        plt.title(f'colony b16-1, ant: {ant} — time spent [frames], 10 frames = 1s')
        plt.savefig(f'{export_dir}/{ant}.png')
        plt.close()
        print (f'{ant} done :)')
####################
#%% run 
plot_heatmaps(colony)
##############
# %%
