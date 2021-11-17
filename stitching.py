from zarr import load

from argparse import ArgumentParser
from glob import glob
import os
from random import sample

from params import *
from image_manipulation import *

import numpy as np
# from scipy.ndimage import shift
import pandas as pd

from matplotlib.cm import get_cmap
from matplotlib.colors import rgb2hex

import napari

from warnings import filterwarnings; filterwarnings("ignore")

my_parser = ArgumentParser(description='Run analysis on a group of images')

my_parser.add_argument('Path',
                       metavar='path',
                       type=str,
                       help='the path to the diretory containing the images')

my_parser.add_argument('--mip', action='store_true')

# Parse arguments
args = my_parser.parse_args()

path = args.Path
ismip = args.mip

cmap = get_cmap('hsv', nR*nC)
# colors = [matplotlib.colors.rgb2hex(c) for c in cmap.colors]
colors = [rgb2hex(cmap(i)) for i in range(cmap.N)]
print(colors)
# colors.reverse()

# stitching_coords = np.array(stitching_coords)
stitching_coords = np.array(stitching_coords) - np.array(stitching_coords[0])
stitching_coords = list(stitching_coords)
montage_w, montage_h = stitching_shape

if ismip:
    filepath = os.path.join(path, '*.xlsx')
    filenames = sorted(glob(filepath), key=os.path.basename)

    file_index = 0
    spots = []
    last_overlap = np.array(stitching_shape) * np.array(img_size) \
        - (np.array(stitching_shape)-1)*np.array(img_size)*stitching_overlap/100 \
        - stitching_size

    for w in range(montage_w):
        if w%2 == 0:
            for h in range(montage_h):
                spot = pd.read_excel(filenames[file_index], index_col=0)
                spot['x-coord'] = spot['x-coord'] + w*img_size[0] - w*img_size[0]*stitching_overlap/100
                spot['y-coord'] = -spot['y-coord'] + h*img_size[1] - h*img_size[1]*stitching_overlap/100
                if w==montage_w-1:
                    spot['x-coord'] = spot['x-coord'] - last_overlap[0]
                if h==montage_h-1:
                    spot['y-coord'] = spot['y-coord'] - last_overlap[1]
                spots.append(spot)
                file_index = file_index + 1
        else:
            for h in range(4, -1, -1):
                spot = pd.read_excel(filenames[file_index], index_col=0)
                spot['x-coord'] = spot['x-coord'] + w*img_size[0] - w*img_size[0]*stitching_overlap/100
                spot['y-coord'] = -spot['y-coord'] + h*img_size[1] - h*img_size[1]*stitching_overlap/100
                if w==montage_w-1:
                    spot['x-coord'] = spot['x-coord'] - last_overlap[0]
                if h==montage_h-1:
                    spot['y-coord'] = spot['y-coord'] - last_overlap[1]
                spots.append(spot)
                file_index = file_index + 1

    for spot, coord in zip(spots, stitching_coords):
        spot['x-coord'] = spot['x-coord'] + coord[0]
        spot['y-coord'] = spot['y-coord'] - coord[1]
    
    spots = pd.concat(spots, axis=0, ignore_index=True)

    spots_assigned_allrounds = []
    for r in range(nR):
        spots_assigned = []
        for c in range(nC):
            coords = spots.loc[(spots['round']==r+1) & (spots['channel']==c+1), ['y-coord', 'x-coord']].to_numpy()
            spots_assigned.append(coords)
        spots_assigned_allrounds.append(spots_assigned)


    viewer = napari.Viewer()
    for r, spots_assigned in enumerate(spots_assigned_allrounds):
        for c, spots in enumerate(spots_assigned):
            spots = np.array(spots)
            # color = colors[nC*r+c]
            color = sample(colors,1)
            if spots.shape[0] > 0:
                viewer.add_points(
                    spots,
                    face_color=color,
                    # edge_color=color,
                    # symbol='ring',
                    size=20,
                    name=f'{r+1} round, {c+1} ch',
                    # visible=False
                )

    napari.run()

    
    
