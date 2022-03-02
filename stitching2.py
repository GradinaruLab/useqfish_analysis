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
import zarr
import scanpy as sc

from warnings import filterwarnings; filterwarnings("ignore")

my_parser = ArgumentParser(description='Run analysis on a group of images')

my_parser.add_argument('Path',
                       metavar='path',
                       type=str,
                       help='the path to the diretory containing the images')

my_parser.add_argument('--mip', action='store_true')
my_parser.add_argument('--cell', action='store_true')

# Parse arguments
args = my_parser.parse_args()

path = args.Path
ismip = args.mip
iscell = args.cell

nR = 5
cmap = get_cmap('rainbow', nR*nC)
# colors = [matplotlib.colors.rgb2hex(c) for c in cmap.colors]
colors = [rgb2hex(cmap(i)) for i in range(cmap.N)]
resolution = 0.151
# stitching_coords = np.array(stitching_coords)
stitching_coords = np.multiply(np.array(stitching_coords), resolution)
stitching_coords = np.array(stitching_coords) - np.array(stitching_coords[0])
stitching_coords = list(stitching_coords)
montage_w, montage_h = stitching_shape

no_overlapped_size = np.array(img_size)-np.array(img_size)*stitching_overlap/100


if ismip:
    filepath = os.path.join(path, 'result/*.xlsx')
    filenames = sorted(glob(filepath), key=os.path.basename)
    n_files = len(filenames)

    stitching_composition = np.zeros((montage_w, montage_h))

    img_index = 0
    for w in range(montage_w):
        if w%2 == 0:
            for h in range(montage_h):
                stitching_composition[w, h] = img_index
                img_index = img_index+1
        else:
            for h in range(montage_h-1, -1, -1):
                stitching_composition[w, h] = img_index
                img_index = img_index+1

    # # print(stitching_composition)
    # # print(stitching_composition[2, 3])
    # # print(np.argwhere(stitching_composition==15).ravel())
    
    # file_index = 0
    # spots = []

    # global_grid = []
    # spots = []
    # for file_index, file in enumerate(filenames):
    #     tile_w, tile_h = np.argwhere(stitching_composition==file_index).ravel()

    #     global_w = tile_w * no_overlapped_size[0] + stitching_coords[file_index][0]
    #     global_h = -(tile_h * no_overlapped_size[1] + stitching_coords[file_index][1])

    #     # global_grid.append([global_w, global_h])

    #     spot = pd.read_excel(filenames[file_index], index_col=0)
    #     spot_notoverlapped = spot
    #     if (tile_w == montage_w-1) & (tile_h < montage_h-1):
    #         spot_notoverlapped = spot.loc[spot['y-coord']<no_overlapped_size[1]]
    #     if (tile_w < montage_w-1) & (tile_h == montage_h-1):
    #         spot_notoverlapped = spot.loc[spot['x-coord']<no_overlapped_size[0]]
    #     if (tile_w < montage_w-1) & (tile_h < montage_h-1):
    #         spot_notoverlapped = spot.loc[(spot['x-coord']<no_overlapped_size[0]) & (spot['y-coord']<no_overlapped_size[1])]

    #     spot_notoverlapped['x-coord'] = spot_notoverlapped['x-coord'] + global_w
    #     spot_notoverlapped['y-coord'] = -(spot_notoverlapped['y-coord'] + global_h)

    #     spots.append(spot_notoverlapped)

    # spots = pd.concat(spots, axis=0, ignore_index=True)

    # spots_assigned_allrounds = []
    # for r in range(nR):
    #     spots_assigned = []
    #     for c in range(nC):
    #         coords = spots.loc[(spots['round']==r+1) & (spots['channel']==c+1), ['y-coord', 'x-coord']].to_numpy()
    #         spots_assigned.append(coords)
    #     spots_assigned_allrounds.append(spots_assigned)


    if iscell:
        adata_endo = sc.read_h5ad(os.path.join(path, 'h5ad_clustered/endo.h5ad'))
        cell_filter = adata_endo.uns['cell_subset']
        cell_labels_stitched = np.zeros(stitching_size, dtype=np.uint)

        cell_index = 1
        label = 1
        zarr_folders = sorted([f.path for f in os.scandir(os.path.join(path, 'cell_labels')) if f.is_dir()])
        for file_index, zarr_folder in enumerate(zarr_folders):
            cell_labels = zarr.load(zarr_folder)
            n_cells = cell_labels.max()
            
            tile_w, tile_h = np.argwhere(stitching_composition==file_index).ravel()

            global_w = tile_w * no_overlapped_size[0] + stitching_coords[file_index][0]
            global_h = -(tile_h * no_overlapped_size[1] + stitching_coords[file_index][1])

            for cell in range(n_cells):
                if cell_filter[cell_index]:
                    cell_coords = np.argwhere(cell_labels==cell+1)[:,:-1]
                    
                    

                    label = label+1
                cell_index = cell_index+1
                    
                    
                    exit()



        

    viewer = napari.Viewer()
    cind = 0
    for r, spots_assigned in enumerate(spots_assigned_allrounds):
        for c, spots in enumerate(spots_assigned):
            spots = np.array(spots)
            # color = colors[nC*r+c]
            # color = sample(colors,1)
            color = colors[cind]
            # print(f'cind: {cind}')
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
                cind = cind + 1
                  

    napari.run()

    
    
