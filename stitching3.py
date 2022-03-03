
import scanpy as sc

import numpy as np

from glob import glob
import os
from argparse import ArgumentParser

from params import *

from converting_anndata import *

from warnings import filterwarnings; filterwarnings("ignore")

sc.settings.verbosity = 3

import napari
import zarr


my_parser = ArgumentParser(description="Run stitching for tiled image")

my_parser.add_argument('Path',
                       metavar='path',
                       type=str,
                       help='the path to the directory containing spot detection result (xlsx files)')

args = my_parser.parse_args()
path = args.Path

# if len(os.path.join(path, 'cell_labels_stitched.zarr'))==0:
    
filepath = os.path.join(path, 'expression_matrix.h5ad')

if len(filepath)==0:
    xlsx2h5ad(path)

adata = sc.read_h5ad(filepath)


# ## stitching
# nR = 5
# resolution = 0.151
# stitching_coords = np.multiply(np.array(stitching_coords), resolution)
# stitching_coords = list(np.array(stitching_coords) - np.array(stitching_coords[0]))
# montage_h, montage_w = stitching_shape

# stitching_composition = np.zeros(stitching_shape)
# img_index = 0
# for w in range(montage_w):
#     if w%2 == 0:
#         for h in range(montage_h):
#             stitching_composition[h, w] = img_index
#             img_index = img_index+1
#     else:
#         for h in range(montage_h-1, -1, -1):
#             stitching_composition[h, w] = img_index
#             img_index = img_index+1


# no_overlapped_size = np.array(img_size) - np.array(img_size)*stitching_overlap/100

# ## cell label stitching
# cell_labels_stitched = np.zeros(stitching_size, dtype=np.uint)
# cell_overlapped = []
# cell_locations = []
# label = 1
# zarr_folders = sorted([f.path for f in os.scandir(os.path.join(path, 'cell_labels')) if f.is_dir()])
# for tile_index, zarr_folder in enumerate(zarr_folders):
#     print(zarr_folder)
#     cell_labels = zarr.load(zarr_folder)
#     cell_labels = np.flip(cell_labels, axis=1)
#     n_cells = cell_labels.max()

#     tile_h, tile_w = np.argwhere(stitching_composition==tile_index).ravel()

#     global_h = tile_h * no_overlapped_size[0] + stitching_coords[tile_index][0]
#     global_w = tile_w * no_overlapped_size[1] + stitching_coords[tile_index][1]

#     for cell in range(n_cells):
#         # print(f'cell index: {cell+1}')
#         isoverlap = False
#         cell_coords = np.argwhere(cell_labels==cell+1)[:,1:]
#         orig_coords = cell_coords.copy()
#         if (tile_w == montage_w-1) & (tile_h < montage_h-1):
#             isoverlap = np.any(cell_coords[:,0] > no_overlapped_size[0])
#         if (tile_w < montage_w-1) & (tile_h == montage_h-1):
#             isoverlap = np.any(cell_coords[:,1] > no_overlapped_size[1])
#         if (tile_w < montage_w-1) & (tile_h < montage_h-1):
#             isoverlap = np.any(cell_coords[:,0] > no_overlapped_size[0]) | np.any(cell_coords[:,1] > no_overlapped_size[0])
#         cell_overlapped.append(isoverlap)
#         # print(f'overlapped? : {isoverlap}')
#         cell_location = [np.NaN, np.NaN]

#         if ~isoverlap:
#             cell_coords[:,0] = cell_coords[:,0] + global_h
#             cell_coords[:,1] = cell_coords[:,1] + global_w

#             subind = (cell_coords[:,0] >= 0) & (cell_coords[:,0] < stitching_size[0])
#             subind2 = (cell_coords[:,1] >= 0) & (cell_coords[:,1] < stitching_size[1])
#             subind = subind & subind2
#             cell_labels_stitched[cell_coords[subind,0], cell_coords[subind,1]] = label
#             cell_location = [np.mean(cell_coords, axis=0)]
#             label = label + 1

#         cell_locations.append(cell_location)
#         print(f'tile index: {tile_index}, cell index: {cell+1}, overlapped? {isoverlap}, label: {label-1}')
#         print(f'orginal coords: {np.mean(orig_coords, axis=0)}, stitched coords: {cell_location}')
# adata.obs['overlap_in_stitched'] = cell_overlapped
# adata.obs['cell_location'] = cell_locations
# adata.write_h5ad(os.path.join(path, 'expression_matrix_stitched.h5ad'))
# zarr.save(os.path.join(path, 'cell_labels_stitched.zarr'), cell_labels_stitched)

cell_labels_stitched = zarr.load(os.path.join(path, 'cell_labels_stitched.zarr'))
# cell_labels_stitched = np.random.random_integers(0, 10, stitching_size)
# cell_labels_stitched = zarr.load(os.path.join(path, 'cell_labels/position00.zarr'))
# cell_labels_stitched = np.flip(cell_labels_stitched, axis=1)
viewer = napari.Viewer()
viewer.add_labels(
    cell_labels_stitched,
    name='cell_labels',
    multiscale=False
)

napari.run()

