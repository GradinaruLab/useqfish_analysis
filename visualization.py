import zarr
import argparse

import os

from params import *

import pandas as pd

import napari
from pylab import *

my_parser = argparse.ArgumentParser(description='Run analysis on a group of images')

my_parser.add_argument('Path',
                       metavar='path',
                       type=str,
                       help='the path to the diretory containing the images')


# Parse arguments
args = my_parser.parse_args()

path = args.Path

cells = zarr.load(os.path.join(path, 'result_cellsegments.zarr'))
zToXYRatioReal = cells['zToXYRatioReal']
imgCells = cells['imgCells']
cellLabels = cells['cellLabels']
cellOutlines = cells['cellOutlines']
nR = cells['nR']
spots = pd.read_excel(os.path.join(path, 'result.xlsx'), index_col=0)

spots_assigned_allrounds = []
for r in range(nR):
    spots_assigned = []
    for c in range(nC):
        coords = spots.loc[(spots['round']==r) & (spots['channel']==c), ['z-coord', 'y-coord', 'x-coord']].to_numpy()
        spots_assigned.append(coords)
    spots_assigned_allrounds.append(spots_assigned)

zyxScale = (zToXYRatioReal, 1, 1)
# with napari.gui_qt():

viewer = napari.Viewer()
viewer.add_image(
    imgCells[:,:,:,1],
    name='dapi reference',
    scale=zyxScale,
    contrast_limits=[0, imgCells[:,:,:,1].max()],
    multiscale=False
)
viewer.add_image(
    imgCells[:,:,:,0],
    name='cell reference',
    scale=zyxScale,
    contrast_limits=[0, imgCells[:,:,:,0].max()],
    multiscale=False
)
viewer.add_labels(
    cellLabels,
    name='cell labels',
    scale=zyxScale,
    multiscale=False
)
viewer.add_image(
    cellOutlines,
    name='cell outlines',
    scale=zyxScale,
    multiscale=False
)

# cmap = np.linspace(0, 1, num=(nR-1)*(nC-1))
cmap = cm.get_cmap('turbo', (nR-1)*(nC-1))
colors = [matplotlib.colors.rgb2hex(c) for c in cmap.colors]
cIndex = 0
for r, spots_assigned in enumerate(spots_assigned_allrounds):
    for c, spots in enumerate(spots_assigned):
        spots = np.array(spots)
        if spots.shape[0] > 0:
            # print(spots.shape)
            # pointProperties = {
            #     'good_point': np.ones((spots.shape[0],), dtype=bool),
            #     'confidence': np.full((spots.shape[0],), cmap[cIndex])
            # }
            # print(cmap[cIndex])
            viewer.add_points(
                spots,
                face_color=colors[cIndex],
                size=10,
                n_dimensional=True,
                name=f'{r} round, {c+1} ch',
                scale=zyxScale,
                blending='additive'
            )
            cIndex = cIndex + 1
napari.run()

    