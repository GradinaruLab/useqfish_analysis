import zarr
import napari
import argparse

import os
import glob

from image_manipulation import *
from params import *

my_parser = argparse.ArgumentParser(description='Run analysis on a group of images')

my_parser.add_argument('Path',
                       metavar='path',
                       type=str,
                       help='the path to the diretory containing the images')


# Parse arguments
args = my_parser.parse_args()

path = args.Path
filepath = os.path.join(path, '*.tif')    
filenames = glob.glob(filepath)

imgReference = image_read(filenames[roundRef])

imgs = zarr.load(os.path.join(path, 'result_images.zip'))
zToXYRatioReal = imgs['zToXYRatioReal']
cellLabels = imgs['cellLabels']
cellOutlines = imgs['cellOutlines']
spots_assigned_allrounds = imgs['spots_assigned_allrounds']

print(f'>> STEP 5. visualization')
zyxScale = (zToXYRatioReal, 1, 1)
with napari.gui_qt():
    viewer = napari.Viewer()
    # for img in results:
    #     viewer.add_image(img)

    viewer.add_image(
        imgReference[0],
        name='dapi reference',
        scale=zyxScale,
        contrast_limits=[0, imgReference[0].max()],
        multiscale=False
    )
    viewer.add_image(
        imgReference[cellch],
        name='cell reference',
        scale=zyxScale,
        contrast_limits=[0, imgReference[cellch].max()],
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
                spots = spots[:,:-1]
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