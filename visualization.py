import zarr
import argparse
import glob

import os
import random

from params import *

import pandas as pd

import napari
# from pylab import *
import matplotlib 

from cell_detection import *

from scipy import ndimage

import warnings
warnings.filterwarnings("ignore")

my_parser = argparse.ArgumentParser(description='Run analysis on a group of images')

my_parser.add_argument('Path',
                       metavar='path',
                       type=str,
                       help='the path to the diretory containing the images')

my_parser.add_argument('--mip', action='store_true')


# Parse arguments
args = my_parser.parse_args()

path = args.Path
ismip = args.mip


cells = zarr.load(os.path.join(path, 'result/result_images.zarr'))
zToXYRatioReal = cells['zToXYRatioReal']
imgCells = cells['imgCells']
cellLabels = cells['cellLabels']
cellOutlines = cells['cellOutlines']
shifts_allrounds = cells['shifts_allrounds']
dapis_shifted = cells['dapis_shifted']
nR = cells['nR']

# print(f'shifts_allrounds: {shifts_allrounds}')

spots = pd.read_excel(os.path.join(path, 'result/result.xlsx'), index_col=0).astype(np.uint16)
spots_assigned_allrounds = []
for r in range(nR):
    spots_assigned = []
    for c in range(nC):
        coords = spots.loc[(spots['round']==r+1) & (spots['channel']==c), ['z-coord', 'y-coord', 'x-coord']].to_numpy()
        spots_assigned.append(coords)
    spots_assigned_allrounds.append(spots_assigned)

filepath = os.path.join(path, '*.tif')    
filenames = sorted(glob.glob(filepath), key=os.path.basename)


if ismip:

    img_first = image_read(filenames[0])
    
    nC, size_z, size_y, size_x = img_first.shape

    mips = [image_mip(img_first, axis=1)]

    for filename in filenames[1:-1]:
        mips.append(image_mip(image_read(filename), axis=1))

    intensity_range = []
    for ch in range(nC):
        ch_min = np.min(np.array([np.min(mip[ch]) for mip in mips]))
        ch_max = np.max(np.array([np.max(mip[ch]) for mip in mips]))
        # ch_min = np.round(np.mean(np.array([np.quantile(mip[ch], 0.25) for mip in mips])))
        # ch_max = np.round(np.mean(np.array([np.quantile(mip[ch], 0.75) for mip in mips])))
        intensity_range.append([ch_min, ch_max])

    mips_warped = []
    for r, mip in enumerate(mips):
        mip_warped = mip.copy()
        for ch in range(nC):
            mip_warped[ch] = ndimage.shift(mip[ch], shift=(shifts_allrounds[r][ch][-2], shifts_allrounds[r][ch][-1]))
            
        mips_warped.append(mip_warped)

    viewer = napari.Viewer()
    viewer.add_image(
        image_mip(imgCells[:,:,:,1], axis=0),
        name='dapi reference',
        multiscale=False,
        visible=False
    )
    viewer.add_image(
        image_mip(imgCells[:,:,:,0], axis=0),
        name='cell reference',
        multiscale=False,
        visible=False
    )

    # for r, mip in enumerate(mips_warped):
    #     for ch in range(nC):
    #         viewer.add_image(
    #             mip[ch],
    #             name=f'origianl, {r+1} round, {ch+1} ch',
    #             multiscale=False,
    #             visible=False
    #         )

    viewer.add_labels(
        np.max(cellLabels, axis=0),
        name='cell labels',
        multiscale=False
    )

    cmap = matplotlib.cm.get_cmap('turbo', (nR-1)*(nC-1))
    colors = [matplotlib.colors.rgb2hex(c) for c in cmap.colors]
    for r, (mip, spots_assigned) in enumerate(zip(mips_warped, spots_assigned_allrounds)):
        for c, spots in enumerate(spots_assigned):
            spots = np.array(spots)
            viewer.add_image(
                mip[c],
                name=f'original, {r+1} round, {c+1} ch',
                multiscale=False,
                visible=False,
                contrast_limits=[intensity_range[c][0], intensity_range[c][1]]
            )
            if spots.shape[0] > 0:
                # print(spots.shape)
                # pointProperties = {
                #     'good_point': np.ones((spots.shape[0],), dtype=bool),
                #     'confidence': np.full((spots.shape[0],), cmap[cIndex])
                # }
                # print(cmap[cIndex])
                color = random.sample(colors,1)
                viewer.add_points(
                    spots[:,1:],
                    face_color=color,
                    # edge_color=color,
                    # symbol='ring',
                    size=10,
                    name=f'{r+1} round, {c+1} ch',
                    # visible=False
                )

    napari.run()
    

else:   # 3d visualization
    zyxScale = (zToXYRatioReal, 1, 1)
    single_round_for_original_img = -1

    img = image_read(filenames[single_round_for_original_img])
    img_warped = [image_warp(img[ch], shift=shifts_allrounds[single_round_for_original_img][ch]) for ch in range(5)]

    viewer = napari.Viewer()
    # check registration across rounds by using cropped dapi
    for r in range(dapis_shifted.shape[0]):
        viewer.add_image(
            dapis_shifted[r],
            name=f'{r+1} round dapi',
            scale=zyxScale,
            # contrast_limits=[dapis_shifted[r].min(), dapis_shifted[r].max()],
            multiscale=False,
            visible=False
        )
    # napari.run()

    # viewer = napari.Viewer()
    viewer.add_image(
        imgCells[:,:,:,1],
        name='dapi reference',
        scale=zyxScale,
        contrast_limits=[0, imgCells[:,:,:,1].max()],
        multiscale=False,
        visible=False
    )
    viewer.add_image(
        imgCells[:,:,:,0],
        name='cell reference',
        scale=zyxScale,
        contrast_limits=[0, imgCells[:,:,:,0].max()],
        multiscale=False,
        visible=False
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

    for ch, img in enumerate(img_warped):
        viewer.add_image(
            img,
            name=f'original, {ch} ch',
            scale=zyxScale,
            multiscale=False,
            visible=False
        )

    cmap = matplotlib.cm.get_cmap('turbo', (nR-1)*(nC-1))
    colors = [matplotlib.colors.rgb2hex(c) for c in cmap.colors]
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
                color = random.sample(colors, 1)
                viewer.add_points(
                    spots,
                    face_color=color,
                    edge_color=color,
                    size=5,
                    n_dimensional=True,
                    name=f'{r+1} round, {c+1} ch',
                    scale=zyxScale,
                    blending='additive',
                    # visible=False
                )
                cIndex = cIndex + 1

    napari.run()