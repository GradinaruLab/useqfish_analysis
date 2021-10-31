# unmute when the error "libgcc_s.so.1 must be installed for pthread_cancel to work" came up
# import ctypes
# libgcc_s = ctypes.CDLL('libgcc_s.so.1')

from cell_detection import *
from spot_detection import *
from image_manipulation import *
from params import *

import pandas as pd
from pylab import *
# import napari

import argparse
import os
import sys
import glob
import dask
import dask.array as da
from dask.diagnostics import ProgressBar

import zarr
import numcodecs

import warnings

from pickling import *

warnings.filterwarnings("ignore")

my_parser = argparse.ArgumentParser(description='Run analysis on a group of images')

my_parser.add_argument('Path',
                       metavar='path',
                       type=str,
                       help='the path to the diretory containing the images')


# Parse arguments
args = my_parser.parse_args()

path = args.Path
filepath = os.path.join(path, '*.tif')    
filenames = sorted(glob.glob(filepath), key=os.path.basename)
print(filenames)
nR = len(filenames) # number of rounds

# 1. 3d spot detection
imgReference = image_read(filenames[roundRef])
dapi_reference = img_as_float32(imgReference[0])
size_z, size_y, size_x = dapi_reference.shape

imgMetadata = read_ome_metadata(filenames[roundRef])
zToXYRatioReal = imgMetadata['zRealSize']/imgMetadata['xRealSize']


print(f'>> STEP 1. Cell detection -')
img_cells = np.stack((imgReference[cellch], imgReference[0]), axis=3)
cellLabels, cellOutlines, zToXYRatioReal, nCells = cell_detection(
    img_cells, 
    zToXYRatioReal=zToXYRatioReal,
    resizeFactor=0.2, 
    color_shift=color_shifts[0]
)

print(f'>> STEP 2. registration - ')
spots_assigned_allrounds = []
for filename in filenames[:roundRef]:
    print(filename)
    img = image_read(filename)
    img = [img[c] for c in range(nC)]
    img = [img_as_float32(ch) for ch in img]

    round_shift = image_shift(
        img[0][:, int(size_y/2-shift_window_size):int(size_y/2), int(size_x/2-shift_window_size):int(size_x/2)],
        dapi_reference[:, int(size_y/2-shift_window_size):int(size_y/2), int(size_x/2-shift_window_size):int(size_x/2)]
    )
    
    shifts = [np.array(round_shift) + np.array(color_shift) for color_shift in color_shifts[1:]]
    # shifts = [round_shift for c in range(nC-1)]
    # img_warpped = [image_warp(ch, shift=shift) for ch, shift in zip(img, shifts)]
    
    # exit()
    print(shifts)


    print(f'>> STEP 2-2. Spot detection -')
    # set up dask for running in parallel
    daimg = [da.from_array(ch, chunks=(1, -1, -1)) for ch in img[1:]]
    
    daimg = [da.map_blocks(median_filter, ch) for ch in daimg]
    daimg = [da.map_blocks(background_subtraction, ch, size=100) for ch in daimg]

    with ProgressBar():
        daimg = [ch.compute() for ch in daimg]
    # print(daimg[0].shape)

    img_delayed = [dask.delayed(ch) for ch in daimg]

    spots = [
        dask.delayed(blob_detection)(
            ch, shift=shift, minSigma=sigma, maxSigma=sigma, numSigma=1, threshold=0.005
        )
        for ch, shift in zip(img_delayed, shifts)
    ]

    spots_assigned = [dask.delayed(spot_assignment)(spot, cellLabels) for spot in spots]

    with ProgressBar():
        # spots = list(dask.compute(*spots))
        spots_assigned = list(dask.compute(*spots_assigned))

    spots_assigned_allrounds.append(spots_assigned)


print(f'>> STEP 4. Save results -')
spots_results = []
for r, spots_assigned in enumerate(spots_assigned_allrounds):
    for c, spots in enumerate(spots_assigned):
        for spot in spots:
            spots_results.append(np.append(np.array([r+1, c+1]), spot))

# for r, (imgRound, spotRound) in enumerate(zip(imgSpotsAssigned, spotAssigned)):
#     for c, (img, spots) in enumerate(zip(imgRound, spotRound)):
#         imgResults.append(img>0)
#         chNames.append(f'round #{r}, channel #{c+1}')
#         for spot in spots:
#             spotResults.append(np.append(np.array([r, c+1]), spot))

# image save for visualization on local 

zarr.save(
    os.path.join(path, 'result/result_cellsegments.zarr'),
    imgCells=img_cells,
    cellLabels=cellLabels,
    cellOutlines=cellOutlines,
    zToXYRatioReal=zToXYRatioReal,
    nR=nR
)

# spotCoords = np.array([np.array(coords) for coords in spotCoords])
spots_results = np.array(spots_results)
print(f'>>>> Total {spots_results.shape[0]} spots detected')
# print(f'{spotCoords.shape}')

resultDf = pd.DataFrame({
    'round': spots_results[:,0],
    'channel': spots_results[:,1],
    'z-coord': spots_results[:,2],
    'y-coord': spots_results[:,3],
    'x-coord': spots_results[:,4],
    'cell id': spots_results[:,5]
})
resultDf.to_excel(excel_writer=os.path.join(path, 'result/result.xlsx'))

target_index = np.zeros((nR, nC), dtype=np.int)
index = 0
for r in range(nR-1):
    for c in range(nC-1):
        target_index[r, c] = index
        index = index + 1

spots_per_cell = np.zeros((nCells+1, (nR-1)*(nC-1)), dtype=np.int)
for spot_index in range(spots_results.shape[0]):
    cell_id = int(spots_results[spot_index,-1])
    r = spots_results[spot_index,0]
    c = spots_results[spot_index,1]
    spots_per_cell[cell_id, target_index[r-1, c-1]] = spots_per_cell[cell_id, target_index[r-1, c-1]]+1

resultDf2 = pd.DataFrame(data=spots_per_cell)
resultDf2.to_excel(excel_writer=os.path.join(path, 'result/result_spots_per_cell.xlsx'))