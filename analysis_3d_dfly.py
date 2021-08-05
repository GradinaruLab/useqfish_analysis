import os
import sys
import warnings
import glob
import time
import gc
import functools
import operator
# import numba
from tqdm import tqdm

import numpy as np
import zarr
import pandas as pd

import dask
import dask.array as da
from dask.diagnostics import ProgressBar

import tifffile
from skimage import filters, feature, img_as_float32, registration, transform, morphology, segmentation
from scipy.ndimage import white_tophat
from scipy.spatial import distance

from cellpose import models, utils, transforms
import mxnet as mx

from oxdls import OMEXML

import napari
from pylab import *


def measure_run_time(func):
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        print(f'func({func.__name__}) running time: {end-start} sec')
        return result
    return wrapper

def read_ome_metadata(filePath):
    """
    Args: file path of ome-tiff

    Returns: dictionary of parsed metadata
    """
    with tifffile.TiffFile(filePath) as tif:
        imgMetadata = OMEXML(tif.ome_metadata)
    
    dictMetadata = {'nChannels': imgMetadata.image().Pixels.channel_count,
                    'xPixels': imgMetadata.image().Pixels.SizeX,
                    'yPixels': imgMetadata.image().Pixels.SizeY,
                    'zPixels': imgMetadata.image().Pixels.SizeZ,
                    'xRealSize': imgMetadata.image().Pixels.PhysicalSizeX,
                    'yRealSize': imgMetadata.image().Pixels.PhysicalSizeY,
                    'zRealSize': imgMetadata.image().Pixels.PhysicalSizeZ,
                    'channelNames': [imgMetadata.image().Pixels.Channel(ch).Name.split('_')[0] for ch in range(imgMetadata.image().Pixels.channel_count)]
    }

    return dictMetadata

# @measure_run_time
def image_read(filename):
    """
    read ome-tiff file as numpy array, converted to [0.0 1.0] float32 type
    if dapi channel is not the first, bring the dapi channel to the first
    """
    # if chDapi != 0:
    #     img = np.zeros_like(imgLoaded)
    #     img[0] = imgLoaded[chDapi]
    #     newIdx = 1
    #     for c in range(imgLoaded.shape[0]):
    #         if c != chDapi:
    #             img[newIdx] = imgLoaded[c]
    #             newIdx = newIdx + 1
    #     imgLoaded = img

    #     del img
    #     gc.collect()

    # return imgLoaded
    return tifffile.imread(filename)

# @numba.jit
def median_filter(img):
    """
    """
    img = img[0,...]
    filtered = filters.median(img)
    # imgMedianFiltered = np.array(
    #     [filters.median(img[z]) for z in range(img.shape[0])]
    # )
    # imgMedianFiltered = np.zeros_like(img)
    # for z in range(img.shape[0]):
    #     imgMedianFiltered[z] = filters.median(img[z])

    # return imgMedianFiltered[None,...]
    # return imgMedianFiltered
    return filtered[None,...]

# @numba.stencil
# def median_for_numba(x):
#     return np.median(x[-1, -1] + x[-1, 0] + x[-1, 1] +
#             x[ 0, -1] + x[ 0, 0] + x[ 0, 1] +
#             x[ 1, -1] + x[ 1, 0] + x[ 1, 1])

def background_subtraction(img, size=10, mode='nearest'):
    """
    """
    img = img[0,...]
    filtered = white_tophat(img, size=size, mode=mode)
    # imgTophat = np.array(
    #     [white_tophat(img[z], size=size, mode=mode) for z in range(img.shape[0])]
    # )

    # imgTophat = np.zeros_like(img)
    # for z in range(img.shape[0]):
    #     imgTophat[z] = white_tophat(img[z], size=size, mode=mode)
    # imgTophat = imgTophat[None, ...]
    
    # return imgTophat[None,...]
    # return imgTophat
    return filtered[None,...]

def blob_detection(img, shift=None, minSigma=1, maxSigma=10, numSigma=10, threshold=0.1, overlap=0):
    """
    """
    ## 2d to 3d stitch
    if shift is None:
        shift = np.array([0, 0, 0])

    spotCoordinates = []
    for z in range(img.shape[0]):
        spots = feature.blob_log(
            np.squeeze(img[z]),
            min_sigma=minSigma,
            max_sigma=maxSigma,
            num_sigma=numSigma,
            threshold=threshold,
            overlap=overlap
        ).astype(np.uint16)

        if spots.size > 0:
            spots = np.pad(spots, ((0,0),(1,0)), mode='constant', constant_values=z)

        spotCoordinates.append(list(spots))
            
    spotCoordinates = functools.reduce(operator.iconcat, spotCoordinates, [])   # make a flatten list
    spotCoordinates = np.array(spotCoordinates)
    # print(f'spotCoordinates.shape: {spotCoordinates.shape}')
    if spotCoordinates.shape[0] > 0:
        spotNewCoordinates = spot_stitch_3d(spotCoordinates, np.sqrt(2)*maxSigma)
        spotCoordsMoved = spot_warp(spotNewCoordinates, shift=shift)
        spotCoordsNoBoundary = spot_remove_out_of_boundary(img.shape, spotCoordsMoved, maxSigma)
    else:
        spotCoordsNoBoundary = spotCoordinates
    # print(spotNewCoordinates.shape, spotCoordsMoved.shape)
    # labelSpots = make_spots_label(img.shape, spotCoordsMoved, maxSigma)
    
    # return imgSpots, spotCoordinates

    # del spots, spotCoordinates, spotNewCoordinates, spotCoordsMoved
    del spots, spotCoordinates
    gc.collect()

    # return labelSpots, spotCoordsMoved
    return spotCoordsNoBoundary

def spot_stitch_3d(spotCoordinates, radius):
    # nZ = spotCoordinates[:,0].max()
    zrange = list(np.unique(spotCoordinates[:,0]))

    spotLabels = np.zeros((spotCoordinates.shape[0],), dtype=np.uint16)
    # for z in range(nZ):
    for ind, z in enumerate(zrange):
        spotIdx = np.where(spotCoordinates[:,0] == z)[0]
        # print(f'spotIdx.shape: {spotIdx.shape}')
        if spotIdx.size > 0:
            if ind==0:
                spotLabels[spotIdx] = np.arange(1, spotIdx.size+1)
                maxLabel = spotIdx.size + 1
                preLabel = np.arange(1, spotIdx.size+1)
            else:
                preIdx = np.where(spotCoordinates[:,0] == z-1)[0]
                if preIdx.size > 0:
                    distances = distance.cdist(spotCoordinates[preIdx, 1:3], spotCoordinates[spotIdx, 1:3])
                    # print(f'distances.shape: {distances.shape}')
                    minDist = np.amin(distances, axis=0)
                    # print(f'minDist.shape: {minDist.shape}')
                    minIdx = np.argmin(distances, axis=0)
                    newLabel = np.zeros((spotIdx.size,))
                    # print(minDist, minIdx)
                    for i, (md, mi) in enumerate(zip(list(minDist), list(minIdx))):
                        if md <= radius:
                            newLabel[i] = preLabel[mi]
                        else:
                            newLabel[i] = maxLabel
                            maxLabel = maxLabel + 1
                    spotLabels[spotIdx] = newLabel
                    preLabel = newLabel
                else:
                    spotLabels[spotIdx] = np.arange(maxLabel, spotIdx.size+maxLabel)
                    maxLabel = maxLabel+spotIdx.size
                    preLabel = np.arange(maxLabel, spotIdx.size+maxLabel)

    spotNewCoordinates = _stitch_3d(spotCoordinates, spotLabels)

    del spotLabels
    gc.collect()

    return spotNewCoordinates

# @numba.njit(parallel=True)
def _stitch_3d(spotCoords, spotLabels):
    nLabels = np.int64(np.max(spotLabels))
    spotNewCoords = np.zeros((nLabels,4), dtype=np.uint16)
    for i in range(nLabels):
        spotIdx = np.where(spotLabels==i)[0]
        spotNewCoords[i,:-1] = np.floor(np.median(spotCoords[spotIdx,:3], axis=0)).astype(np.uint16)
        spotNewCoords[i, -1] = i

    return spotNewCoords

# @measure_run_time
def make_spots_mask(imgShape, spotCoordinates):
    _, nY, nX = imgShape
    maskSpots = np.zeros(imgShape, dtype=bool)

    if not isinstance(spotCoordinates, list):
        spotCoordinates = list(spotCoordinates)

    for coord in spotCoordinates:
        xyCoord = coord[1:-1]
        radius = coord[-1]

        xyCoord = np.where(xyCoord>=radius, xyCoord, radius)
        xyCoord = np.where(
            xyCoord < (np.array([nY, nX])-radius),
            xyCoord,
            np.array([nY, nX])-radius-1
        )

        maskSpots[
            coord[0],
            xyCoord[0]-radius:xyCoord[0]+radius+1,
            xyCoord[1]-radius:xyCoord[1]+radius+1
        ] = True

    return maskSpots

# @numba.njit(parallel=True)
def make_spots_label(imgShape, spotCoords, radius):
    labelSpots = np.zeros(imgShape, dtype=np.uint16)
    
    for coord in spotCoords:
        labelSpots[
            coord[0]-radius:coord[0]+radius+1,
            coord[1]-radius:coord[1]+radius+1,
            coord[2]-radius:coord[2]+radius+1
        ] = coord[-1]

    return labelSpots

@measure_run_time
def image_downsample_shape(img, zToXYRatioReal, resizeFactor=0.2):
    # nC, nZ, nX, nY = img.shape
    nZ, nY, nX, nC = img.shape
    # imgdTFirst = np.zeros((nZ, nX, nY, nC), dtype=img.dtype)
    # imgdTFirst[:,:,:,0] = img[1].copy()
    # imgdTFirst[:,:,:,1] = img[0].copy()
    # imgFloat = img_as_float32(imgdTFirst)
    imgResized = transforms.resize_image(img, rsz=resizeFactor)
    # zToXYRatio = zToXYRatioReal*resizeFactor
    
    _, nXResized, nYResized, _  = imgResized.shape
    imgResizedShaped = np.zeros((nZ, nC, nXResized, nYResized), dtype=imgResized.dtype)
    for z in range(nZ):
        for c in range(nC):
            imgResizedShaped[z, c, :, :] = imgResized[z, :, :, c]
    
    # print(imgResizedShaped.shape)
    return imgResizedShaped

@measure_run_time
# @numba.njit
def image_normalize_layers(img):
    imgNormalized = np.zeros_like(img)
    nZ, nC, _, _ = img.shape

    meanTotal = [np.mean(img[:,c,:,:].ravel()[np.nonzero(img[:,c,:,:].ravel())]) for c in range(nC)]
    stdTotal = [np.std(img[:,c,:,:].ravel()[np.nonzero(img[:,c,:,:].ravel())]) for c in range(nC)]

    for z in range(nZ):
        for c in range(nC):
            flattenLayer = img[z,c,:,:].ravel()
            meanLayer = np.mean(flattenLayer[np.nonzero(flattenLayer)])
            stdLayer = np.std(flattenLayer[np.nonzero(flattenLayer)])
            imgNormalized[z,c,:,:] = meanTotal[c] + (img[z,c,:,:] - meanLayer) * (stdTotal[c]/stdLayer)

    return imgNormalized

@measure_run_time
def image_gaussian_filter(img, sigma=1):
    imgFiltered = np.zeros_like(img)
    nZ, nC, _, _ = img.shape

    for z in range(nZ):
        for c in range(nC):
            imgFiltered[z,c,:,:] = filters.gaussian(img[z,c,:,:], sigma=sigma)
    
    return imgFiltered

def run_cellpose(
    img, 
    channels, 
    gpu=True, 
    model='cyto', 
    diameter=30, 
    device=None, 
    anisotropy=None,
    stitch_threshold=0.25,
    cellprob_threshold=0.0,
    flow_threshold=0.4,
    min_size=15,
    do_3D=False):
    model = models.Cellpose(gpu=gpu, model_type=model, device=device, torch=False)
    mask, flow, style, diam = model.eval(
        img, 
        channels=channels, 
        do_3D=do_3D,
        diameter=diameter,
        anisotropy=anisotropy,
        stitch_threshold=stitch_threshold,
        cellprob_threshold=cellprob_threshold,
        flow_threshold=flow_threshold,
        min_size=min_size
    )
    return mask, flow, style, diam

def image_with_outlines(img, mask):
    outlines = utils.masks_to_outlines(mask)
    outZ, outY, outX = np.nonzero(outlines)
    imgOutlined = np.zeros((img.shape[0], img.shape[1], img.shape[2], 3))
    imgOutlined[outZ, outY, outX] = np.array([1, 1, 1])
    return imgOutlined

@measure_run_time
def mask_closing(mask, selem=morphology.ball(2), expand=True, expandDist=2):
    maskClosed = np.zeros_like(mask)

    # for i in range(1, mask.max()+1):
    #     maskI = (mask==i)
    #     maskIClosed = morphology.binary_closing(maskI, selem=selem)
    #     # maskIClosed = morphology.binary_erosion(maskIClosed, selem=selem)
    #     # maskIClosed = morphology.binary_dilation(maskIClosed, selem=expand)

    #     z, y, x = np.nonzero(maskIClosed)
    #     # for c in range(z.size):
    #     #     if mask[z[c],y[c],x[c]]==0 or mask[z[c],y[c],x[c]]==i:
    #     #         maskClosed[z[c],y[c],x[c]] = i
                
    #     maskClosed[z,y,x] = i

    for i in tqdm(range(1, mask.max()+1)):
        maskI = (mask==i)
        for z in range(maskI.shape[0]):
            maskIHull = morphology.convex_hull_image(maskI[z])
            y, x = np.nonzero(maskIHull)
            maskClosed[z,y,x] = i
        for y in range(maskI.shape[1]):
            maskIHull = morphology.convex_hull_image(maskI[:,y,:])
            z, x = np.nonzero(maskIHull)
            maskClosed[z,y,x] = i
        for x in range(maskI.shape[2]):
            maskIHull = morphology.convex_hull_image(maskI[:,:,x])
            z, y = np.nonzero(maskIHull)
            maskClosed[z,y,x] = i

    if expand:
        maskClosed = segmentation.expand_labels(maskClosed, distance=expandDist)

    return maskClosed

@measure_run_time
def mask_upsample(mask, finalShape=None):
    if finalShape is None:
        raise ValueError('please provide the final shape of the mask')
    elif mask.ndim != len(finalShape):
        raise ValueError('dimension of final shape is not matched with input mask')
    else:
        maskUpsampled = transform.resize(mask, finalShape, order=0, preserve_range=True)
        
        return maskUpsampled.astype(np.int)

@measure_run_time
def cell_detection(filename, cellCh=1, resizeFactor=0.2, color_shift=None):
    """
    main function for cell detection using cellpose
    """
    img = image_read(filename)
    imgMetadata = read_ome_metadata(filename)
    zToXYRatioReal = imgMetadata['zRealSize']/imgMetadata['xRealSize']
    zToXYRatio = zToXYRatioReal * resizeFactor
    
    img_cells = np.stack((img[cellCh], img[0]), axis=3)
    img_cells = img_as_float32(img_cells)
   
    if color_shift is not None:
        img_cells[1] = image_warp(img_cells[1], color_shift)

    imgResizedShaped = image_downsample_shape(img_cells, zToXYRatioReal, resizeFactor=resizeFactor)
    imgNormalized = image_normalize_layers(imgResizedShaped)
    imgFiltered = image_gaussian_filter(imgNormalized, sigma=2)

    # print(imgFiltered.shape)

    # _, _, _, diamEstimated = run_cellpose(
    #     imgFiltered,
    #     [0,0],
    #     gpu=True,
    #     # device=mx.gpu(1),
    #     anisotropy=zToXYRatio,
    #     diameter=None,
    #     cellprob_threshold=-5,
    #     flow_threshold=0.6,
    #     stitch_threshold=0.1
    # )
    # medianDiameters = np.median(diamEstimated)
    medianDiameters = 30

    mask, _, _, _ = run_cellpose(
        # imgResizedShaped,
        imgFiltered,
        [0,0],
        gpu=True,
        # device=mx.gpu(1),
        anisotropy=zToXYRatio,
        diameter=medianDiameters,
        cellprob_threshold=-5,
        flow_threshold=0.6,
        min_size=1000, #min_size is not actually working well. how does it work with 3d, stitch_threshold
        do_3D=True,
        # stitch_threshold=0.1 # not necessarily critical
    )    
    
    nCells = np.unique(mask).size - 1
    print(f'>>>> {nCells} of cells detected')
    print(f'>>>> median of diameters: {medianDiameters}')

    # outlineRGB = image_with_outlines(np.squeeze(imgResizedShaped[:,0,:,:]), mask)

    ## nuclei detection
    # maskNuclei, flowNuclei, _, diamsNuclei = run_cellpose(
    #     np.stack((imgFiltered[:,1,:,:], np.zeros(imgFiltered[:,1,:,:].shape)), axis=1),
    #     [0,0],
    #     gpu=True,
    #     device=mx.gpu(1),
    #     anisotropy=zToXYRatio,
    #     diameter=24,
    #     cellprob_threshold=-5,
    #     flow_threshold=0.6,
    #     # stitch_threshold=0.25,
    #     model='nuclei',
    #     do_3D=True
    # )
    
    # outlineRGBNuclei = image_with_outlines(np.squeeze(imgResizedShaped[:,1,:,:]), maskNuclei)

    # print(f'zToXYRatio: {zToXYRatio}')
    # nNuclei = np.unique(maskNuclei).size - 1
    # print(f'>>>> {nNuclei} of nuclei detected')
    # medianDiametersNuclei = np.median(diamsNuclei)
    # print(f'>>>> median of nucleus diameters: {medianDiametersNuclei}')

    # centroids, maskMarker = make_nuclei_markers(maskNuclei)

    ## convex_hull
    # maskConvex = mask_convex_hull(mask)
    maskClosed = mask_closing(mask.copy())
    maskUpsampled = mask_upsample(maskClosed, (img.shape[1], img.shape[2], img.shape[3]))
    maskUpsampledOutlined = image_with_outlines(img[0], maskUpsampled)

    del imgResizedShaped, imgNormalized, imgFiltered, mask, maskClosed
    gc.collect()

    return maskUpsampled, maskUpsampledOutlined, zToXYRatioReal, nCells

@measure_run_time
def image_shift(movImg, refImgDapi):
    """
    return shift coordinates
    """
    # refImg = refImg[0,...]
    shift, _, _ = registration.phase_cross_correlation(refImgDapi, movImg)

    # return shift[None,...]
    return shift

@measure_run_time
def image_warp(img, shift=None):
    # img = img[0,...]
    if shift is None:
        shift = np.array([0,0,0])
    # else:
    #     shift = shift[0,...]

    z, y, x = img.shape
    newZ, newY, newX = np.mgrid[:z, :y, :x]
    newCoordinates = np.array([newZ-shift[0], newY-shift[1], newX-shift[2]])
    imgMoved = transform.warp(img, newCoordinates, mode='constant', cval=0)

    del newZ, newY, newX, newCoordinates
    gc.collect()

    # return imgMoved[None,...]
    return imgMoved

# @numba.njit
def spot_warp(coords, shift=None):
    if shift is None:
        shift = np.array([0,0,0])

    newCoords = np.column_stack((coords[:,0]-shift[0], coords[:,1]-shift[1], coords[:,2]-shift[2], coords[:,3]))
    return newCoords

def spot_remove_out_of_boundary(imgShape, spotCoords, radius):
    nZ, nY, nX = imgShape

    if not isinstance(spotCoords, list):
        spotCoords = list(spotCoords)
    
    newCoords = []
    for coord in spotCoords:
        zyxCoord = coord[:-1]
        label = coord[-1]

        zyxCoord = np.where(
            zyxCoord >= radius, 
            zyxCoord, 
            np.NaN
            # radius
        )
        zyxCoord = np.where(
            zyxCoord < (np.array([nZ, nY, nX])-radius),
            zyxCoord,
            np.NaN
            # np.array([nZ, nY, nX])-radius-1
        )
        if np.any(np.isnan(zyxCoord)):
            continue
        newCoords.append(np.concatenate((zyxCoord, label), axis=None).astype(np.uint16))

    return newCoords
    

def spot_assignment(spotCoords, cellLabels):
    # spotAssigned = [np.concatenate((coord[:-1], cellLabels[coord[0], coord[1], coord[2]]), axis=None) for coord in spotCoords]
    spotAssigned = []
    for coord in spotCoords:
        # print(coord[:-1])
        z, y, x = coord[:-1]
        spotAssigned.append(np.array([z, y, x, cellLabels[z,y,x]]))
    return spotAssigned

if __name__=="__main__":

    warnings.filterwarnings("ignore")

    path = sys.argv[1]
    # z0 = sys.argv[2]
    # zend = sys.argv[3]

    filepath = os.path.join(path, '*.tif')    
    filenames = glob.glob(filepath)
    nR = len(filenames)
    nC = 5
    roundRef = -1    # assuming the last round is the reference of cells and nuclei
    cellch = 3

    sigma = 3
    shift_window_size = 200

    color_shifts = [[8, 0, 0],          # color aberration shifts to 647 aquired by imaging focalcheck #1, 500nm
                    [-1, 0, 1],
                    [0, 0, 0],
                    [0, 0, 0],
                    [1, -3, 2]]

    # 1. 3d spot detection
    imgReference = image_read(filenames[roundRef])
    dapi_reference = img_as_float32(imgReference[0])
    size_z, size_y, size_x = dapi_reference.shape

    print(f'>> STEP 1. Cell detection -')
    cellLabels, cellOutlines, zToXYRatioReal, nCells = cell_detection(filenames[roundRef], cellCh=cellch, resizeFactor=0.2, color_shift=color_shifts[0])
  
    # zyxScale = (zToXYRatioReal, 1, 1)
    # with napari.gui_qt():
    #     viewer = napari.Viewer()
    #     # for img in results:
    #     #     viewer.add_image(img)

    #     viewer.add_image(
    #         imgReference[0],
    #         name='dapi reference',
    #         scale=zyxScale,
    #         contrast_limits=[0, imgReference[0].max()],
    #         multiscale=False
    #     )
    #     viewer.add_image(
    #         imgReference[cellch],
    #         name='cell reference',
    #         scale=zyxScale,
    #         contrast_limits=[0, imgReference[cellch].max()],
    #         multiscale=False
    #     )
    #     viewer.add_labels(
    #         cellLabels,
    #         name='cell labels',
    #         scale=zyxScale,
    #         multiscale=False
    #     )
    #     viewer.add_image(
    #         cellOutlines,
    #         name='cell outlines',
    #         scale=zyxScale,
    #         multiscale=False
    #     )

    # exit()

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
        
        print(f'>> STEP 2-2. Spot detection -')
        with ProgressBar():
            # spots = list(dask.compute(*spots))
            spots_assigned = list(dask.compute(*spots_assigned))

        spots_assigned_allrounds.append(spots_assigned)
        
        # with napari.gui_qt():
        #     viewer = napari.Viewer()

        #     for ch, spot in zip(daimg, spots):
        #         viewer.add_image(ch)
        #         viewer.add_points(np.array(spot)[:,:-1], size=10, n_dimensional=True, blending='additive')

        # exit()
    
    # 4. data save
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
    resultDf.to_excel(excel_writer=os.path.join(path, 'result.xlsx'))
    
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
    resultDf2.to_excel(excel_writer=os.path.join(path, 'result_spots_per_cell.xlsx'))
    # zarr.save_group(
    #     os.path.join(path, 'results.zarr'), 
    #     celllabels=cellLabels,
    #     celloutlines=cellOutlines,
    #     spots=np.array([np.array(img) for img in imgSpotsAssigned])
    # )

    # del imgSpotsAssigned, imgSpotsLabeled
    # gc.collect()

    # 5. visualization:
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