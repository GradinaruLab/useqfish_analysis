import numpy as np
from skimage.feature import blob_log
import faulthandler

faulthandler.enable()
from functools import reduce
import operator
from scipy.spatial.distance import cdist
import gc
import math


def blob_detection(
    img, shift=None, minSigma=1, maxSigma=10, numSigma=10, threshold=0.1, overlap=0
):
    """ """
    ## 2d to 3d stitch
    if shift is None:
        shift = np.array([0, 0, 0])

    spotCoordinates = []
    for z in range(img.shape[0]):
        spots = blob_log(
            np.squeeze(img[z]),
            min_sigma=minSigma,
            max_sigma=maxSigma,
            num_sigma=numSigma,
            threshold=threshold,
            # threshold_rel=threshold_rel,
            overlap=overlap,
        ).astype(np.uint16)

        if spots.size > 0:
            spots = np.pad(spots, ((0, 0), (1, 0)), mode="constant", constant_values=z)

        spotCoordinates.append(list(spots))

    spotCoordinates = reduce(
        operator.iconcat, spotCoordinates, []
    )  # make a flatten list
    spotCoordinates = np.array(spotCoordinates)
    # print(f'spotCoordinates.type: {spotCoordinates.dtype}')
    if spotCoordinates.shape[0] > 0:
        spotNewCoordinates = spot_stitch_3d(spotCoordinates, 1)
        spotCoordsMoved = spot_warp(spotNewCoordinates, shift=shift)
        spotCoordsNoBoundary = spot_remove_out_of_boundary(
            img.shape, spotCoordsMoved, maxSigma
        )
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
    if radius < 1:
        radius = 1

    zrange = list(np.unique(spotCoordinates[:, 0]))

    spotLabels = np.zeros((spotCoordinates.shape[0],), dtype=np.uint16)
    # for z in range(nZ):
    for ind, z in enumerate(zrange):
        spotIdx = np.where(spotCoordinates[:, 0] == z)[0]
        # print(f'spotIdx.shape: {spotIdx.shape}')
        if spotIdx.size > 0:
            if ind == 0:
                spotLabels[spotIdx] = np.arange(1, spotIdx.size + 1)
                maxLabel = spotIdx.size + 1
                preLabel = np.arange(1, spotIdx.size + 1)
            else:
                preIdx = np.where(spotCoordinates[:, 0] == z - 1)[0]
                if preIdx.size > 0:
                    distances = cdist(
                        spotCoordinates[preIdx, 1:3], spotCoordinates[spotIdx, 1:3]
                    )
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
                    spotLabels[spotIdx] = np.arange(maxLabel, spotIdx.size + maxLabel)
                    maxLabel = maxLabel + spotIdx.size
                    preLabel = np.arange(maxLabel, spotIdx.size + maxLabel)

    spotNewCoordinates = _stitch_3d(spotCoordinates, spotLabels)

    del spotLabels
    gc.collect()

    return spotNewCoordinates


# @numba.njit(parallel=True)
def _stitch_3d(spotCoords, spotLabels):
    nLabels = np.int64(np.max(spotLabels))
    spotNewCoords = np.zeros((nLabels, 4), dtype=np.uint16)
    for i in range(nLabels):
        spotIdx = np.where(spotLabels == i)[0]
        spotNewCoords[i, :-1] = np.floor(
            np.median(spotCoords[spotIdx, :3], axis=0)
        ).astype(np.uint16)
        spotNewCoords[i, -1] = i

    return spotNewCoords


def spot_remove_out_of_boundary(imgShape, spotCoords, radius):
    nZ, nY, nX = imgShape

    if not isinstance(spotCoords, list):
        spotCoords = list(spotCoords)

    newCoords = []
    for coord in spotCoords:
        zyxCoord = coord[:-1]
        label = coord[-1]

        zyxCoord = np.where(
            # zyxCoord >= radius,
            zyxCoord >= 0,
            zyxCoord,
            np.NaN
            # radius
        )
        zyxCoord = np.where(
            # zyxCoord < (np.array([nZ, nY, nX])-radius),
            zyxCoord < np.array([nZ, nY, nX]),
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
        spotAssigned.append(np.array([z, y, x, cellLabels[z, y, x]]))
    return spotAssigned


def spot_warp(coords, shift=None):
    if shift is None:
        shift = np.array([0, 0, 0])

    newCoords = np.column_stack(
        (
            coords[:, 0] + shift[0],
            coords[:, 1] + shift[1],
            coords[:, 2] + shift[2],
            coords[:, 3],
        )
    )
    return newCoords


def isnoise(pre_coords, cur_coords, dist_threshold=1):
    if math.dist(pre_coords, cur_coords) <= dist_threshold:
        return True
    else:
        return False
