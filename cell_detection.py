from skimage.morphology import convex_hull_image
from skimage.transform import resize
from skimage.segmentation import expand_labels
from skimage import img_as_float32
import numpy as np
from cellpose.models import Cellpose
import gc
from tqdm import tqdm

from image_manipulation import (
    image_downsample_shape,
    image_normalize_layers,
    image_gaussian_filter,
)


def cell_detection(
    img, zToXYRatioReal=1, resizeFactor=0.2, medianDiameters=np.float32(30.0)
):
    """
    Main function for cell detection using cellpose.
    """

    zToXYRatio = zToXYRatioReal * resizeFactor
    img_cells = img_as_float32(img)

    imgResizedShaped = image_downsample_shape(img_cells, resizeFactor=resizeFactor)
    imgNormalized = image_normalize_layers(imgResizedShaped)
    imgFiltered = image_gaussian_filter(imgNormalized, sigma=2)

    mask, _, _, _ = run_cellpose(
        imgFiltered,
        [0, 0],
        gpu=True,
        # device=mx.gpu(1),
        anisotropy=zToXYRatio,
        diameter=medianDiameters,
        # diameter=None,
        cellprob_threshold=-5,
        flow_threshold=0.6,
        min_size=1000,  # min_size is not actually working well. how does it work with 3d, stitch_threshold
        do_3D=True,
    )

    nCells = np.unique(mask).size - 1
    print(f">>>> {nCells} of cells detected")
    print(f">>>> median of diameters: {medianDiameters}")

    maskClosed = mask_close(mask.copy())
    maskUpsampled = mask_upsample(
        maskClosed, (img.shape[0], img.shape[1], img.shape[2])
    )

    del imgResizedShaped, imgNormalized, imgFiltered, mask, maskClosed
    gc.collect()

    return maskUpsampled, zToXYRatioReal, nCells


def run_cellpose(
    img,
    channels,
    gpu=True,
    model_type="cyto2",
    diameter=30.0,
    device=None,
    anisotropy=None,
    stitch_threshold=0.25,
    cellprob_threshold=0.0,
    flow_threshold=0.4,
    min_size=15,
    do_3D=False,
):
    """
    Uses cellpose to segment cells. Returns mask, flow, style and diam of cells.
    """
    model = Cellpose(gpu=gpu, model_type=model_type)
    mask, flow, style, diam = model.eval(
        img,
        channels=channels,
        do_3D=do_3D,
        diameter=diameter,
        anisotropy=anisotropy,
        stitch_threshold=stitch_threshold,
        cellprob_threshold=cellprob_threshold,
        flow_threshold=flow_threshold,
        min_size=min_size,
    )

    return mask, flow, style, diam


def mask_close(mask):
    """
    Closes the mask. TODO: What does this actually do?
    """
    maskClosed = np.zeros_like(mask)
    for i in tqdm(range(1, mask.max() + 1)):
        maskI = mask == i
        for z in range(maskI.shape[0]):
            maskIHull = convex_hull_image(maskI[z])
            y, x = np.nonzero(maskIHull)
            maskClosed[z, y, x] = i
        for y in range(maskI.shape[1]):
            maskIHull = convex_hull_image(maskI[:, y, :])
            z, x = np.nonzero(maskIHull)
            maskClosed[z, y, x] = i
        for x in range(maskI.shape[2]):
            maskIHull = convex_hull_image(maskI[:, :, x])
            z, y = np.nonzero(maskIHull)
            maskClosed[z, y, x] = i

    return maskClosed


def mask_upsample(mask, finalShape=None):
    """
    Resizes the mask to the finalShape with datatype np.unit16.
    """
    if finalShape is None:
        raise ValueError("please provide the final shape of the mask")
    elif mask.ndim != len(finalShape):
        raise ValueError("dimension of final shape is not matched with input mask")
    else:
        maskUpsampled = resize(mask, finalShape, order=0, preserve_range=True)

        return maskUpsampled.astype(np.uint16)


def mask_expanded(mask, expandDist=1):
    """
    Expands a mask so that regions within expandDist of a label get the same label
    """
    maskExpanded = np.zeros_like(mask)
    for z in tqdm(range(mask.shape[0])):
        maskExpanded[z] = expand_labels(mask[z], distance=expandDist)
    return maskExpanded.astype(np.uint16)
