# %%
# from cell_detection import *
# from spot_detection import *
# from image_manipulation import *
from params_smhcr import *

import argparse
import os
import sys
import glob
import dask
import dask.array as da
from dask.diagnostics import ProgressBar

import tifffile

from matplotlib import pyplot as plt
import numpy as np
from skimage import data, util, transform, feature, measure, filters, metrics, exposure, img_as_float32

import warnings
warnings.filterwarnings("ignore")

# my_parser = argparse.ArgumentParser(description='Run analysis on a group of images')

# my_parser.add_argument('Path',
#                        metavar='path',
#                        type=str,
#                        help='the path to the diretory containing the images')


# # Parse arguments
# args = my_parser.parse_args()

# path = args.Path
# filepath = os.path.join(path, '*.tif')    
# filenames = sorted(glob.glob(filepath), key=os.path.basename)
# print(filenames)
# nR = len(filenames) # number of rounds
def match_locations(img0, img1, coords0, coords1, radius=5, sigma=3):
    """Match image locations using SSD minimization.

    Areas from `img0` are matched with areas from `img1`. These areas
    are defined as patches located around pixels with Gaussian
    weights.

    Parameters:
    -----------
    img0, img1 : 2D array
        Input images.
    coords0 : (2, m) array_like
        Centers of the reference patches in `img0`.
    coords1 : (2, n) array_like
        Centers of the candidate patches in `img1`.
    radius : int
        Radius of the considered patches.
    sigma : float
        Standard deviation of the Gaussian kernel centered over the patches.

    Returns:
    --------
    match_coords: (2, m) array
        The points in `coords1` that are the closest corresponding matches to
        those in `coords0` as determined by the (Gaussian weighted) sum of
        squared differences between patches surrounding each point.
    """
    y, x = np.mgrid[-radius:radius + 1, -radius:radius + 1]
    weights = np.exp(-0.5 * (x ** 2 + y ** 2) / sigma ** 2)
    weights /= 2 * np.pi * sigma * sigma

    match_list = []
    for r0, c0 in coords0:
        roi0 = img0[r0 - radius:r0 + radius + 1, c0 - radius:c0 + radius + 1]
        roi1_list = [img1[r1 - radius:r1 + radius + 1,
                          c1 - radius:c1 + radius + 1] for r1, c1 in coords1]
        # sum of squared differences
        ssd_list = [np.sum(weights * (roi0 - roi1) ** 2) for roi1 in roi1_list]
        match_list.append(coords1[np.argmin(ssd_list)])

    return np.array(match_list)


def compare(*images, **kwargs):
    """
    Utility function to display images side by side.
    
    Parameters
    ----------
    image0, image1, image2, ... : ndarrray
        Images to display.
    labels : list
        Labels for the different images.
    """
    f, axes = plt.subplots(1, len(images), **kwargs)
    axes = np.array(axes, ndmin=1)
    
    labels = kwargs.pop('labels', None)
    if labels is None:
        labels = [''] * len(images)
    
    for n, (image, label) in enumerate(zip(images, labels)):
        axes[n].imshow(image, interpolation='nearest', cmap='gray')
        axes[n].set_title(label)
        axes[n].axis('off')
    
    f.tight_layout()


def image_adjust(img):
    vmin, vmax = np.percentile(img, q=(0.5, 99.5))

    img_clipped = exposure.rescale_intensity(
        img,
        in_range=(vmin, vmax),
        out_range=np.float32
    )

    return img_clipped



montage_overlap = 1     # % overlap between images in the montage
# montage_shape = (3, 5)  # montage shape w x h
montage_shape = (1, 4)

## read reference image and metadata
# imgReference = image_read(filenames[roundRef])
# dapi_reference = imgReference[0]

# imgMetadata = read_ome_metadata(filenames[roundRef])
# zToXYRatioReal = imgMetadata['zRealSize']/imgMetadata['xRealSize']
# nC, size_z, size_y, size_x = imgReference.shape

# imgs = [dask.delayed(image_read)(filename) for filename in filenames]
# imgs_mip = [dask.delayed(img, axis=0) for img in imgs]
# dapis = [da.from_delayed(img[0], shape=(size_y, size_x), dtype=imgReference.type) for img in imgs_mip]
# dapis = da.stack(dapis)

montage_shape_w, montage_shape_h = montage_shape

filenames = sorted(glob.glob('./stitching/*.tif'), key=os.path.basename)
dapis = [img_as_float32(tifffile.imread(filename)) for filename in filenames]
y_size, x_size = dapis[0].shape
overlap_size = int(np.round(y_size * montage_overlap * 2/100))
compare(*dapis)

# %%
dapis = [image_adjust(dapi) for dapi in dapis]
compare(*dapis)


# %%
pos = 1
img_ref = dapis[pos-1]
img_cur = dapis[pos]
ref = img_ref[y_size - overlap_size:, :]
cur = img_cur[:overlap_size, :]
# ref = dapis[pos-1]
# cur = dapis[pos]

img_list = [ref,cur]

# %%
min_dist = 5
corner_list = [
    feature.corner_peaks(
        feature.corner_harris(img), 
        threshold_rel=0.001, 
        min_distance=min_dist
    )
    for img in img_list
]

img0 = img_list[0]
coords0 = corner_list[0]
matching_corners = [
    match_locations(img0, img1, coords0, coords1, min_dist)
    for img1, coords1 in zip(img_list, corner_list)
]

src = np.array(coords0)
trfm_list = [
    measure.ransac(
        (dst, src),
        transform.EuclideanTransform, 
        min_samples=3,
        residual_threshold=2, 
        max_trials=100)
        for dst in matching_corners
]
print(trfm_list)

# %%
fig, ax_list = plt.subplots(1, 2, figsize=(6, 9), sharex=True, sharey=True)
for idx, (im, trfm, (ax0, ax1)) in enumerate(zip(img_list, trfm_list, ax_list)):
    ax0.imshow(im, cmap="gray", vmin=0, vmax=1)
    ax1.imshow(transform.warp(im, trfm), cmap="gray", vmin=0, vmax=1)

    if idx == 0:
        ax0.set_title(f"Tilted images")
        ax0.set_ylabel(f"Reference Image")
        ax1.set_title(f"Registered images")

    ax0.set(xticklabels=[], yticklabels=[], xticks=[], yticks=[])
    ax1.set_axis_off()

fig.tight_layout()
# %%
