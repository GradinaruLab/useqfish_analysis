# unmute when the error "libgcc_s.so.1 must be installed for pthread_cancel to work" came up
# import ctypes
# libgcc_s = ctypes.CDLL('libgcc_s.so.1')

from tkinter import N
from cell_detection import *
from spot_detection import *
from image_manipulation import *
from params import *

import pandas as pd

# from pylab import *

from argparse import ArgumentParser
import os
from glob import glob

import dask
import dask.array as da
from dask.diagnostics import ProgressBar
from dask_image import ndfilters
import zarr

from warnings import filterwarnings

filterwarnings("ignore")

my_parser = ArgumentParser(description="Run analysis on a group of images")

my_parser.add_argument(
    "Path",
    metavar="path",
    type=str,
    help="the path to the diretory containing the images",
)


# Parse arguments
args = my_parser.parse_args()

path = args.Path
filepath = os.path.join(path, "*.tif")
filenames = sorted(glob(filepath), key=os.path.basename)
print(filenames)
nR = len(filenames)  # number of rounds

imgReference = image_read(filenames[roundRef])

imgMetadata = read_ome_metadata(filenames[roundRef])
zToXYRatioReal = imgMetadata["PhysicalSizeZ"] / imgMetadata["PhysicalSizeX"]
nC, size_z, size_y, size_x = imgReference.shape



print(f">> STEP 1. Cell detection -")
img_cells = img_as_float32(np.stack((imgReference[cellch], imgReference[0]), axis=3))
cellLabels, zToXYRatioReal, nCells = cell_detection(
    img_cells, zToXYRatioReal=zToXYRatioReal, resizeFactor=0.2
)


save_to_zarr(
    os.path.join(path, "result/result_images.zarr"),
    cellLabels=cellLabels,
)
