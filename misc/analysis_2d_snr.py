from math import sqrt
import numpy as np
import pandas as pd

from skimage import io, filters, morphology, measure, color, util, segmentation, transform, restoration, exposure
from skimage.feature import blob_log
from skimage.draw import circle
# from skimage.external import tifffile
from scipy.ndimage import convolve

import matplotlib.pyplot as plt
import napari

import sys
import os
from tqdm import tqdm
from numba import njit
from math import sqrt


## $ python3 analysis_2d.py (filepath)/*.tif


file = sys.argv[1]
path = file.split('/')[0]
name = file.split('/')[-1]

for root, dirs, fname in os.walk(os.path.join(path)):
    fnames = fname
files = [os.path.join(path, fname) for fname in fnames if fname[0] is not '.']
imgs = [io.imread(file) for file in files]
print('-- images loaded')

if len(files) > 1:
    img = imgs[1]
    img_40 = imgs[0]
else:
    img = imgs[0]
    img_40 = img

bg = morphology.area_opening(img_40, area_threshold=2000)
fg = img_40 - bg
ero = morphology.erosion(fg)
ero = morphology.erosion(ero)
mask_fg = morphology.erosion(ero) > 0
mask_bg = util.invert(mask_fg)

img_fg = img * mask_fg.astype(np.int)
img_bg = img * mask_bg.astype(np.int)

fginten = img_fg.ravel()
fginten = fginten[np.nonzero(fginten)]
bgmean = img_bg.ravel()
bgmean = bgmean[np.nonzero(bgmean)].mean()
print(bgmean)
result_df = pd.DataFrame({files[0] : fginten.flatten()})
result_df.to_excel(excel_writer=files[0] + "_result.xlsx")
result_snr = pd.DataFrame({files[0] : fginten.flatten()/bgmean})
result_snr.to_excel(excel_writer=files[0] + "_result_snr.xlsx")

with napari.gui_qt():
    viewer = napari.Viewer()
    
    viewer.add_image(img, name='original', is_pyramid=False)
    viewer.add_image(img_40, name='40perc laser', is_pyramid=False)
    viewer.add_image(mask_fg, name='mask fg', is_pyramid=False,)
    viewer.add_image(mask_bg, name='mask bg', is_pyramid=False)

