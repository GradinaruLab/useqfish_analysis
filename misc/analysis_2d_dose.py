from math import sqrt
import numpy as np
import pandas as pd

from skimage import io, filters, morphology, measure, color, util, segmentation
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

def blobdetection(img, threshold):
    # print('-- blob detection')
    img_med = filters.median(img, morphology.disk(3))
    img_bg = img_med - morphology.area_opening(img_med)
    img_spot = np.zeros_like(img)
    spots = blob_log(img_bg, max_sigma=10, threshold=threshold, overlap=0).astype(np.uint16)
    # spots[:,2] = spots[:,2] * sqrt(2)
    img_spot[spots[:,0], spots[:,1]] = 1
    img_spot = convolve(img_spot, morphology.disk(5), mode='constant')
    return spots, img_spot

# @njit
def spotassigment(spots, dapi_rps, img, disk):
    # print('-- spot assigning to cells')
    spotnum = spots.shape[0]
    spotid = np.zeros((spotnum,))
    spot_label = np.zeros_like(img)
    disk_r = int(np.floor(disk.shape[0]/2))
    for i in range(spotnum):
        min_dists = np.array([distance(spots[i,:2], rp).min() for rp in dapi_rps])
        if np.min(min_dists) < 100:
            spotid[i] = np.argmin(min_dists)+1
            spot_y, spot_x = spot_label[spots[i,0]-disk_r:spots[i,0]+disk_r+1, spots[i,1]-disk_r:spots[i,1]+disk_r+1].shape
            spot_label[spots[i,0]-disk_r:spots[i,0]+disk_r+1, spots[i,1]-disk_r:spots[i,1]+disk_r+1] = disk[:spot_y, :spot_x]*spotid[i]
    # return spotid, spot_label
    return spot_label, spotid
    # return spot_label  

@njit
def distance(pt, refpts):
    # pt = pt.astype(np.float32)
    # refpts = refpts.astype(np.float32)
    dists = np.zeros((refpts.shape[0],))
    for i in range(refpts.shape[0]):
        dists[i] = sqrt((pt[0]-refpts[i,0])**2+(pt[1]-refpts[i,1])**2)
    return dists

def celloutline(merge_label):
    rps = measure.regionprops(merge_label)
    celloutlined = np.zeros_like(merge_label)
    for rp in rps:
        min_row, min_col, max_row, max_col = rp.bbox
        convex_img = rp.convex_image
        celloutlined[min_row:max_row, min_col:max_col] = convex_img.astype(np.int)*rp.label
    return celloutlined

@njit
def spotcount(num_cells, spotid):
    counts = np.zeros((num_cells,))
    for i in range(spotid.size):
        if spotid[i] > 0:
            counts[int(spotid[i]-1)] += 1
    return counts

file = sys.argv[1]
path = file.split('/')[0]
name = file.split('/')[-1]

for root, dirs, fname in os.walk(path):
    fnames = fname
files = [os.path.join(path, fname) for fname in fnames if fname[0] is not '.']
# files = files[1:]
# print(files)
imgs = [io.imread(file) for file in files]
channels = [list(file.split('.')[-2])[-1] for file in files]

imgs_rna = []
channels_rna = []
dapi_ind = 0
for c, (img, ch) in enumerate(zip(imgs, channels)):
    if ch is 'i':
        dapi_ind = c
    if (ch is not 'i') and (ch is not '0'):
        imgs_rna.append(img)
        channels_rna.append(int(ch))
ch_rna_int = np.array(channels_rna).argsort()


imgs_rna_cp = []
for cri in list(ch_rna_int):
    imgs_rna_cp.append(imgs_rna[cri])

imgs_rna = imgs_rna_cp
channels_rna = list(np.array(channels_rna)[ch_rna_int])

# for img, ch in zip(imgs[:-1], channels[:-1]):
#     if int(ch) in range(1,5):
#         imgs_rna.append(img)
#         channels_rna.append(ch)

print('-- images loaded')

## dapi 
dapi = imgs[dapi_ind]
dapi_mask = util.invert(util.img_as_bool(dapi))
dapi_label = measure.label(dapi_mask)
dapi_labelrgb = color.label2rgb(dapi_label, dapi, bg_label=0)

## blob detection and assigning to each cell
spots = []
imgs_spots = []
spots_label = []
spots_id = []

num_cells = dapi_label.max()
dapi_rps = measure.regionprops(dapi_label)
dapi_coords = [dapi_rp.coords for dapi_rp in dapi_rps]

spots_counts = []
thresholds = [0.02, 0.01, 0.04, 0.02]
result = {'cell id': range(1,num_cells+1)}
for img, threshold, ch in tqdm(zip(imgs_rna, thresholds, channels_rna)):
    spot, img_spot = blobdetection(img, threshold)
    spot_label, spotid = spotassigment(spot, dapi_coords, dapi_label, morphology.disk(3))
    spot_labelrgb = color.label2rgb(spot_label, bg_label=0)
    
    counts = spotcount(num_cells, spotid)

    spots.append(spot)
    imgs_spots.append(img_spot)
    spots_label.append(spot_label)
    spots_id.append(spotid)
    spots_counts.append(counts)
    result.update({ch : counts})

spots_label.append(dapi_label)
spots_merge = np.array(spots_label).max(axis=0)
spots_merge_rgb = color.label2rgb(spots_merge, bg_label=0)

celloutlined = celloutline(spots_merge)
merge_bndry = segmentation.mark_boundaries(spots_merge_rgb, celloutlined)
io.imsave("data/segmented.jpg", merge_bndry)

result_df = pd.DataFrame(result)
result_df.to_excel(excel_writer = "data/results.xlsx")

with napari.gui_qt():
    viewer = napari.Viewer()
    viewer.add_image(dapi, name='dapi')
    viewer.add_image(dapi_label, name='dapi detected')
    for img_spot, img_rna, ch in zip(imgs_spots, imgs_rna, channels_rna):
        viewer.add_image(img_rna, name='ch'+str(ch))
        viewer.add_image(img_spot, name='ch'+str(ch)+' detected')
        # viewer.add_image(spot_label)
    # viewer.add_image(spots_merge_rgb)
    viewer.add_image(merge_bndry)