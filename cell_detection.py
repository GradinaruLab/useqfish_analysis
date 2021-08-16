from cellpose import models, utils
from image_manipulation import image_downsample_shape, image_gaussian_filter, image_normalize_layers, mask_upsample, image_with_outlines
from tqdm import tqdm
import numpy as np
from skimage import morphology, transform, img_as_float32, segmentation

import gc
import tifffile
from oxdls import OMEXML

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