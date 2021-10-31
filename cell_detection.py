from cellpose import models, utils
from image_manipulation import *
import numpy as np
from skimage import img_as_float32


def run_cellpose(
    img, 
    channels, 
    gpu=True, 
    model_type='cyto', 
    diameter=30, 
    device=None, 
    anisotropy=None,
    stitch_threshold=0.25,
    cellprob_threshold=0.0,
    flow_threshold=0.4,
    min_size=15,
    do_3D=False):
    # model = models.Cellpose(gpu=gpu, model_type=model, device=device, torch=False)
    model = models.Cellpose(gpu=gpu, model_type=model_type)
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

def cell_detection(img, zToXYRatioReal=1, resizeFactor=0.2, color_shift=None):
    """
    main function for cell detection using cellpose
    """
    
    zToXYRatio = zToXYRatioReal * resizeFactor
    img_cells = img_as_float32(img)
   
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
    maskUpsampled = mask_upsample(maskClosed, (img.shape[0], img.shape[1], img.shape[2]))
    maskUpsampledOutlined = utils.masks_to_outlines(maskUpsampled)
    # maskUpsampledOutlined = image_with_outlines(img[0], maskUpsampled)

    # del imgResizedShaped, imgNormalized, imgFiltered, mask, maskClosed
    # gc.collect()

    return maskUpsampled, maskUpsampledOutlined, zToXYRatioReal, nCells


