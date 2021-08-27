from skimage import filters, feature, img_as_float32, registration, transform, morphology, segmentation
import numpy as np
from cellpose import transforms, utils
from scipy.ndimage import white_tophat
import tifffile


def image_read(filename):
    """
    read ome-tiff file as numpy array, converted to [0.0 1.0] float32 type
    if dapi channel is not the first, bring the dapi channel to the first
    """
    return tifffile.imread(filename)


def image_gaussian_filter(img, sigma=1):
    imgFiltered = np.zeros_like(img)
    nZ, nC, _, _ = img.shape

    for z in range(nZ):
        for c in range(nC):
            imgFiltered[z,c,:,:] = filters.gaussian(img[z,c,:,:], sigma=sigma)
    
    return imgFiltered


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


def mask_upsample(mask, finalShape=None):
    if finalShape is None:
        raise ValueError('please provide the final shape of the mask')
    elif mask.ndim != len(finalShape):
        raise ValueError('dimension of final shape is not matched with input mask')
    else:
        maskUpsampled = transform.resize(mask, finalShape, order=0, preserve_range=True)
        
        return maskUpsampled.astype(np.int)


def image_with_outlines(img, mask):
    outlines = utils.masks_to_outlines(mask)
    outZ, outY, outX = np.nonzero(outlines)
    imgOutlined = np.zeros((img.shape[0], img.shape[1], img.shape[2], 3))
    imgOutlined[outZ, outY, outX] = np.array([1, 1, 1])
    return imgOutlined


def image_with_outlines(img, mask):
    outlines = utils.masks_to_outlines(mask)
    outZ, outY, outX = np.nonzero(outlines)
    imgOutlined = np.zeros((img.shape[0], img.shape[1], img.shape[2], 3))
    imgOutlined[outZ, outY, outX] = np.array([1, 1, 1])
    return imgOutlined


def image_shift(movImg, refImgDapi):
    """
    return shift coordinates
    """
    # refImg = refImg[0,...]
    shift, _, _ = registration.phase_cross_correlation(refImgDapi, movImg)

    # return shift[None,...]
    return shift