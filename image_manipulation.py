import numpy as np

from skimage.filters import gaussian, median
from skimage.morphology import square
from skimage.transform import warp, EuclideanTransform, AffineTransform
from skimage.registration import phase_cross_correlation
from skimage.util import img_as_float32
from skimage.exposure import rescale_intensity


from cellpose.transforms import resize_image
from cellpose.utils import masks_to_outlines

from scipy.ndimage import white_tophat

from tifffile import imread, TiffFile

import gc
from tqdm import tqdm
from oxdls import OMEXML


def image_crop(img, window_size):
    """
    Crops image along xy plane, starting along the middle of the plane and shifting left and down.
    """
    _, size_y, size_x = img.shape
    return img[
        :,
        int(size_y / 2 - window_size) : int(size_y / 2),
        int(size_x / 2 - window_size) : int(size_x / 2),
    ]


def image_read(filename):
    """
    Read ome-tiff file as numpy array.
    """
    return imread(filename)


def image_gaussian_filter(img, sigma=1):
    """
    Applies gaussian filter for each xy plane in the z stack and for each channel.
    """
    imgFiltered = np.zeros_like(img)
    nZ, nC, _, _ = img.shape

    for z in range(nZ):
        for c in range(nC):
            imgFiltered[z, c, :, :] = gaussian(img[z, c, :, :], sigma=sigma)

    return imgFiltered


def image_mip(img, axis=0):
    """
    Returns the Maximum Intensity Projection along a given axis.
    """
    return img.max(axis)


def image_rescale_intensity(img, out_range="dtype"):
    """
    Return image after stretching or shrinking its intensity levels.
    """
    return rescale_intensity(img, out_range=out_range).astype(np.float32)


def image_threshold(img, percentile=99):
    """
    Returns the mean intensity of pixels whose intensity is above a certain percentile.
    """
    percentile_intensity = np.percentile(img, percentile)
    img_above_percentile_intensity = img[img > percentile_intensity]
    return np.mean(img_above_percentile_intensity)


def median_filter(img):
    """
    Returns the local median of the image. 
    """
    img = img[0, ...]
    filtered = median(img, square(3))

    return filtered[None, ...]


def image_normalize_layers(img):
    imgNormalized = np.zeros_like(img)
    nZ, nC, _, _ = img.shape

    meanTotal = [
        np.mean(img[:, c, :, :].ravel()[np.nonzero(img[:, c, :, :].ravel())])
        for c in range(nC)
    ]
    stdTotal = [
        np.std(img[:, c, :, :].ravel()[np.nonzero(img[:, c, :, :].ravel())])
        for c in range(nC)
    ]

    for z in range(nZ):
        for c in range(nC):
            flattenLayer = img[z, c, :, :].ravel()
            meanLayer = np.mean(flattenLayer[np.nonzero(flattenLayer)])
            stdLayer = np.std(flattenLayer[np.nonzero(flattenLayer)])
            imgNormalized[z, c, :, :] = meanTotal[c] + (img[z, c, :, :] - meanLayer) * (
                stdTotal[c] / stdLayer
            )

    return imgNormalized


def image_downsample_shape(img, resizeFactor=0.2):
    nZ, nY, nX, nC = img.shape
    imgResized = resize_image(img, rsz=resizeFactor)

    _, nXResized, nYResized, _ = imgResized.shape
    imgResizedShaped = np.zeros((nZ, nC, nXResized, nYResized), dtype=imgResized.dtype)
    for z in range(nZ):
        for c in range(nC):
            imgResizedShaped[z, c, :, :] = imgResized[z, :, :, c]

    return imgResizedShaped


def background_subtraction(img, size=10, mode="nearest"):
    """ """
    img = img[0, ...]
    filtered = white_tophat(img, size=size, mode=mode)
    return filtered[None, ...]


def image_with_outlines(img, mask):
    """
    Given a mask, returns an image of the mask with oulines.
    """
    outlines = masks_to_outlines(mask)
    outZ, outY, outX = np.nonzero(outlines)
    imgOutlined = np.zeros((img.shape[0], img.shape[1], img.shape[2], 3))
    imgOutlined[outZ, outY, outX] = np.array([1, 1, 1])
    return imgOutlined


def image_shift(refImg, movImg):
    """
    return shift coordinates
    """

    if refImg.shape != movImg.shape and refImg.size > movImg.size:
        rz, ry, rx = refImg.shape
        mz, my, mx = movImg.shape
        movImg = np.pad(movImg, ((0, rz - mz), (0, ry - my), (0, rx - mx)), mode="mean")
        print(f"refImg.shape: {refImg.shape}, movImg.shape: {movImg.shape}")

    shift_zx, _, _ = phase_cross_correlation(
        img_as_float32(image_mip(refImg, axis=1)),
        img_as_float32(image_mip(movImg, axis=1)),
        upsample_factor=10,
    )
    shift_zy, _, _ = phase_cross_correlation(
        img_as_float32(image_mip(refImg, axis=2)),
        img_as_float32(image_mip(movImg, axis=2)),
        upsample_factor=10,
    )
    shift = np.array([(shift_zx[0] + shift_zy[0]) / 2, shift_zy[1], shift_zx[1]])

    if (np.abs(shift_zx[1]) > 50) | (np.abs(shift_zy[1]) > 50):
        shift, _, _ = phase_cross_correlation(refImg, movImg)

    return shift.astype(np.float32)


def image_warp(img, shift=None):
    """
    Warp an image according to a given coordinate transformation.
    """
    if shift is None:
        shift = np.array([0, 0, 0])

    z, y, x = img.shape
    newZ, newY, newX = np.mgrid[:z, :y, :x]
    newCoordinates = np.array([newZ - shift[0], newY - shift[1], newX - shift[2]])
    imgMoved = warp(img, newCoordinates, mode="constant", cval=0)

    del newZ, newY, newX, newCoordinates
    gc.collect()

    return imgMoved


def read_ome_metadata(filePath):
    """
    Args: file path of ome-tiff

    Returns: dictionary of parsed metadata
    """
    with TiffFile(filePath) as tif:
        imgMetadata = OMEXML(tif.ome_metadata)

    dictMetadata = {
        "nChannels": imgMetadata.image().Pixels.channel_count,
        "xPixels": imgMetadata.image().Pixels.SizeX,
        "yPixels": imgMetadata.image().Pixels.SizeY,
        "zPixels": imgMetadata.image().Pixels.SizeZ,
        "xRealSize": imgMetadata.image().Pixels.PhysicalSizeX,
        "yRealSize": imgMetadata.image().Pixels.PhysicalSizeY,
        "zRealSize": imgMetadata.image().Pixels.PhysicalSizeZ,
        "channelNames": [
            imgMetadata.image().Pixels.Channel(ch).Name.split("_")[0]
            for ch in range(imgMetadata.image().Pixels.channel_count)
        ],
    }

    return dictMetadata
