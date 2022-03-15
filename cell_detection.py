from image_manipulation import *
from cellpose.models import Cellpose

# import numpy as np
from skimage import img_as_float32


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
    # model = models.Cellpose(gpu=gpu, model_type=model, device=device, torch=False)
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
        # z_axis=0,
        # channel_axis=1
    )
    return mask, flow, style, diam


def cell_detection(img, zToXYRatioReal=1, resizeFactor=0.2):
    """
    main function for cell detection using cellpose
    """

    zToXYRatio = zToXYRatioReal * resizeFactor
    img_cells = img_as_float32(img)

    imgResizedShaped = image_downsample_shape(img_cells, resizeFactor=resizeFactor)
    imgNormalized = image_normalize_layers(imgResizedShaped)
    imgFiltered = image_gaussian_filter(imgNormalized, sigma=2)

    medianDiameters = np.float32(30.0)

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

    maskClosed = mask_closing(mask.copy())
    maskUpsampled = mask_upsample(
        maskClosed, (img.shape[0], img.shape[1], img.shape[2])
    )

    maskExpanded = maskUpsampled

    del imgResizedShaped, imgNormalized, imgFiltered, mask, maskClosed, maskUpsampled
    gc.collect()

    return maskExpanded, zToXYRatioReal, nCells
