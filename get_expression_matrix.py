# Copyright 2022 California Institute of Technology
from tkinter import N
from webbrowser import get
from cell_detection import *
from spot_detection import *
from image_manipulation import *
import params
from params import *

import pandas as pd

from argparse import ArgumentParser
import os
from glob import glob
import logging

from warnings import filterwarnings

filterwarnings("ignore")


def get_expression_matrix(path):
    # start log
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    fh = logging.FileHandler(f"{path}/get_expression_matrix.log")
    fh.setLevel(logging.DEBUG)
    # create formatter
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    # add formatter to ch
    fh.setFormatter(formatter)
    # add ch to logger
    logger.addHandler(fh)

    # get filenames
    filepath = os.path.join(path, "*.tif")
    filenames = sorted(glob(filepath), key=os.path.basename)

    # print paramters to log file
    parameters = globals().get("params", None)
    if parameters:
        [
            logger.info(f"{key} {value}")
            for key, value in parameters.__dict__.items()
            if not (key.startswith("__") or key.startswith("_"))
        ]
    print(filenames)

    nR = len(filenames)  # number of rounds

    imgReference = image_read(filenames[roundRef])

    imgMetadata = read_ome_metadata(filenames[roundRef])

    zToXYRatioReal = imgMetadata["PhysicalSizeZ"] / imgMetadata["PhysicalSizeX"]

    nC, size_z, size_y, size_x = imgReference.shape

    print(f">> STEP 1. Cell detection -")
    logger.info(">> STEP 1. Cell detection -")
    img_cells = img_as_float32(
        np.stack((imgReference[cellch], imgReference[0]), axis=3)
    )
    cellLabels, zToXYRatioReal, nCells = cell_detection(
        img_cells, zToXYRatioReal=zToXYRatioReal, resizeFactor=0.2
    )
    logger.info(f"{nCells} cells detected")

    dapi_reference_cropped = image_crop(imgReference[0], shift_window_size)

    spots_assigned_allrounds = []
    shifts_allrounds = []
    dapis_shifted = []

    logger.debug(">> STEPS 2 and 3: \n Registration and Spot Detection for each round")
    for filename in filenames[:roundRef]:

        print("\n")
        print(filename)
        logger.debug(f"Analyzing file {filename}")

        img = dask.delayed(image_read)(filename)

        print(f">> STEP 2. registration - ")

        dapi = img_as_float32(img[0].compute())
        dapi_cropped = image_crop(dapi, shift_window_size)
        round_shift = image_shift(dapi_reference_cropped, dapi_cropped)

        shifts = [
            np.array(round_shift) + np.array(color_shift)
            for color_shift in color_shifts
        ]
        print(shifts)
        logger.debug(f"shits are {shifts}")

        shifts_allrounds.append(np.array(shifts).astype(np.float32))
        dapis_shifted.append(
            image_warp(image_crop(dapi_cropped, 500), shift=np.array(shifts[0]))
        )

        print(f">> STEP 3. Spot detection -")

        spots, spots_assigned = spot_detection(
            img, absolute_thresholds, cellLabels, shifts=shifts, shape=dapi.shape
        )
        spots_assigned_allrounds.append(spots_assigned)

        print(
            f"# of spots detected for this round: {[len(spots) for spots in spots_assigned]}"
        )
        logger.debug(
            f"{[len(spots) for spots in spots_assigned]} spots detected in this round"
        )

    dapis_shifted.append(
        image_warp(image_crop(dapi_reference_cropped, 500), shift=color_shifts[0])
    )

    print(f">> STEP 4. Save results -")
    logger.debug(">> STEP 4. Save results -")
    save_to_zarr(
        os.path.join(path, "result/result_images.zarr"),
        imgCells=img_cells,
        cellLabels=cellLabels,
        zToXYRatioReal=zToXYRatioReal,
        shifts_allrounds=np.array(shifts_allrounds),
        dapis_shifted=np.array(dapis_shifted),
        nR=nR,
        thresholds=np.array(absolute_thresholds),
    )

    spots_results = []
    for r, spots_assigned in enumerate(spots_assigned_allrounds):
        for c, spots in enumerate(spots_assigned):
            for spot in spots:
                spots_results.append(np.append(np.array([r + 1, c + 1]), spot))

    spots_results = np.array(spots_results)
    print(f">>>> Total {spots_results.shape[0]} spots detected")
    print(f">>>> intensity threshold: {absolute_thresholds}")

    resultDf = pd.DataFrame(
        {
            "round": spots_results[:, 0],
            "channel": spots_results[:, 1],
            "z-coord": spots_results[:, 2],
            "y-coord": spots_results[:, 3],
            "x-coord": spots_results[:, 4],
            "cell id": spots_results[:, 5],
        }
    )
    resultDf.to_excel(excel_writer=os.path.join(path, "result/result.xlsx"))
    logger.info(f"Cell Detection Results:\n {resultDf.to_string()}")

    target_index = np.zeros((nR, nC), dtype=np.int)
    index = 0
    for r in range(nR - 1):
        for c in range(nC - 1):
            target_index[r, c] = index
            index = index + 1

    spots_per_cell = np.zeros((nCells + 1, (nR - 1) * (nC - 1)), dtype=np.int)
    for spot_index in range(spots_results.shape[0]):
        cell_id = int(spots_results[spot_index, -1])
        r = spots_results[spot_index, 0]
        c = spots_results[spot_index, 1]
        spots_per_cell[cell_id, target_index[r - 1, c - 1]] = (
            spots_per_cell[cell_id, target_index[r - 1, c - 1]] + 1
        )

    resultDf2 = pd.DataFrame(data=spots_per_cell)
    resultDf2.to_excel(
        excel_writer=os.path.join(path, "result/result_spots_per_cell.xlsx")
    )
    logger.info(f"Results Matrix:\n {resultDf2.to_string()}")


if __name__ == "__main__":
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

    get_expression_matrix(path)
