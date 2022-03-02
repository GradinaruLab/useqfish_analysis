import numpy as np
import zarr
import os

path = './result/result_images.zarr'
filename = 'position6.zarr'

cells = zarr.load(path)
cellLabels = cells['cellLabels']

print(cellLabels.shape)

zarr.save(
    os.path.join('./expression_matrices/211112/cell_labels', filename),
    cellLabels
)

# np.save(
#     os.path.join('./expression_matrices/220116/cell_labels', filename),
#     cellLabels
# )

# read = np.load(os.path.join('./expression_matrices/220116/cell_labels', filename))
# print(read.shape)