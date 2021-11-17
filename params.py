nC = 5      # number of channels / colors
nR = 5
# nR = 11     # number of rounds (except the last dT labeling)
roundRef = -1    #  round containing the cell reference; use -1 for the last round
cellch = 3      # channel containing cell reference 
img_size = (2048, 2048)

sigma = 2      # estimate size of spots to detect

n_variants = 6  # number of variants in pool

# shift_window_size = 200     # ? TODO
shift_window_size = 1000

color_shifts = [[8, 0, 0],          # color aberration shifts for each channel to 647 aquired by imaging focalcheck #1, 500nm
                [-1, 0, 1],         # from 210507_focalcheck
                [0, 0, 0],
                [0, 0, 0],
                [1, -3, 2]]

# color_shifts = [[ 8, 0, 0],      # from 211031_focalcheck
#                 [-1, 1, 0],
#                 [-1, 1, 0],
#                 [ 0, 0, 0],
#                 [ 1, -5, 2]]

stitching_shape = (3, 5)                        # montage shape (w, h)
stitching_overlap = 1                           # % overlap between tiles
stitching_size = (4550, 9719)                   # the pixel size of stiched image in fusion/imaris
stitching_coords = [[0.795, -20.030, 0.000],    # stitching coordinates of each position
                    [1.986, -15.122, 0.000],    # calculated with FUSION stitching
                    [3.177, -10.215, 0.000],    # stored in XX_global.xml
                    [4.368, -5.307, 0.000],     # [x, y, z]
                    [5.367, -1.189, 0.000],
                    [0, 0, 0],
                    [-1.041, -4.117, 0.078],
                    [-2.190, -9.026, 0.000],
                    [-3.381, -13.933, 0.000],
                    [-4.572, -18.841, 0.000],
                    [-5.873, -18.359, -0.157],
                    [-4.656, -13.859, -0.157],
                    [-3.481, -8.972, -0.137],
                    [-2.343, -4.021, -0.146],
                    [-1.274, 0.104, -0.250]]

# gene_list = ['PHP.N', 'PHP.eB', 'PHP.B8', None,             # 210430
#              'PHP.Astro', 'CAP-B10', 'PHP.V1', None,
#              'gad1', 'slc17a7', 'mbp', 'gja1',
#              'cldn5', 'hexb', 'acta2', 'mrc1',
#              'sst', 'pvalb', None, 'vip']

gene_list = ['PHP.N', 'PHP.eB', None, 'PHP.B8',             # 210828
             'PHP.Astro', 'CAP-B10', None, 'PHP.V1',
             'gad1', 'slc17a7', 'gja1', 'mbp',
             'cldn5', 'hexb', 'mrc1', 'acta2',
             'sst', 'pvalb', 'vip', None]

gene_list_ordered = [gene for gene in gene_list if gene is not None]

# gene_list = ['PHP.N', 'PHP.eB', None, 'PHP.B8',
#              'PHP.Astro', 'CAP-B10', None, 'PHP.V1',
#              'cux2', 'slc30a3', 'osr1', 'rorb',
#             #  'trh', 'chrna6', 'hsd11b1', 'foxp2',
#              'trh', 'chrna6', 'hsd11b1', None,          # foxp2, tshz2 excluded due to (putative) non-specific binding
#             #  'rprm', 'tshz2', 'crym', 'gad1',
#              'rprm', None, 'crym', 'gad1',
#              'sncg', 'lamp5', 'serpinf7', 'krt73',
#              'tac1', 'chodl', 'crh', 'th',
#              'calb2', 'gpc3', 'calb1', 'tpbg',
#              'sst', 'reln', 'pvalb', 'vip',
#              'opalin', 'ccnb1', 'enpp6', 'dcn',
#              'pecam1', 'kcnj8', 'mrc1', 'hexb']

# gene_list_ordered = ['PHP.eB', 'CAP-B10', 'PHP.N', 
#                      'PHP.Astro', 'PHP.V1', 'PHP.B8',
#                      'slc30a3', 'cux2', 'osr1', 'rorb',
#                      'chrna6', 'trh', 'hsd11b1',
#                      'rprm', 'crym', 'gad1', 'lamp5',
#                      'sncg', 'serpinf7', 'krt73', 'chodl',
#                      'tac1', 'crh', 'th', 'gpc3', 
#                      'calb2', 'calb1', 'tpbg', 'reln',
#                      'sst', 'pvalb', 'vip',
#                      'ccnb1', 'opalin', 'enpp6',
#                      'dcn',
#                      'kcnj8',
#                      'pecam1',
#                      'mrc1',
#                      'hexb']

# marker_genes_dict = {
#     'Neuronal': ['slc30a3'],
#     'Layer 2': ['cux2'],
#     'Layer 3': ['cux2'],
#     'Layer 4': ['rorb', 'hsd11b1'],
#     'Layer 5': ['rorb', 'chrna6', 'hsd11b1', 'rprm', 'crym'],
#     'Layer 6': ['trh', 'rprm'],
#     'Inhibitory': ['gad1'],
# }