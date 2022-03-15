nC = 5  # number of channels / colors
# nR = 5
# nR = 13     # number of rounds (except the last dT labeling)
# nR = 10
roundRef = -1  #  round containing the cell reference; use -1 for the last round
cellch = 3  # channel containing cell reference
# cellch = 1

# nC = 3          # smHCR (211112)
# nR = 1
# roundRef = -1
# cellch = 2
img_size = (2048, 2048)

sigma = [0.4, 0.7, 1.0]  # estimate size of spots to detect

# shift_window_size = 200     # ? TODO
shift_window_size = 1000
# absolute_thresholds = [300., 500., 200., 500.]        # 210430
# absolute_thresholds = [350., 450., 250., 350.]    # 210828
# absolute_thresholds = [250., 700., 300., 500.]    # 211229
# absolute_thresholds = [400., 700., 250., 600.]      # 220116
# absolute_thresholds = [500., 600., 200., 400.]      # 211229 - cortex
# absolute_thresholds = [500., 700., 400., 500.]        # 211229 - striatum
# absolute_thresholds = [150., 200., 150., 200.]          # 211229 - thalamus
# absolute_thresholds = [200., 300., 150., 300.]          # 211229 - cerebellum
# absolute_thresholds = [200., 600.]                          # 211112 - smhcr
absolute_thresholds = [200.0, 600.0, 250.0, 600.0]  # 211024

virus_list = ["PHP.eB", "CAP-B10", "PHP.N", "PHP.Astro", "PHP.V1", "PHP.B8"]
n_variants = len(virus_list)  # number of variants in pool

# leiden_resolution = .5  # 210430, 210828
leiden_resolution = 0.5

# color_shifts = [[8, 0, 0],          # color aberration shifts for each channel to 647 aquired by imaging focalcheck #1, 500nm
#                 [-1, 0, 1],         # from 210507_focalcheck
#                 [0, 0, 0],
#                 [0, 0, 0],
#                 [1, -3, 2]]

color_shifts = [
    [8, 0, 0],  # from 211031_focalcheck
    [-1, 1, 0],
    [-1, 1, 0],
    [0, 0, 0],
    [1, -5, 2],
]

# color_shifts = [[ 8, 0, 0],      # smHCR (211112)
#                 [-1, 1, 0],
#                 [ 0, 0, 0]]

# stitching_shape = (3, 5)                        # montage shape (w, h)
# stitching_overlap = 1                           # % overlap between tiles
# stitching_size = (4550, 9719)                   # the pixel size of stiched image in fusion/imaris
# stitching_coords = [[0.795, -20.030, 0.000],    # stitching coordinates of each position
#                     [1.986, -15.122, 0.000],    # calculated with FUSION stitching
#                     [3.177, -10.215, 0.000],    # stored in XX_global.xml
#                     [4.368, -5.307, 0.000],     # [x, y, z]
#                     [5.367, -1.189, 0.000],
#                     [0, 0, 0],
#                     [-1.041, -4.117, 0.078],
#                     [-2.190, -9.026, 0.000],
#                     [-3.381, -13.933, 0.000],
#                     [-4.572, -18.841, 0.000],
#                     [-5.873, -18.359, -0.157],
#                     [-4.656, -13.859, -0.157],
#                     [-3.481, -8.972, -0.137],
#                     [-2.343, -4.021, -0.146],
#                     [-1.274, 0.104, -0.250]]

######################
# stitching parameters - 220116
# stitching_shape = (4, 6)                        # montage shape (w, h)
stitching_shape = (6, 4)  # (h, w)
stitching_overlap = 10  # % overlap between tiles
# stitching_size = (7530, 11134)
stitching_size = (11134, 7530)
stitching_coords = [
    [-7.615, -22.138, 0.261],
    [-6.120, -17.701, 0.125],
    [-4.598, -13.090, 0.071],
    [-3.059, -8.732, 0.107],
    [-1.517, -4.413, 0.049],
    [0.0, 0.0, 0.0],
    [-4.876, 1.226, -0.081],
    [-6.135, -3.616, -0.117],
    [-7.392, -8.382, -0.220],
    [-9.012, -12.471, 0.129],
    [-10.267, -17.106, 0.000],
    [-12.130, -22.385, -0.612],
    [-16.937, -21.038, -0.701],
    [-15.563, -16.348, -0.694],
    [-14.174, -11.718, -0.563],
    [-12.777, -7.307, -0.491],
    [-11.356, -2.713, -0.533],
    [-9.933, 1.870, -0.534],
    [-14.822, 3.189, -0.544],
    [-16.165, -1.452, -0.559],
    [-17.562, -6.150, -0.524],
    [-19.007, -10.466, -0.621],
    [-20.426, -15.106, -0.790],
    [-21.791, -19.762, -0.847],
]


# gene_list = ['PHP.N', 'PHP.eB', 'PHP.B8', None,             # 210430
#              'PHP.Astro', 'CAP-B10', 'PHP.V1', None,
#              'gad1', 'slc17a7', 'mbp', 'gja1',
#              'cldn5', 'hexb', 'acta2', 'msr1',
#              'sst', 'pvalb', None, 'vip']

# gene_list = ['PHP.N', 'PHP.eB', None, 'PHP.B8',             # 210828
#              'PHP.Astro', 'CAP-B10', None, 'PHP.V1',
#              'gad1', 'slc17a7', 'gja1', 'mbp',
#              'cldn5', 'hexb', 'msr1', 'acta2',
#              'sst', 'pvalb', 'vip', None]

# gene_list_ordered = [gene for gene in gene_list if gene is not None]

# gene_list = ['PHP.N', 'PHP.eB', None, 'PHP.B8',         # 211024
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
#                      'slc30a3',
#                      'cux2', 'rorb', 'hsd11b1', 'chrna6', 'rprm', 'crym', 'osr1', 'trh',
#                      'gad1',
#                      'lamp5', 'sncg', 'krt73',
#                      'pvalb', 'tpbg', 'sst', 'chodl', 'crh', 'reln',
#                      'vip', 'gpc3']

# gene_list = ['PHP.N', 'PHP.eB', None, 'PHP.B8',         # 211202
#              'PHP.Astro', 'CAP-B10', None, 'PHP.V1',
#              'gad1', 'slc17a7', 'gja1', 'mbp',
#              'cldn5', 'hexb', 'msr1', 'acta2',
#              'sst', 'pvalb', 'vip', 'slc6a4',
#              'chat', 'th', None, None,
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
#              'pecam1', 'kcnj8', 'mrc1', 'hexb',
#              'aqp4', 'gja1', 'gad2', 'slc17a7',
#              'ly6a', 'pcp4', 'ly6c1', 'car4',
#              'gfap', 'mbp', None, None]

# gene_list = ['PHP.N', 'PHP.eB', None, 'PHP.B8',         # 211202
#              'PHP.Astro', 'CAP-B10', None, 'PHP.V1',
#             #  'gad1', 'slc17a7', 'gja1', 'mbp',
#              None, 'slc17a7', 'gja1', 'mbp',
#              'cldn5', 'hexb', 'msr1', 'acta2',
#              'sst', 'pvalb', 'vip', 'slc6a4',
#              'chat', 'th', None, None,
#              'cux2', 'slc30a3', 'osr1', 'rorb',
#             #  'trh', 'chrna6', 'hsd11b1', 'foxp2',
#              'trh', 'chrna6', 'hsd11b1', None,          # foxp2, tshz2 excluded due to (putative) non-specific binding
#             #  'rprm', 'tshz2', 'crym', 'gad1',
#              'rprm', None, 'crym', 'gad1',
#              'sncg', 'lamp5', 'serpinf1', 'krt73',
#             #  'tac1', 'chodl', 'crh', 'th',
#              'tac1', 'chodl', 'crh', None,
#              'calb2', 'gpc3', 'calb1', 'tpbg',
#             #  'sst', 'reln', 'pvalb', 'vip',
#              None, 'reln', None, None,
#              'opalin', 'ccnb1', 'enpp6', 'dcn',
#             #  'pecam1', 'kcnj8', 'mrc1', 'hexb',
#              'pecam1', 'kcnj8', 'mrc1', None,
#             #  'aqp4', 'gja1', 'gad2', 'slc17a7',
#              'aqp4', None, 'gad2', None,
#              'ly6a', 'pcp4', 'ly6c1', 'car4',
#             #  'gfap', 'mbp', None, None]
#             'gfap', None, None, None]

# gene_list_ordered = ['PHP.eB', 'CAP-B10', 'PHP.N',
#                      'PHP.Astro', 'PHP.V1', 'PHP.B8',
#                      'slc30a3', 'slc17a7', 'pcp4',
#                      'cux2', 'rorb', 'hsd11b1', 'chrna6', 'rprm', 'crym', 'osr1', 'trh',
#                      'gad1', 'gad2',
#                      'lamp5', 'sncg', 'krt73',
#                      'pvalb', 'tpbg', 'sst', 'chodl', 'crh', 'reln',
#                      'vip', 'gpc3',
#                      'gja1', 'serpinf1', 'gfap', 'aqp4',
#                      'mbp', 'ccnb1', 'opalin', 'enpp6',
#                      'dcn', 'acta2',
#                      'kcnj8', 'msr1', 'mrc1',
#                      'pecam1', 'cldn5',
#                      'hexb',
#                      'th', 'chat', 'slc6a4',
#                      'calb1', 'calb2',
#                      'ly6a', 'ly6c1', 'car4']

# gene_list_ordered = ['PHP.eB', 'CAP-B10', 'PHP.N',                  # striatum
#                      'PHP.Astro', 'PHP.V1', 'PHP.B8',
#                      'tac1', 'gad2', 'trh', 'tpbg', 'sst', 'reln', 'calb1', 'crym', 'gad1', 'chat', 'pcp4']

# gene_list_ordered = ['PHP.eB', 'CAP-B10', 'PHP.N',                  # midbrain
#                      'PHP.Astro', 'PHP.V1', 'PHP.B8',
#                      'calb2', 'tac1', 'pvalb', 'gad2', 'th', 'slc6a4', 'trh', 'sst', 'chrna6', 'calb1', 'gpc3', 'chat', 'sncg', 'rprm', 'gad1', 'pcp4']

# gene_list_ordered = ['PHP.eB', 'CAP-B10', 'PHP.N',                  # cerebellum
#                      'PHP.Astro', 'PHP.V1', 'PHP.B8',
#                      'pvalb', 'gad2', 'enpp6', 'slc30a3', 'hsd11b1', 'sst', 'calb1', 'gad1', 'pcp4']

# gene_list = ['PHP.N', 'PHP.eB', 'slc30a3', 'PHP.B8',                  # 211229
#              'PHP.Astro', 'CAP-B10', 'gad2', 'PHP.V1',
#              'cux2', 'slc17a7', 'ctgf', 'sulf2',
#              'foxp1', 'foxp2', 'calb2', 'calb1',
#              'sst', 'pvalb', 'vip', 'drd1',
#              'tac1', 'lamp5', 'drd2', 'gad1',
#              'ly6a', 'pcp4', 'crym', 'rorb',
#             #  None, 'pcp4', 'crym', 'rorb',
#              'rprm', 'reln', 'hsd11b1', 'tpbg',
#              'necab1', 'tnnt1', 'crh', 'prkcd',
#              'gdf10', 'lgi2', 'ppp1r17', 'gabra6',
#              'trh', 'chrna6', 'enpp6', 'th',
#              'sncg', 'chodl', 'serpinf1', 'krt73',
#              'gpc3', 'ccnb1', None, None]

# gene_list_ordered = ['PHP.eB', 'CAP-B10', 'PHP.N', 'PHP.Astro', 'PHP.V1', 'PHP.B8',
#                      'slc30a3', 'gad2',
#                      'cux2', 'slc17a7', 'ctgf', 'sulf2',
#                      'foxp1', 'foxp2', 'calb2', 'calb1',
#                      'sst', 'pvalb', 'vip', 'drd1',
#                      'tac1', 'lamp5', 'drd2', 'gad1',
#                      'ly6a', 'pcp4', 'crym', 'rorb',
#                      'rprm', 'reln', 'hsd11b1', 'tpbg',
#                      'necab1', 'tnnt1', 'crh', 'prkcd',
#                      'gdf10', 'lgi2', 'ppp1r17', 'gabra6',
#                      'trh', 'chrna6', 'enpp6', 'th',
#                      'sncg', 'chodl', 'serpinf1', 'krt73',
#                      'gpc3', 'ccnb1']

# gene_list = [                                   # 220116
#     'PHP.N', 'PHP.eB', 'slc30a3', 'PHP.B8',
#     'PHP.Astro', 'CAP-B10', 'gad2', 'PHP.V1',
#     'cux2', 'slc17a7', 'ctgf', 'sulf2',
#     'foxp1', 'foxp2', 'calb2', 'calb1',
#     'sst', 'pvalb', 'vip', 'gad1',
#     'tac1', 'lamp5', 'crym', 'rorb',
#     'rprm', 'reln', 'hsd11b1', 'tpbg',
#     'necab1', 'pcp4', 'crh', 'th',
#     'trh', 'lgi2', 'ppp1r17', 'krt73',
#     'sncg', 'chrna6', None, None
# ]

# gene_list_ordered = [
#     'PHP.eB', 'CAP-B10', 'PHP.N', 'PHP.Astro', 'PHP.V1', 'PHP.B8',
#     'slc30a3', 'slc17a7',
#     'cux2', 'calb1', 'lamp5', 'foxp1',
#     'rorb', 'sulf2', 'th',
#     'pcp4', 'rprm', 'crym', 'foxp2', 'necab1', 'hsd11b1',
#     'tpbg',
#     'ctgf', 'trh', 'chrna6','ppp1r17',
#     'gad1', 'gad2',
#     'pvalb', 'reln',
#     'sst', 'tac1', 'calb2', 'crh',
#     'vip',
#     'sncg','krt73',
#     'lgi2',
# ]
# gene_list_ordered = [gene for gene in gene_list if gene is not None]

# marker_genes_dict = {
#     'Neuronal': ['slc30a3'],
#     'Layer 2': ['cux2'],
#     'Layer 3': ['cux2'],
#     'Layer 4': ['rorb', 'hsd11b1'],
#     'Layer 5': ['rorb', 'chrna6', 'hsd11b1', 'rprm', 'crym'],
#     'Layer 6': ['trh', 'rprm'],
#     'Inhibitory': ['gad1'],
# }
