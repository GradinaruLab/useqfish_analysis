nC = 5      # number of channels / colors
roundRef = -1    #  round containing the cell reference; use -1 for the last round
cellch = 3      # channel containing cell reference 

sigma = 3       # estimate size of spots to detect
shift_window_size = 200     # ? TODO

color_shifts = [[8, 0, 0],          # color aberration shifts for each channel to 647 aquired by imaging focalcheck #1, 500nm
                [-1, 0, 1],
                [0, 0, 0],
                [0, 0, 0],
                [1, -3, 2]]