import os
import numpy as np
import wrapper
import data_engine as engine

dirpath = "/Users/nikan/Documents/ma_lab/ColabFoldWrapper/graphing-utils/templates/iteration3"
dir_contents = os.listdir(dirpath)

distances = {}
for file in dir_contents:
    if file.endswith(".pdb"):
        distances[file] = float(wrapper.run_distance_finder(f"{dirpath}/{file}", "100", "473"))

distances_to_convert = np.array(list(distances.values()))
e_conversions = engine.compute_E(distances_to_convert)
i = 0
for filename, distance in distances.items():
    distances[filename] = e_conversions[i]
    i += 1

y_exp = 0.291
sigma = 0.083
n = 90
bin_width = 0.025
bins = np.arange(0.0, 1.0, bin_width)
bin_centers = bins[:-1] + bin_width / 2
wrapper.plot_fret_efficiencies(distances, 0, bin_centers, n)
