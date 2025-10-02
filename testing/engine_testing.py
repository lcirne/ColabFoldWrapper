"""
Test suite for data_engine.py
"""
import os
import sys
import subprocess
from Bio.PDB import PDBParser
import numpy as np
import sys

# Path of this file (engine_testing.py)
here = os.path.dirname(os.path.abspath(__file__))
# Add the parent directory of ColabFoldWrapper to sys.path
sys.path.append(os.path.abspath(os.path.join(here, "..", "..")))
from ColabFoldWrapper.data_engine import graph_output_accuracy, build_distribution, compute_E


def run_distance_finder(structure_file, p1, p2): #helper function for build_distribution_test
    """
    Runs an external script to calculate the distance between two residues
    in a protein structure.

    Args:
        structure_file (str): Path to the .pdb file.
        p1 (str): Residue index 1.
        p2 (str): Residue index 2.

    Returns:
        str or None: Distance in angstroms as a string, or None if failed.
    """
    distance_finder_filepath = os.path.join(here, "..", "distance_finder", "DistanceFinder.py")
    distance = subprocess.run(
        ["python3", distance_finder_filepath, structure_file, p1, p2],
        capture_output=True,
        text=True
    )
    if distance.returncode != 0:
        print("Error:", distance.stderr)
        return None
    return distance.stdout.strip()


def build_distribution_test(file_eff_dict=None):
    print("build_distribution_test")
    if not file_eff_dict:
        file_eff_dict = {}
        wrapper_testing_dir = os.path.join(here, "..", "..", "ColabFold_output", "wrapper_testing")

        test_files = os.listdir(wrapper_testing_dir)
        for file in test_files:
            if file.endswith(".pdb"):
                distance = run_distance_finder(f"{wrapper_testing_dir}/{file}", "100", "473")
                if distance:
                    file_eff_dict[file] = float(distance)

        efficiencies = list(file_eff_dict.values())
        efficiencies = compute_E(efficiencies)
        i = 0
        for file, distance in file_eff_dict.items():
            file_eff_dict[file] = efficiencies[i]
            i += 1

    # Histogram of initial data pre-processing
    graph_output_accuracy(file_eff_dict, graph_name="pre-processing-graph")

    y_exp = 0.291
    sigma = 0.083
    for i in range(5):
        gaussian_dist_efficiencies = build_distribution(file_eff_dict, y_exp, sigma, seed=i+42)
        dist_clean = {str(k): float(v) for k, v in gaussian_dist_efficiencies.items()}
        graph_output_accuracy(dist_clean, graph_name=f"post-processing-{i}")

    print("Test complete")




if __name__ == '__main__':
    build_distribution_test()
