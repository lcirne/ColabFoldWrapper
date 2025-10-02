"""
Luke Cirne
Data Engine for ColabFold Wrapper
Functions for graphing and data analysis
Ma Lab
"""

import matplotlib.pyplot as plt
from matplotlib import ticker
import pandas as pd
import numpy as np
import random
import os
from collections import defaultdict


def compute_E(distances, R_0=51):
    for i in range(len(distances)):
        distances[i] = float(1 / (1 + (distances[i] / R_0)**6))
    return distances


def graph_output_accuracy(distances: dict, bins=0.05, graph_name=None) -> str:
    # Collect and convert distances
    dist_values = [float(d) for d in distances.values()]

    # Define bin range: from 5 below min to 5 above max, in steps of 10
    min_d = min(dist_values)
    max_d = max(dist_values)
    bin_start = 0
    bin_end = np.ceil(max_d)
    bin_edges = np.arange(bin_start, bin_end, bins)

    # Plot
    plt.figure(figsize=(8, 5))
    plt.hist(dist_values, bins=bin_edges, edgecolor="black", color="skyblue")
    plt.title("CF Output Distances (Å)")
    plt.xlabel("Distance (Å)")
    plt.ylabel("Frequency")
    xticks = np.arange(0, bin_end, 0.1)
    plt.xticks(xticks)
    plt.tight_layout()

    # Save
    plot_name = "iteration_distances_hist"
    if graph_name:
        plot_name = graph_name
    plt.savefig(f"{plot_name}.png")
    return plot_name


def build_distribution(
    file_eff_dict: dict,
    mean: float,
    std: float,
    bin_width: float = 0.05,
    seed: int = None
) -> dict:

    """
    Selects file-efficiency pairs such that their histogram best follows a Gaussian distribution defined by the provided mean and standard deviation.

    Parameters:
        file_eff_dict : dict
            Dictionary of {filename: efficiency}, where efficiency is a float.
        mean : float
            Mean value for the target Gaussian distribution.
        std : float
            Standard deviation for the target Gaussian distribution.
        bin_width : float, optional
            Width of histogram bins. Default is 5.
        seed : int, optional
            Random seed for reproducibility. Default is 42.

    Returns
        dict: Dictionary of {filename: efficiency} containing the selected
            file-efficiency pairs adjusted to match the Gaussian distribution.
    """

    # Set seeds for reproducibility
    if seed is not None:
        np.random.seed(seed)
        random.seed(seed)

    # Extract efficiencies
    efficiencies = np.array(list(file_eff_dict.values()))
    filenames = np.array(list(file_eff_dict.keys()))
    N = len(efficiencies)

    # Define bin edges across observed range
    min_val, max_val = efficiencies.min(), efficiencies.max()
    bins = np.arange(min_val, max_val + bin_width, bin_width)

    # Bin assignments for each efficiency
    bin_indices = np.digitize(efficiencies, bins) - 1

    # Compute bin centers
    bin_centers = bins[:-1] + bin_width / 2

    # --- Step 3: Compute Gaussian-based target counts ---
    gauss_probs = np.exp(-0.5 * ((bin_centers - mean) / std) ** 2)
    gauss_probs /= gauss_probs.sum()  # normalize
    target_counts = np.round(gauss_probs * N).astype(int)

    # Adjust for rounding (so sum == N)
    diff = N - target_counts.sum()
    if diff != 0:
        # fix by adjusting the bin with the highest probability
        target_counts[np.argmax(gauss_probs)] += diff

    # Group files by bin
    bin_to_files = defaultdict(list)
    for fname, eff, bidx in zip(filenames, efficiencies, bin_indices):
        bin_to_files[bidx].append((fname, eff))

    # Collect selected filename-efficiency pairs
    selected = {}
    dupe_counter = defaultdict(int)

    for bidx, desired_count in enumerate(target_counts):
        available_files = bin_to_files.get(bidx, [])

        if desired_count == 0:
            continue

        if len(available_files) >= desired_count:
            # Too many files, sample down
            chosen = random.sample(available_files, desired_count)
        else:
            # Too few files, duplicate as needed
            multiplier = -(-desired_count // len(available_files))  # ceiling division
            extended = available_files * multiplier
            chosen = random.sample(extended, desired_count)

        # Add chosen pairs with dupe suffixes if needed
        for fname, eff in chosen:
            if fname in selected:
                dupe_counter[fname] += 1
                name, ext = os.path.splitext(fname)
                new_fname = f"{name}_dupe{dupe_counter[fname]}{ext}"
                selected[new_fname] = eff
            else:
                selected[fname] = eff

    return selected
