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


def graph_output_accuracy(efficiencies: dict, bins=0.05, graph_name=None) -> str:
    # Collect and convert distances
    effs = np.array([float(d) for d in efficiencies.values()])

    # If bins is a float, treat it as bin width and generate edges
    if isinstance(bins, float) or isinstance(bins, int):
        min_d = effs.min()
        max_d = effs.max()
        bin_edges = np.arange(min_d, max_d + bins, bins)
    else:
        # If bins is an array (from build_distribution), use it directly
        bin_edges = bins

    # Debug
    #print(effs.min(), effs.max())
    #print(bin_edges[:5], bin_edges[-5:])
    # Plot
    plt.figure(figsize=(8, 5))
    plt.hist(effs, bins=bin_edges, edgecolor="black", color="skyblue")
    plt.title("CF Output Distances (Å)")
    plt.xlabel("Distance (Å)")
    plt.ylabel("Frequency")
    xticks = np.arange(bin_edges.min(), bin_edges.max()+0.1, 0.1)
    plt.xticks(xticks)
    plt.tight_layout()

    # Save
    plot_name = "iteration_distances_hist"
    if graph_name:
        plot_name = graph_name
    plt.savefig(f"{plot_name}.png")
    return plot_name


def graph_output_accuracy_bar(efficiencies: dict, bins=0.05, graph_name=None) -> str:
    """
    Plots a bar chart where each bar corresponds to a histogram bin.
    X values are bin centers, and Y values are counts of efficiencies in each bin.
    """

    # Convert dictionary values to numpy array
    effs = np.array([float(d) for d in efficiencies.values()])

    # Determine bin edges and centers
    if isinstance(bins, float) or isinstance(bins, int):
        min_d = effs.min()
        max_d = effs.max()
        bin_edges = np.arange(min_d, max_d + bins, bins)
        bin_centers = bin_edges[:-1] + bins / 2
    else:
        # If bins provided as array
        bin_edges = bins
        bin_centers = bin_edges[:-1] + (bin_edges[1] - bin_edges[0]) / 2

    # Count how many values fall into each bin
    counts, _ = np.histogram(effs, bins=bin_edges)

    # --- Plot ---
    plt.figure(figsize=(8, 5))
    plt.bar(bin_centers, counts, width=(bin_edges[1] - bin_edges[0]) * 0.9,
            color="mediumseagreen", edgecolor="black")
    plt.title("CF Output Distances (Å) — Bar Plot")
    plt.xlabel("Distance (Å)")
    plt.ylabel("Frequency")

    # Set x-axis ticks
    xticks = np.arange(bin_edges.min(), bin_edges.max() + 0.1, 0.1)
    plt.xticks(xticks)
    plt.tight_layout()

    # --- Save ---
    plot_name = "iteration_distances_bar"
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
    # TODO: Set N to an arbitrary number to reduce the influence of templates on the output
    N = len(efficiencies)

    # Define bin edges across observed range
    min_val, max_val = efficiencies.min(), efficiencies.max()
    print(f"min_vale: {min_val} max_val: {max_val}")
    bins = np.arange(min_val, max_val + bin_width, bin_width)

    # Bin assignments for each efficiency
    bin_indices = np.digitize(efficiencies, bins) - 1
    print(f"bin_indices: {bin_indices}")

    # Compute bin centers
    bin_centers = bins[:-1] + bin_width / 2

    # --- Step 3: Compute Gaussian-based target counts ---
    gauss_probs = np.exp(-0.5 * ((bin_centers - mean) / std) ** 2)
    gauss_probs /= gauss_probs.sum()  # normalize
    target_counts = np.round(gauss_probs * N).astype(int)

    """
    TODO: Test without code chunk to see if highest bin will be reduced
    # Adjust for rounding (so sum == N)
    diff = N - target_counts.sum()
    if diff != 0:
        # fix by adjusting the bin with the highest probability
        target_counts[np.argmax(gauss_probs)] += diff
    """

    # Group files by bin
    bin_to_files = defaultdict(list)
    for fname, eff, bidx in zip(filenames, efficiencies, bin_indices):
        bin_to_files[bidx].append((fname, eff))

    # Collect selected filename-efficiency pairs
    selected = {}
    dupe_counter = defaultdict(int)

    for bidx, desired_count in enumerate(target_counts):
        print(f"desired_count: {desired_count}")
        available_files = bin_to_files.get(bidx, [])

        if desired_count == 0 or len(available_files) == 0:
            continue

        """
        TODO: Create a count variable to monitor how many changes need to be made
        to the outputted distribution to build a gaussian distribution.
        """
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
        print(len(chosen))
        print(len(selected))

    return selected, bins, bin_centers
