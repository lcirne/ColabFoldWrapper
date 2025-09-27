"""
Luke Cirne
Data Engine for ColabFold Wrapper
Functions for graphing and data analysis
Ma Lab
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math
from typing import Dict, List, Union, Tuple

def graph_output_accuracy(distances: dict) -> str:
    # Collect and convert distances
    dist_values = [float(d) for d in distances.keys()]

    # Define bin range: from 5 below min to 5 above max, in steps of 10
    min_d = min(dist_values)
    max_d = max(dist_values)
    bin_start = np.floor(min_d - 5)
    bin_end = np.ceil(max_d + 5)
    bin_edges = np.arange(bin_start, bin_end + 10, 10)  # +10 to include final edge

    # Plot
    plt.figure(figsize=(8, 5))
    plt.hist(dist_values, bins=bin_edges, edgecolor="black", color="skyblue")
    plt.title("CF Output Distances (Å)")
    plt.xlabel("Distance (Å)")
    plt.ylabel("Frequency")
    plt.xticks(bin_edges)
    plt.tight_layout()

    # Save
    plot_name = "iteration_distances_hist"
    plt.savefig(f"{plot_name}.png")
    return plot_name


def _normal_cdf(x: float, mu: float, sigma: float) -> float:
    if sigma <= 0:
        raise ValueError("sigma must be > 0")
    z = (x - mu) / (sigma * math.sqrt(2))
    return 0.5 * (1.0 + math.erf(z))

def _bucket_index(x: float, mean: float, bucket_size: float) -> int:
    # Central bucket is [mean - b/2, mean + b/2)
    start = mean - bucket_size / 2.0
    return math.floor((x - start) / bucket_size)

def _bucket_bounds(idx: int, mean: float, bucket_size: float) -> Tuple[float, float]:
    start = mean - bucket_size / 2.0
    lo = start + idx * bucket_size
    hi = lo + bucket_size
    return lo, hi

def _as_multimap(pool: Dict[float, Union[str, List[str]]]) -> Dict[float, List[str]]:
    multi = {}
    for dist, val in pool.items():
        if isinstance(val, list):
            files = val
        else:
            files = [val]
        multi.setdefault(dist, []).extend(files)
    return multi


def _largest_remainder_apportion(N: int,
                                 probs: Dict[int, float],
                                 caps: Dict[int, int]) -> Dict[int, int]:
    # Initial floor allocation
    quotas = {i: min(caps[i], int(math.floor(N * probs[i]))) for i in probs}
    allocated = sum(quotas.values())
    remainder = N - allocated
    if remainder <= 0:
        return quotas

    # Largest remainder while honoring caps
    # Tie-breaker: smaller |bucket center - mean| preferred (we’ll inject later)
    remainders = {i: (N * probs[i] - quotas[i]) for i in probs}
    while remainder > 0:
        # Candidates that still have room
        cands = [i for i in probs if quotas[i] < caps[i]]
        if not cands:
            break
        # Pick the one with largest fractional remainder
        cands.sort(key=lambda i: remainders[i], reverse=True)
        for i in cands:
            if remainder == 0:
                break
            if quotas[i] < caps[i]:
                quotas[i] += 1
                remainder -= 1
    return quotas


def choose_gaussian_distances(
    pool: Dict[float, Union[str, List[str]]],
    mean: float,
    bucket_size: float,
    sigma: float = None,
    target_total: int = None,
) -> Dict[float, List[str]]:
    """
    Select items whose distances approximate a Gaussian distribution.

    Parameters
    ----------
    pool : dict[distance -> filename or list[filename]]
        Distances can repeat by providing a list of filenames for that distance.
    mean : float
        Center of the desired Gaussian.
    bucket_size : float
        Width of each bucket (central bucket is [mean - b/2, mean + b/2)).
    sigma : float, optional
        Standard deviation of the Gaussian. Defaults to bucket_size.
    target_total : int, optional
        Desired total number of items to select. If None, the algorithm
        chooses the maximum N that fits Gaussian proportions without
        exceeding any bucket's capacity.

    Returns
    -------
    dict[distance -> list[filenames]]
        Only the chosen items.
    """
    if bucket_size <= 0:
        raise ValueError("bucket_size must be > 0")
    if not pool:
        return {}

    if sigma is None:
        sigma = bucket_size

    multi = _as_multimap(pool)

    # Build per-item list, then bucket them
    items: List[Tuple[float, str]] = []
    for d, files in multi.items():
        for f in files:
            items.append((float(d), f))

    # Determine bucket indices for all items
    buckets: Dict[int, List[Tuple[float, str]]] = {}
    for d, f in items:
        idx = _bucket_index(d, mean, bucket_size)
        buckets.setdefault(idx, []).append((d, f))

    # Capacities per bucket
    caps = {idx: len(v) for idx, v in buckets.items()}

    # Gaussian mass per bucket (only where we have capacity)
    weights = {}
    for idx, cap in caps.items():
        lo, hi = _bucket_bounds(idx, mean, bucket_size)
        w = _normal_cdf(hi, mean, sigma) - _normal_cdf(lo, mean, sigma)
        weights[idx] = max(0.0, w)

    # Normalize over non-empty buckets with positive weight
    positive = {i: w for i, w in weights.items() if w > 0 and caps[i] > 0}
    if not positive:
        # Fallback: equal weights across non-empty buckets
        positive = {i: 1.0 for i in caps.keys() if caps[i] > 0}

    total_w = sum(positive.values())
    probs = {i: positive[i] / total_w for i in positive}

    total_available = sum(caps.values())
    if target_total is None:
        # Max N that doesn't exceed any bucket capacity in expectation
        N_float_limits = [caps[i] / probs[i] for i in probs]
        N = int(math.floor(min(N_float_limits)))
        N = max(0, min(N, total_available))
    else:
        N = max(0, min(int(target_total), total_available))

    if N == 0:
        return {}

    # Apportion counts to buckets
    quotas = _largest_remainder_apportion(N, probs, caps)

    # Within-bucket selection: pick closest-to-mean first
    selected_pairs: List[Tuple[float, str]] = []
    for idx, need in quotas.items():
        if need <= 0:
            continue
        lo, hi = _bucket_bounds(idx, mean, bucket_size)
        # Sort by closeness to mean, then by closeness to bucket center toward mean
        bucket_items = buckets[idx]
        # Secondary tiebreaker: distance from mean, then from central line
        bucket_items.sort(key=lambda t: (abs(t[0] - mean), abs((lo + hi) / 2.0 - mean), t[0], t[1]))
        selected_pairs.extend(bucket_items[:need])

    # Aggregate back to {distance: [filenames]}
    out: Dict[float, List[str]] = {}
    for d, f in selected_pairs:
        out.setdefault(d, []).append(f)
    return out


def explode_distances(selection: Dict[float, List[str]]) -> List[float]:
    """Helper: flatten selection to a list of distances (one per selected filename)."""
    result = []
    for d, files in selection.items():
        result.extend([d] * len(files))
    # Optional: sort for readability
    return sorted(result)
