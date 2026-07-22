"""
Generate a grid-scan heatmap from iteration-10 ColabFold outputs.

Each heatmap cell represents the z-score of a trial's mean FRET efficiency
relative to the experimental FRET distribution.
"""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import PDBParser

import data_engine as engine


EXPERIMENTAL_MEAN = 0.291
EXPERIMENTAL_STDEV = 0.083
PROBE_NUMBER_1 = 100
PROBE_NUMBER_2 = 473
ITERATION_NUMBER = 10


def get_cb_coord(residues, index: int) -> np.ndarray:
    """Return the CB coordinate for a 1-based residue index, using CA as fallback."""
    try:
        residue = residues[index - 1]
    except IndexError as error:
        raise IndexError(f"Residue index {index} out of range.") from error

    if "CB" in residue:
        return residue["CB"].coord
    if "CA" in residue:
        return residue["CA"].coord
    raise ValueError(f"Residue {index} has neither CB nor CA atom.")


def compute_distance(residues, probe_number_1: int, probe_number_2: int) -> float:
    """Compute the Euclidean distance between the two probe residue coordinates."""
    coord1 = get_cb_coord(residues, probe_number_1)
    coord2 = get_cb_coord(residues, probe_number_2)
    return float(np.linalg.norm(coord2 - coord1))


def calculate_distance(structure_file: Path) -> float:
    """Use the same PDB parsing and probe-distance calculation as DistanceFinder.py."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", structure_file)
    residues = [
        residue
        for model in structure
        for chain in model
        for residue in chain
        if residue.id[0] == " "
    ]
    return compute_distance(residues, PROBE_NUMBER_1, PROBE_NUMBER_2)


def read_trial_parameters(jobs_path: Path) -> tuple[int, int]:
    """Return the n and m_e_msa values recorded for a container trial."""
    with jobs_path.open("r", encoding="utf-8") as file:
        jobs = json.load(file)["jobs"]

    if not jobs:
        raise ValueError("jobs.json does not contain any jobs.")

    trial = jobs[-1]
    return int(trial["n"]), int(trial["m_e_msa"])


def calculate_efficiencies(container: Path) -> np.ndarray:
    """Calculate iteration-10 FRET efficiencies for every PDB in a container."""
    iteration_path = container / "output_pool" / f"iteration{ITERATION_NUMBER}"
    pdb_files = sorted(iteration_path.glob("*.pdb"))
    if not pdb_files:
        raise FileNotFoundError(f"No PDB files found in {iteration_path}")

    distances = [
        calculate_distance(pdb_file)
        for pdb_file in pdb_files
    ]
    return np.asarray(engine.compute_E(np.asarray(distances, dtype=float)))


def collect_trial_zscores(containers_path: Path) -> dict[tuple[int, int], float]:
    """Return {(n, m_e_msa): z-score} for complete container trials."""
    trial_zscores = {}
    for container in sorted(path for path in containers_path.iterdir() if path.is_dir()):
        jobs_path = container / "jobs.json"
        if not jobs_path.exists():
            print(f">>> SKIPPING {container.name}: jobs.json not found")
            continue

        try:
            n, m_e_msa = read_trial_parameters(jobs_path)
            efficiencies = calculate_efficiencies(container)
        except (FileNotFoundError, KeyError, ValueError) as error:
            print(f">>> SKIPPING {container.name}: {error}")
            continue

        z_score = (float(np.mean(efficiencies)) - EXPERIMENTAL_MEAN) / EXPERIMENTAL_STDEV
        trial_zscores[(n, m_e_msa)] = z_score
        print(
            f">>> {container.name}: n={n}, m_e_msa={m_e_msa}, "
            f"mean FRET={np.mean(efficiencies):.4f}, z-score={z_score:.3f}"
        )

    if not trial_zscores:
        raise RuntimeError("No complete iteration-10 trials were found.")
    return trial_zscores


def create_heatmap(trial_zscores: dict[tuple[int, int], float], output_path: Path) -> None:
    """Save and display the n-by-m_e_msa z-score heatmap."""
    n_values = sorted({n for n, _ in trial_zscores})
    msa_values = sorted({msa for _, msa in trial_zscores})
    heatmap = np.full((len(msa_values), len(n_values)), np.nan)

    for (n, m_e_msa), z_score in trial_zscores.items():
        heatmap[msa_values.index(m_e_msa), n_values.index(n)] = z_score

    figure, axis = plt.subplots(figsize=(max(6, len(n_values) * 1.5), max(5, len(msa_values) * 1.2)))
    image = axis.imshow(heatmap, cmap="coolwarm", aspect="auto")
    colorbar = figure.colorbar(image, ax=axis)
    colorbar.set_label("Z-score of mean FRET efficiency")

    axis.set_xticks(range(len(n_values)), labels=n_values)
    axis.set_yticks(range(len(msa_values)), labels=msa_values)
    axis.set_xlabel("n")
    axis.set_ylabel("m_e_msa")
    axis.set_title(f"Iteration {ITERATION_NUMBER} FRET efficiency z-scores")

    for row, m_e_msa in enumerate(msa_values):
        for column, n in enumerate(n_values):
            z_score = heatmap[row, column]
            if not np.isnan(z_score):
                axis.text(column, row, f"{z_score:.2f}", ha="center", va="center")

    figure.tight_layout()
    figure.savefig(output_path, dpi=300)
    print(f">>> SAVED HEATMAP: {output_path}")
    plt.show()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--containers", type=Path, default=Path("containers"))
    parser.add_argument("--output", type=Path, default=Path("grid_scan_fret_zscore_heatmap.png"))
    arguments = parser.parse_args()

    trial_zscores = collect_trial_zscores(arguments.containers)
    create_heatmap(trial_zscores, arguments.output)


if __name__ == "__main__":
    main()
