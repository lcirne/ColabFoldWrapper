from Bio.PDB import PDBParser
import numpy as np
import sys

if len(sys.argv) > 1:
    # --- User adjustments ----
    structure_file = sys.argv[1]
    probe_number_1 = int(sys.argv[2])
    probe_number_2 = int(sys.argv[3])
    # -------------------------
else:
    # --- Default Values ------
    structure_file = 'ExamplePDB.pdb'
    probe_number_1 = 100  # 1-based index (as in MATLAB)
    probe_number_2 = 473
    # -------------------------

# Parse the PDB file
parser = PDBParser(QUIET=True)
structure = parser.get_structure('structure', structure_file)

# Get all standard residues across all chains (can be adjusted for specific chains if needed)
residues = [res for model in structure for chain in model for res in chain if res.id[0] == ' ']

# Function to get CB atom coordinates of a given residue index
def get_cb_coord(residues, index):
    """
    Returns coordinates of CB atom for residue at the given 1-based index.
    Falls back to CA if CB is missing (e.g., for Glycine).
    """
    try:
        residue = residues[index - 1]  # Convert 1-based MATLAB index to 0-based Python
        if 'CB' in residue:
            return residue['CB'].coord
        elif 'CA' in residue:
            print(f"Warning: Residue {index} has no CB atom. Using CA instead.")
            return residue['CA'].coord
        else:
            raise ValueError(f"Residue {index} has neither CB nor CA atom.")
    except IndexError:
        raise IndexError(f"Residue index {index} out of range.")


# Function to compute and return distance between given probe numbers
def compute_distance(residues, probe_number_1, probe_number_2):
    # Get CB coordinates for each probe residue
    coord1 = get_cb_coord(residues, probe_number_1)
    coord2 = get_cb_coord(residues, probe_number_2)

    # Compute and return Euclidean distance
    return np.linalg.norm(coord2 - coord1)

line_distance = compute_distance(residues, probe_number_1, probe_number_2)

# Output result
print(line_distance)

