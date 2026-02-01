# ColabFold Wrapper

**A wrapper used in combination with the ColabFold Batch command line utility for implementing experimental `Fluoresence Resonance Energy Transfer (FRET)` data into the base ColabFold protocol.**

The additional influence of experimental data in the AlphaFold2 Evoformer acts as a counterbalance over existing model noise when attempting to visualize proteins in a disordered state. Program is currently in development and only designed to work with the UvrD helicase as of now. This wrapper can be treated as a proof of concept of the potential benefits from implementing experimental data into ColabFold.

## Requirements
- The wrapper can be ran using `python3` on a `linux` environment
- ColabFold Batch (see https://github.com/sokrypton/ColabFold for installation instructions)
- Runtime for ColabFold is heavily accelerated with GPU assistance
- **A pdb file for your amino acid sequence**
- Template data for the protein you are attempting to model (optional)
- An installation of Anaconda for dependencies (optional but recommended)

## Running ColabFold Wrapper
- Have the batch cli tool installed
- Organize file structure in the following format:

                project_root/  
                ├── wrapper.py                     # Main wrapper script  
                ├── data_engine.py                 # Data processing logic
                └── distance_finder/  
                    └── DistanceFinder.py          # Calculate the distance between FRET 'probes'
  
                ├── <your_protein_sequence>.pdb    # Input protein structure  
                └── <custom_templates_dir>/        # Optional custom templates  
                    ├── <template1>.pdb  
                    ├── <template2>.pdb  
                    ├── <template1>.cif  
                    └── <template2>.cif
- Once this file structure has been established, you can install dependencies if not using Anaconda `pip install -r < requirements.txt`
- Run the wrapper `python3 wrapper.py`

## Important Caveats
- Right now DistanceFinder.py is only measuring the distance between two specific amino acids in the UvrD Helicase chain. The probe values will need to be changed to appropriately match the experimental FRET data available for other proteins. Future plans involve a more streamlined input process that will include easy modification of these probe values.
