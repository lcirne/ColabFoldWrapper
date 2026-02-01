# ColabFold Wrapper

**A wrapper used in combination with the ColabFold Batch command line utility for implementing experimental `Fluoresence Resonance Energy Transfer (FRET)` data into the base ColabFold protocol.**

The additional influence of experimental data in the AlphaFold Evoformer acts as a counterbalance over existing model noise when attempting to visualize proteins in a disordered state.

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
                ├── <your_protein_sequence>.pdb    # Input protein structure  
                └── <custom_templates_dir>/        # Optional custom templates  
                    ├── <template1>.pdb  
                    ├── <template2>.pdb  
                    ├── <template1>.cif  
                    └── <template2>.cif
- Once this file structure has been established, you can install dependencies if not using Anaconda `pip install -r < requirements.txt`
- Run the wrapper `python3 wrapper.py`
