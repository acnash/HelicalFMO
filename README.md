# HelicalFMO

## UNDER ACTIVE DEVELOPMENT (27.03.2025)

Fragment Molecular Orbital analysis for interhelical transmembrane protein interactions. 

This code performs several operations:

1. Rotates helical transmembrane dimers and writes to PDB files.
2. Writes nearest residue-residue atomic records to PDF files. 
3. Generates GAMESS input file for FMO analysis based on nearest residue-residue PDB files.
4. Performs 1 to 3 at scale.

## Isolating protein-protein interaction residues

### Command line parameters
- ```--file``` The PDB input file. The file must contain only two chains, A and B. All other chains are ignored.
- ````--mode```` Set to ```contact_distance ``` to isolate protein-protein interaction residues.
- ```--output_folder``` For each input file (whether that's a single file or a directory of files), a corresponding results PDB and data file will be saved to this location. The PDB file contains only those residues within contact distance, and the data file contains two lists of isolated residues. 

For example, to isolate interaction residues within 3 angstrom, while renumbering the residues of chain B (from residue 1): 
```
python helical_FMO.py --file C:\Users\Anthony\PyCharmProjects\HelicalFMO\structures\4auo_single_chains_AB.pdb --mode contact_distance --output_folder C:\Users\Anthony\PyCharmProjects\HelicalFMO\temp\ --distance_cutoff 3 --renum_chains B:1
```
## Generating GAMESS FMO input files

For example...
```
python helical_FMO.py --file C:\Users\Anthony\PyCharmProjects\HelicalFMO\temp\contacts_0.pdb --mode fmo --output_folder C:\Users\Anthony\PyCharmProjects\HelicalFMO\temp\ --basis STO-3G --theory HF
```
