# HelicalFMO

## UNDER ACTIVE DEVELOPMENT (27.03.2025)

Fragment Molecular Orbital analysis for interhelical transmembrane protein interactions. 

This code performs several operations:

1. Rotates helical transmembrane dimers and writes to PDB files.
2. Writes nearest residue-residue atomic records to PDF files. 
3. Caps isolated peptide fragments with hydrogen atoms to ensure the overall formal charge is consistent with the original structure. 
3. Generates GAMESS input file for FMO analysis based on nearest residue-residue PDB files.
4. Performs 1 to 4 at scale.

## Isolating protein-protein interaction residues

### Command line parameters
- ```--file``` The PDB input file. The file must contain only two chains, A and B. All other chains are ignored.
- ```--mode``` Set to ```contact_distance ``` to isolate protein-protein interaction residues.
- ```--output_folder``` For each input file (whether that's a single file or a directory of files), a corresponding results PDB and data file will be saved to this location. The PDB file contains only those residues within contact distance, and the data file contains two lists of isolated residues. 
- ```--distance_cutoff``` (default 3 angstrom) An atom-to-atom distance measure between residues. Residue pairs within (<=) 
the cutoff distance are kept as interaction residues.
- ```--renum_chains``` Renumbers the residues on chain X from residue N e.g., A:1 renumbers residues, begging from "1", on chain A. For both chains, e.g., A:1 B:5, renumbers chain A beginning with "1", and chain B, beginning with "5". This operations
ensures hydrogen atom caps are added to the correct FMO fragment. Problems will occur if a chain has disordered or missing residue numbers e.g., 3, 4, 6, 11, 12, or 1, 3, 4, 7, 5.
The code uses the sequential order of residue numbers to determine whether there is a break in the peptide backbone, and therefore, a hydrogen atom cap is required. 
- ```--ignore_num_start_res``` Ignores all residues up to this residue e.g., ```5``` ensures the first four residues are ignored. Uses residue number. Renumber with ```--renum_chains``` if necessary.  
- ```--ignore_num_end_res``` Ignores all residues after this residue e.g., ```26``` ensures all residues after the 26th are ignored. Uses residue number. Renumber with ```--renum_chains``` if necessary.  

For example, to isolate interaction residues within 3 angstrom, while renumbering the residues of chain B (from residue 1): 
```
python helical_FMO.py --file C:\Users\Anthony\PyCharmProjects\HelicalFMO\structures\4auo_single_chains_AB.pdb --mode contact_distance --output_folder C:\Users\Anthony\PyCharmProjects\HelicalFMO\temp\ --distance_cutoff 3 --renum_chains B:1
```

## Preparing PDB files for GAMESS FMO input file generator

Before generating a GAMESS FMO input file, we need to ensure the input PDB file is formatted for such purpose. 
Having isolated the interaction residues at the helix-helix interface, we need to cap those residues that lost their 
peptide-bond neighbouring residue. We do this by adding a hydrogen atom to the cut nitrogen atom and cut carbon-carbonyl atom
on the backbone, then a very short geometry optimisation is performed using PSI4 while keeping the original atoms
restrained. Adding hydrogen caps ensures that the residue retains its typical formal charge. 

## Generating GAMESS FMO input files

### Command line parameters
- ```...``` -
- ```...``` -
- ```...``` - 

For example...
```
python helical_FMO.py --file C:\Users\Anthony\PyCharmProjects\HelicalFMO\temp\contacts_0.pdb --mode fmo --output_folder C:\Users\Anthony\PyCharmProjects\HelicalFMO\temp\ --basis STO-3G --theory HF
```
