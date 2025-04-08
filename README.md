# HelicalFMO

## UNDER ACTIVE DEVELOPMENT (27.03.2025)

Fragment Molecular Orbital analysis for interhelical transmembrane protein interactions. 

<p align="center">
    <img src="doc/MMP14_TM.jpg" alt="Description" width="40%">
</p>

This code performs several operations:

1. Rotates helical transmembrane dimers and writes to PDB files.
2. Writes nearest residue-residue atomic records to PDF files. 
3. Caps isolated peptide fragments with hydrogen atoms to ensure the overall formal charge is consistent with the original structure. 
4. Generates GAMESS input file for FMO analysis based on nearest residue-residue PDB files.
5. Performs 1 to 4 at scale.

## Installation
I've tried to keep the number of dependencies to a minimum. 

- Clone HelicalFMO main. 
- Install MDAnalysis inside a base environment or an environment of your choosing. Obviously HelicalFMO needs access. 
- Install PSI4 `conda create -n psi4_env psi4 -c conda-forge`
- Install <a href="https://pypi.org/project/biopython/">Biopython.</a> 
- Install <a href="https://github.com/clauswilke/PeptideBuilder">PeptideBuilder</a> with `pip install PeptideBuilder`
- Install <a href="https://github.com/leeping/geomeTRIC">geomeTRIC</a> by changing to the psi4 environment `conda activate psi4_env` then `conda install -c conda-forge geometric`

---

## Generating a helix-helix dimer from sequence

HelicalFMO calls PeptideBuilder to build homo and hetero all-atom helical dimers.  

---

## Generating helix-helix starting structures as a function of rotation 

HelicalFMO can take any helix-helix (two chain) PDB file and generate an ensemble of starting structures as a function of helical rotation. 

---


## Isolating protein-protein interaction residues

### Command line parameters
- ```--mode``` Set to ```contact_distance ``` to isolate protein-protein interaction residues.
- ```--file``` The PDB input file. The file must contain only two chains, A and B. All other chains are ignored.
- ```--folder``` Path to a folder containing one or more PDB input files. This is an alternative choice to the single input ```--file``` parameter.
- ```--output_folder``` For each input file (whether that's a single file or a directory of files), a corresponding results PDB and data file will be saved to this location. The PDB file contains only those residues within contact distance, and the data file contains two lists of isolated residues. 
- ```--distance_cutoff``` (default 3 angstrom) An atom-to-atom distance measure between residues. Residue pairs within (<=) 
the cutoff distance are kept as interaction residues.
- ```--renum_chains``` Renumbers the residues on chain X from residue N e.g., A:1 renumbers residues, begging from "1", on chain A. For both chains, e.g., A:1 B:5, renumbers chain A beginning with "1", and chain B, beginning with "5". This operation
ensure hydrogen atom caps are added to the correct FMO fragment. Problems will occur if a chain has disordered or missing residue numbers e.g., 3, 4, 6, 11, 12, or 1, 3, 4, 7, 5.
The code uses the sequential order of residue numbers to determine whether there is a break in the peptide backbone, and therefore, a hydrogen atom cap is required. 
- ```--ignore_num_start_res``` Ignores all residues up to this residue e.g., ```5``` ensures the first four residues are ignored. Uses residue number. Renumber with ```--renum_chains``` if necessary.  
- ```--ignore_num_end_res``` Ignores all residues after this residue e.g., ```26``` ensures all residues after the 26th are ignored. Uses residue number. Renumber with ```--renum_chains``` if necessary.  

For example, to isolate interaction residues within 3 angstrom, while renumbering the residues of chain B (from residue 1): 
```
python helical_FMO.py --file C:\Users\Anthony\PyCharmProjects\HelicalFMO\structures\4auo_single_chains_AB.pdb --mode contact_distance --output_folder C:\Users\Anthony\PyCharmProjects\HelicalFMO\temp\ --distance_cutoff 3 --renum_chains B:1
```

### Output

Two output files are generated per PDB input file. Firstly, the residues (resname and resid) within interhelical contact distance are written to a text file, for example:
```
Chain A Residues:
GLN 46
ARG 47
SER 49
LEU 51
THR 52
ILE 55
SER 56
VAL 59
GLY 60
LEU 63
VAL 64
Chain B Residues:
GLN 6
ARG 7
SER 9
LEU 11
THR 12
ILE 15
SER 16
VAL 19
GLY 20
LEU 23
VAL 24
``` 
The second, a PDB file containing only those residues. 

---
## Preparing PDB files for GAMESS FMO input file generator

Before generating a GAMESS FMO input file, we need to ensure the input PDB file is formatted for such purpose. This is tricky, and getting it wrong will add artificial energy contributions to the FMO analysis. Read carefully. 

Having isolated the interaction residues at the helix-helix interface (using the steps above), we need to cap those residues that lost their 
peptide-bond neighbouring residue e.g., of ARG-LYS-PRO-CYS, only LYS and CYS were within the contact cutoff, therefore they have both lost their peptide-bond to a neighbouring residue. 

Capping involves adding a hydrogen atom to the cut nitrogen atom on the N-terms and a hydrogen atom to the cut carbon-carbonyl atom
on the C-term, both, part of the peptide backbone. Then a very short geometry optimisation is performed using PSI4 while keeping the original atoms
restrained. Adding hydrogen caps ensures that the residue retains its typical formal charge. 

Note: the first residue on each chain must have started with a resid of 1 during protein-protein interaction residue steps detailed above. Also, the N and C terms of each chain are down to the User's discretion. 
Setting the residue IDs to 1 instructs the program to avoiding adding a hydrogen atom cap to original N and C terms. If, for example, the protein-protein
interaction step rejects the first residue on a chain (because it wasn't within a cutoff of a residue from the neighbouring chain), the first residue read will have a residue ID of 2. 
The software, will know to add a hydrogen cap to that residue. 

### Command line parameters
- ```--mode``` Set to ```cap ``` to add hydrogen cap atoms to cut residues.
- ```--file``` The PDB input file. The file must contain only two chains, A and B. It's likely this file will have been generated having first performed ```--mode contact_distance```.
- ```--folder``` Path to a folder containing one or more PDB input files. This is an alternative choice to the single input ```--file``` parameter.
- ```--output_folder``` For each input file (whether that's a single file or a directory of files), new capped PDB files will be saved to this location.

For example, we cap a file refined to only include residues on both chains within a 3 angstrom interchain cutoff distance. A new file is saved to the ```--output_folder``` location.  
```
--file C:\Users\Anthony\PyCharmProjects\HelicalFMO\temp\contacts_0.pdb --mode cap --output_folder C:\Users\Anthony\PyCharmProjects\HelicalFMO\temp\
```
---

## Generating GAMESS FMO input files

### Command line parameters
- ```--mode``` Set to ```fmo ``` to build GAMESS FMO input files.
- ```...``` -
- ```...``` - 

For example...
```
python helical_FMO.py --file C:\Users\Anthony\PyCharmProjects\HelicalFMO\temp\contacts_0.pdb --mode fmo --output_folder C:\Users\Anthony\PyCharmProjects\HelicalFMO\temp\ --basis STO-3G --theory HF
```
