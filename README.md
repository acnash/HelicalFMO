# HelicalFMO


## UNDER ACTIVE DEVELOPMENT (22.04.2025)

__To do__
1. Check whether PSI4 coordinates are retrieved (I don't think they are).
2. Go through and implement fault tolerance.
3. Go through and implement Logger.
4. Implement a cross-angle and Left/Right handedness. From _n_ points along the backbone, adjust the crossing angle by optimising the energy (use an average). 
5. Implement ANI-X2 force field to measure potential energy (of the full TM peptide-peptide).
6. Add minimisation parameters.
7. Clean the code, add function/method documentation. 
---

Fragment Molecular Orbital analysis for inter-helical transmembrane dimer residue interactions. 

<p align="center">
    <img src="doc/MMP14_TM.jpg" alt="Description" width="40%">
</p>

This code performs several operations:

1. Generated _de novo_ transmembrane dimers.
2. Rotates helical transmembrane dimers and writes to PDB files.
3. Writes nearest residue-residue atomic records to PDF files. 
4. Caps isolated peptide fragments with hydrogen atoms to ensure the overall formal charge is consistent with the original structure. 
5. Generates GAMESS input file for FMO analysis based on nearest residue-residue PDB files.
6. Performs 1 to 5 at scale.

## Installation
I've tried to keep the number of dependencies to a minimum. 

- Clone HelicalFMO main. 
- Install MDAnalysis inside a base environment or an environment of your choosing. Obviously HelicalFMO needs access. 
- Install PSI4 `conda create -n psi4_env psi4 -c conda-forge`
- Install yaml ```pip install yaml``` and PyYaml ```pip install PyYaml```.
- Install <a href="https://pypi.org/project/biopython/">Biopython</a> with ```pip install biopython```
- Install <a href="https://github.com/clauswilke/PeptideBuilder">PeptideBuilder</a> with `pip install PeptideBuilder`
- Install <a href="https://github.com/leeping/geomeTRIC">geomeTRIC</a> by changing to the psi4 environment `conda activate psi4_env` then `conda install -c conda-forge geometric`
- Install RDKit for Python

---

## Run command line parameters

```--input_file``` All input parameters are stored in a YAML file. Provide absolute path and file name. Example YAML files is in ```./examples```.

Each software operation is defined by its YAML section, and all options are held there in:
```
operation:
    option: value
    option: value
    # comment
```

You can execute multiple operations, which execute in the order in which they are presented in the YAML file:
```
operation_1:
    option: value
    option: value
operation_2:
    option: value
    option: value
```

The operations are:

- ```generate```
- ```rotation```
- ```contact_distance```
- ```cap```
- ```fmo```

And their input options are defined below.

---

## Generating a helix-helix dimer from sequence

HelicalFMO calls PeptideBuilder to build homo and hetero all-atom helical dimers.  

### Command line parameters
- ```seq_a``` Single amino acid codes e.g., ATLLIF for the first helical peptide.
- ```seq_b``` Single amino acid codes e.g., CFIAVG for the second helical peptide. If provided, a heterodimer is built. If not provided, a homodimer is built based on ```seq_a```.
- ```output_file``` File path and file name to save the generated hetero/homodimer as a PDB. 

For example, to create a homodimer:

```
generate:
  seq_a: KKKLALITFFIALVGLGAVGLAWTRKKK
  output_file: C:\Users\Anthony\PycharmProjects\HelicalFMO\temp\homodimer.pdb
```

To create a heterodimer:

```
generate:
  seq_a: KKKLALITFFIALVGLGAVGLAWTRKKK
  seq_b: KKTLALISVRGGLIAIYAVGLAWTRKYK
  output_file: C:\Users\Anthony\PycharmProjects\HelicalFMO\temp\homodimer.pdb
```

---

## Generating helix-helix starting structures as a function of rotation 

HelicalFMO can take any helix-helix (two chain) PDB file and generate an ensemble of starting structures as a function of helical rotation. 

### Command line parameters
- ```file``` The PDB input file. The file must contain only two chains, A and B. All other chains are ignored.
- ```output_folder``` The folder path output PDB files are saved to.
- ```rotation_angle``` The rotation angle. Default is 20 degrees. 
- ```both_helices``` If set to False, helix A is held fixed while helix A will rotate. When true (default) both helices rotate.

For example:

```
rotation:
  file: C:\Users\Anthony\PycharmProjects\HelicalFMO\temp\homodimer.pdb
  output_folder: C:\Users\Anthony\PycharmProjects\HelicalFMO\temp\
  rotation_angle: 20
  both_helices: True
```

---

## Isolating protein-protein interaction residues

### Command line parameters
- ```file``` The PDB input file. The file must contain only two chains, A and B. All other chains are ignored.
- ```folder``` Path to a folder containing one or more PDB input files. This is an alternative choice to the single input ```file``` parameter.
- ```output_folder``` For each input file (whether that's a single file or a directory of files), a corresponding results PDB and data file will be saved to this location. The PDB file contains only those residues within contact distance, and the data file contains two lists of isolated residues. 
- ```distance_cutoff``` (default 3 angstrom) An atom-to-atom distance measure between residues. Residue pairs within (<=) 
the cutoff distance are kept as interaction residues.
- ```renum_chains``` Renumbers the residues on chain X from residue N e.g., A:1 renumbers residues, begging from "1", on chain A. For both chains, e.g., A:1 B:5, renumbers chain A beginning with "1", and chain B, beginning with "5". This operation
ensure hydrogen atom caps are added to the correct FMO fragment. Problems will occur if a chain has disordered or missing residue numbers e.g., 3, 4, 6, 11, 12, or 1, 3, 4, 7, 5.
The code uses the sequential order of residue numbers to determine whether there is a break in the peptide backbone, and therefore, a hydrogen atom cap is required. 
- ```ignore_num_start_res``` Ignores all residues up to this residue e.g., ```5``` ensures the first four residues are ignored. Uses residue number. Renumber with ```renum_chains``` if necessary.  
- ```ignore_num_end_res``` Ignores all residues after this residue e.g., ```26``` ensures all residues after the 26th are ignored. Uses residue number. Renumber with ```renum_chains``` if necessary.  

For example, to isolate interaction residues within 3 angstrom, while renumbering the residues of chain B (from residue 1): 

```
contact_distance:
  file: C:\Users\Anthony\PyCharmProjects\HelicalFMO\structures\4auo_single_chains_AG.pdb
  output_folder: C:\Users\Anthony\PyCharmProjects\HelicalFMO\temp\
  distance_cutoff: 3
  renum_chains: B:1
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
- ```file``` The PDB input file. The file must contain only two chains, A and B. It's likely this file will have been generated having first performed ```contact_distance```.
- ```folder``` Path to a folder containing one or more PDB input files. This is an alternative choice to the single input ```file``` parameter.
- ```output_folder``` For each input file (whether that's a single file or a directory of files), new capped PDB files will be saved to this location.

For example, we cap a file refined to only include residues on both chains within a 3 angstrom interchain cutoff distance. A new file is saved to the ```output_folder``` location.  

```
cpa:
    file: C:\Users\Anthony\PyCharmProjects\HelicalFMO\temp\contacts_0.pdb 
    output_folder: C:\Users\Anthony\PyCharmProjects\HelicalFMO\temp\
```

---

## Generating GAMESS FMO input files

This code builds a GAMESS input file for FMO analysis. All residues in ```file``` will be included. It's essential to run ```cap``` before (see above) 
to prepare the PDB file. This ensures the appropriate residues are included from a well-formed PDB file; residues IDs start from 1, and peptides are delimited by Chain A and B. 

The GAMESS input file is stored in ```output_folder```.

### Command line parameters
- ```file``` The path a PDB input file.
- ```folder``` The folder path to multiple PDB input files.
- ```output_folder``` The folder path to the generated GAMESS FMO input file. 
- ```basis``` The basis set added to the GAMESS input file. Options are ```STO-3G``` (default) and ```"6-31G*"```. They represent low and reasonable accuracy, respectively.  
- ```theory``` The level of quantum theory added to the GAMESS input file. Options are ```HF``` (default) and ```MP2```.

For example, to generate an FMO input file with the HF/STO-3G level of theory and basis set:

```
fmo:
  file: C:\Users\Anthony\PyCharmProjects\HelicalFMO\temp\contacts_0.pdb
  output_folder: C:\Users\Anthony\PyCharmProjects\HelicalFMO\temp
  basis: STO-3G
  theory: HF
```
