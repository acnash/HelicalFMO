# HelicalFMO
Fragment Molecular Orbital analysis for interhelical transmembrane protein interactions

This code performs several operations:

1. Rotates helical transmembrane dimers and writes to PDB files.
2. Writes nearest residue-residue atomic records to PDF files. 
3. Generates GAMESS input file for FMO analysis based on nearest residue-residue PDB files.
4. Performs 1 to 3 at scale.

## Calculating protein-protein interaction residues

For example, to calculate interaction residues within 3 angstrom, while renumbering the residues of chain B (from residue 1): 
```
--file C:\Users\Anthony\PyCharmProjects\HelicalFMO\structures\4auo_single_chains_AB.pdb --mode contact_distance --output_folder C:\Users\Anthony\PyCharmProjects\HelicalFMO\temp\ --distance_cutoff 3 --renum_chains B:1
```
