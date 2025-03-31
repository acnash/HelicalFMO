from typing import List

import MDAnalysis as mda


def check_chains_in_pdb(u: mda.Universe) -> bool:
    # Load the PDB file into a Universe object

    # Check if Chain A and Chain B are present
    chain_A = set(u.select_atoms('chainid A'))  # Get the chains from atoms of chain A
    chain_B = set(u.select_atoms('chainid B'))  # Add chains from chain B to the set

    # If both 'A' and 'B' are in the set, return True
    return chain_A and chain_B


def collect_two_chains(u: mda.Universe, ignore_num_start_res, ignore_num_end_res) -> mda.Universe:
    u.trajectory[0]  # Select the first model/frame

    chain_a_first_model = u.select_atoms("chainid A")
    chain_b_first_model = u.select_atoms("chainid B")

    chain_a_residues = list(chain_a_first_model.residues)
    chain_b_residues = list(chain_b_first_model.residues)

    #if ignore_num_start_res == 0 and ignore_num_end_res == 0:
    #    chain_a_selected = chain_a_residues
    #    chain_b_selected = chain_b_residues
    if ignore_num_start_res > 0 and ignore_num_end_res == 0:
        chain_a_selected = chain_a_residues[ignore_num_start_res:]
        chain_b_selected = chain_b_residues[ignore_num_start_res:]
    elif ignore_num_start_res > 0 and ignore_num_end_res > 0:
        chain_a_selected = chain_a_residues[ignore_num_start_res:-ignore_num_end_res]
        chain_b_selected = chain_b_residues[ignore_num_start_res:-ignore_num_end_res]
    elif ignore_num_start_res == 0 and ignore_num_end_res > 0:
        chain_a_selected = chain_a_residues[:-ignore_num_end_res]
        chain_b_selected = chain_b_residues[:-ignore_num_end_res]
    else:  #  ignore_num_start_res == 0 and ignore_num_end_res == 0:
        chain_a_selected = chain_a_residues
        chain_b_selected = chain_b_residues

    # Select atoms of the selected residues
    selected_atoms_a = u.select_atoms("chainid A and resid " + " ".join(str(res.resid) for res in chain_a_selected))
    selected_atoms_b = u.select_atoms("chainid B and resid " + " ".join(str(res.resid) for res in chain_b_selected))

    # Combine the selections (chain A + chain B)
    selected_atoms = selected_atoms_a + selected_atoms_b

    # Create a new Universe with only the selected atoms
    new_universe = mda.Universe.empty(len(selected_atoms))
    new_universe.atoms = selected_atoms.atoms

    return new_universe


def renumber_chain_resids(u: mda.Universe, chain_list: List[str]) -> mda.Universe:
    # Create a copy of the universe to modify
    new_u = mda.Universe.empty(n_atoms=len(u.atoms), trajectory=True)
    new_u.atoms = u.atoms  # Copy atom data (topology)
    new_u.load_new(u.filename)  # Load coordinates (pdb or other structure)

    # Loop through each specified chain
    for chain_data in chain_list:
        chain, resid = chain_data.split(":")
        resid = int(resid)
        # Select all residues from the given chain
        chain_residues = new_u.select_atoms(f"chainid {chain}").residues

        # Renumber residues starting from 1
        for new_resid, residue in enumerate(chain_residues, start=resid):
            residue.resid = new_resid

    return new_u

def sort_universe(u: mda.Universe) -> mda.Universe:
    sorted_atoms = sorted(u.atoms, key=lambda x: (x.segid, x.residue.resid, x.name))
    sorted_universe = mda.Universe(u.filename, atoms=sorted_atoms)
    return sorted_universe
