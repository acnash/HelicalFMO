import MDAnalysis as mda


def check_chains_in_pdb(u: mda.Universe) -> bool:
    """
    Checks if both Chain A and Chain B are present in the PDB file.

    Parameters:
    -----------
    file_path : str
        Path to the PDB file.

    Returns:
    --------
    bool
        True if both Chain A and Chain B are present in the PDB file, False otherwise.
    """
    # Load the PDB file into a Universe object

    # Check if Chain A and Chain B are present
    chain_A = set(u.select_atoms('chainid A'))  # Get the chains from atoms of chain A
    chain_B = set(u.select_atoms('chainid B'))  # Add chains from chain B to the set

    # If both 'A' and 'B' are in the set, return True
    return chain_A and chain_B


def collect_two_chains(u: mda.Universe, ignore_num_start_res, ignore_num_end_res) -> mda.Universe:
    """
    Loads a multi-model PDB file and retains all atoms from chain A and chain B in the first model.

    Parameters:
    -----------
    file_path : str
        Path to the multi-model PDB file.

    Returns:
    --------
    mda.Universe
        MDAnalysis Universe object containing only atoms from chain A and chain B in the first model.
    """


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
