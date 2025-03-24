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
    # Load the PDB file into a Universe object

    # Select the first model (first frame)
    u.trajectory[0]  # Select the first model/frame

    # Now select all atoms from chain A and chain B in the first model
    chain_a_first_model = u.select_atoms("chainid A")  # Select atoms from chain A in the first model
    chain_b_first_model = u.select_atoms("chainid B")  # Select atoms from chain B in the first model

    # Get the residues of each chain
    chain_a_residues = list(chain_a_first_model.residues)
    chain_b_residues = list(chain_b_first_model.residues)

    # Exclude the first and last five residues from both chains by their positions in the residue list
    if ignore_num_start_res == 0 and ignore_num_end_res == 0:
        chain_a_selected = chain_a_residues
        chain_b_selected = chain_b_residues
    elif ignore_num_start_res > 0 and ignore_num_end_res == 0:
        chain_a_selected = chain_a_residues[ignore_num_start_res:]  # Skip the first and last five residues
        chain_b_selected = chain_b_residues[ignore_num_start_res:]  # Skip the first and last five residues
    elif ignore_num_start_res > 0 and ignore_num_end_res > 0:
        chain_a_selected = chain_a_residues[ignore_num_start_res:-ignore_num_end_res]  # Skip the first and last five residues
        chain_b_selected = chain_b_residues[ignore_num_start_res:-ignore_num_end_res]  # Skip the first and last five residues
    elif ignore_num_start_res == 0 and ignore_num_end_res > 0:
        chain_a_selected = chain_a_residues[:-ignore_num_end_res]  # Skip the first and last five residues
        chain_b_selected = chain_b_residues[:-ignore_num_end_res]  # Skip the first and last five residues

    # Select atoms of the selected residues
    selected_atoms_a = u.select_atoms("chainid A and resid " + " ".join(str(res.resid) for res in chain_a_selected))
    selected_atoms_b = u.select_atoms("chainid B and resid " + " ".join(str(res.resid) for res in chain_b_selected))

    # Combine the selections (chain A + chain B)
    selected_atoms = selected_atoms_a + selected_atoms_b

    # Create a new Universe with only the selected atoms
    new_universe = mda.Universe.empty(len(selected_atoms))
    new_universe.atoms = selected_atoms.atoms

    # Add the selected atoms to the new Universe
    #new_universe.add_TopologyAttr('atoms', selected_atoms.atoms)
    #new_universe.add_TopologyAttr('residues', selected_atoms.residues)

    return new_universe
