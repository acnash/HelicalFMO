import MDAnalysis as mda


def check_chains_in_pdb(file_path: str) -> bool:
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
    u = mda.Universe(file_path)

    # Check if Chain A and Chain B are present
    chains = set(u.select_atoms("chain A").chains)  # Get the chains from atoms of chain A
    chains.update(u.select_atoms("chain B").chains)  # Add chains from chain B to the set

    # If both 'A' and 'B' are in the set, return True
    return 'A' in chains and 'B' in chains


def collect_two_chains(file_path: str) -> mda.Universe:
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
    u = mda.Universe(file_path)

    # Select the first model (first frame)
    u.trajectory[0]  # Select the first model/frame

    # Now select all atoms from chain A and chain B in the first model
    chain_a_first_model = u.select_atoms("chain A")  # Select atoms from chain A in the first model
    chain_b_first_model = u.select_atoms("chain B")  # Select atoms from chain B in the first model

    # Combine the selections (chain A + chain B)
    selected_atoms = chain_a_first_model + chain_b_first_model

    new_universe = mda.Universe.empty(len(selected_atoms))

    # Add the selected atoms to the new Universe
    new_universe.add_TopologyAttr('atoms', selected_atoms.atoms)
    new_universe.add_TopologyAttr('residues', selected_atoms.residues)

    return new_universe
