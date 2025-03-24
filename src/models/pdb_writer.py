import MDAnalysis as mda

def write_fragments_pdb(file_name: str, close_atoms) -> None:
    # Write the selected atoms to a new PDB file
    with mda.Writer(file_name, close_atoms.n_atoms) as writer:
        writer.write(close_atoms)