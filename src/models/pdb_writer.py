from typing import List

import MDAnalysis as mda
from MDAnalysis import AtomGroup

from src.logger_config import get_logger


def write_fragments_pdb(file_name: str, close_atoms) -> None:
    # Write the selected atoms to a new PDB file
    with mda.Writer(file_name, close_atoms.n_atoms) as writer:
        writer.write(close_atoms)


def write_found_residue_contacts(file_name: str, chain_list: List[AtomGroup], chain_name_list: List[str]) -> None:
    logger = get_logger(__name__)
    with open(file_name, "w") as f:
        for chain, name in zip(chain_list, chain_name_list):
            f.write(f"Chain {name} Residues:\n")
            chain_resnames = chain.residues.resnames
            chain_resids = chain.residues.resids

            for resname, resid in zip(chain_resnames, chain_resids):
                f.write(f"{resname} {resid}\n")

    print(f"Residues written to {file_name}")
    logger.info(f"Residues written to {file_name}")
