from typing import List

import MDAnalysis as mda
from MDAnalysis import AtomGroup
import warnings

from src.logger_config import get_logger

warnings.filterwarnings("ignore")  # turn off MDAnalysis warning when writing PDB files


def write_fragments_pdb(file_name: str, close_atoms) -> None:
    # Write the selected atoms to a new PDB file
    with mda.Writer(file_name) as writer:  #, close_atoms.n_atoms) as writer:
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


def write_pdb_to_psi4in(u: mda.Universe, output_file: str, charge: int) -> None:
    element_atomic_number = {
        "H": 1,  # Hydrogen
        "He": 2,  # Helium
        "Li": 3,  # Lithium
        "Be": 4,  # Beryllium
        "B": 5,  # Boron
        "C": 6,  # Carbon
        "N": 7,  # Nitrogen
        "O": 8,  # Oxygen
        "F": 9,  # Fluorine
        "Ne": 10,  # Neon
        "Na": 11,  # Sodium
        "Mg": 12,  # Magnesium
        "Al": 13,  # Aluminum
        "Si": 14,  # Silicon
        "P": 15,  # Phosphorus
        "S": 16,  # Sulfur
        "Cl": 17  # Chlorine
    }

    xyz_list = []
    for atom in u.atoms:
        element = atom.name[:2] if atom.name[:2] in element_atomic_number else atom.name[0]  # Extract element symbol from atom name
        #atomic_number = element_atomic_number[element]  # Placeholder (MDAnalysis doesn't store atomic numbers)
        xyz_entry = f"{element} {atom.position[0]:.3f} {atom.position[1]:.3f} {atom.position[2]:.3f} \n"
        xyz_list.append(xyz_entry)

    with open(output_file, "w") as file:
        file.write("molecule {\n")
        file.write(str(charge) + " " + "1\n")
        file.writelines(xyz_list)
        file.write("}\n")
        file.write("\n")
        file.write("memory 8 GB\n")
        file.write("\n")
        file.write("set basis sto-3g\n")
        file.write("\n")
        file.write("gradient('HF')\n")


def update_pdb_from_xyz():
    pass
