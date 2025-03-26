import os
from typing import List, Any

import MDAnalysis as mda

from src.controllers.controller import Controller
from src.logger_config import get_logger


class FMOController(Controller):

    def __init__(self, basis: str, theory: str):
        super().__init__()

        self.logger = get_logger(__name__)
        self.basis = basis
        self.theory = theory
        self.pdb_files = []
        self.universe_list = []
        self.output_location = ""
        self.system_charge_list = []


    def validate_inputs(self, pdb_file, pdb_folder, output_location):
        if pdb_file:
            self.pdb_files.append(pdb_file)
            self.universe_list.append(mda.Universe(pdb_file))
        else:
            self.pdb_files = [os.path.join(pdb_folder, f) for f in os.listdir(pdb_folder) if f.endswith(".pdb")]
            for file in self.pdb_files:
                self.universe_list.append(mda.Universe(file))

        self.output_location = output_location

    def run_controller(self):
        # operate at scale over all universe objects
        # an iteration results in an FMO GAMESS input file
        for universe in self.universe_list:
            chain_A_residues = universe.select_atoms("chainid A").residues
            chain_B_residues = universe.select_atoms("chainid B").residues

            chain_A_uni = self.__reorder_atoms_in_chain(chain_A_residues)
            chain_B_uni = self.__reorder_atoms_in_chain(chain_B_residues)

            chain_A_indat = self.__build_indat(chain_A_residues)
            chain_B_indat = self.__build_indat(chain_B_residues)

            chain_A_frag_boundary_ids = self.__find_fragment_boundaries(chain_A_uni)
            chain_B_frag_boundary_ids = self.__find_fragment_boundaries(chain_B_uni)

            chain_A_fmoxyz_list = self.__build_FMOXYZ(chain_A_uni)
            chain_B_fmoxyz_list = self.__build_FMOXYZ(chain_B_uni)

    def __build_indat(self, chain: mda.AtomGroup) -> List[str]:
        indat_list = []
        count = 1
        for residue in chain.residues:
            for atom in residue.atoms:
                indat_list.append(str(count))
                indat_list.append(",")
            indat_list.append("\n")
            count = count + 1

        return indat_list

    def __build_FMOXYZ(self, u: mda.Universe) -> List[str]:
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

        # Loop through all atoms in the Universe
        for atom in u.atoms:
            element = atom.name[:2] if atom.name[:2] in element_atomic_number else atom.name[0]  # Extract element symbol from atom name
            atomic_number = element_atomic_number[element]  # Placeholder (MDAnalysis doesn't store atomic numbers)
            xyz_entry = f"{element} {atomic_number} {atom.position[0]:.3f} {atom.position[1]:.3f} {atom.position[2]:.3f}"
            xyz_list.append(xyz_entry)

        return xyz_list

    def __find_fragment_boundaries(self, u: mda.Universe) -> List[tuple[Any, Any]]:
        residues = u.residues

        # List to store recorded atom indices
        c_ca_indices = []

        # Loop through residues except the last one
        for i in range(len(residues) - 1):
            current_res = residues[i]
            next_res = residues[i + 1]

            # Check if the next residue follows sequentially
            if current_res.resid + 1 == next_res.resid:
                # Get "C" and "CA" atom indices from the current residue
                c_atom = current_res.atoms.select_atoms("name C")
                ca_atom = current_res.atoms.select_atoms("name CA")

                # Record atom indices if found
                if len(c_atom) > 0 and len(ca_atom) > 0:
                    c_ca_indices.append((c_atom[0].index, ca_atom[0].index))

        # Print results
        #print("Recorded C and CA atom indices:", c_ca_indices)
        return c_ca_indices


    def __reorder_atoms_in_chain(self, residues: mda.AtomGroup) -> mda.Universe:

        # List to store reordered AtomGroups
        reordered_atom_groups = []

        # Iterate through each residue in Chain A
        for res in residues:
            atoms = res.atoms  # Get all atoms in the residue

            # Select CA and C atoms
            ca_atom = atoms.select_atoms("name CA")
            c_atom = atoms.select_atoms("name C")

            # Select all other atoms except CA and C
            other_atoms = atoms.select_atoms("not name CA C")

            # Merge atoms in new order (other atoms first, then CA and C)
            reordered_atom_group = other_atoms + ca_atom + c_atom

            # Store the reordered AtomGroup
            reordered_atom_groups.append(reordered_atom_group)

        # Merge all reordered AtomGroups into a new universe
        new_u = mda.Merge(*reordered_atom_groups)  # Unpacking AtomGroups

        return new_u

    def __add_gamess_instructions(self):
        instructions = []
        for system_charge in self.system_charge_list:

            if self.basis == "6-31G*":
                basis_str = "$BASIS GBASIS=N31 NGAUSS=6 NDFUNC=1 $END\n"
            else:
                basis_str = "$BASIS GBASIS=STO NGAUSS=3 $END\n"

            if self.theory == "HF":
                control_str = f"$CONTROL SCFTYP=RHF RUNTYP=ENERGY ICHARG={str(system_charge)} $END\n"
            else:
                control_str = f"$CONTROL SCFTYP=RHF RUNTYP=ENERGY MPLEVL=2 ICHARG={str(system_charge)} $END\n"

            memory_str = "$SYSTEM MEMORY=1000000 $END\n"
            instructions.append("".join([basis_str, control_str, memory_str]))


    def __write_hybrid_orbitals(self, basis_set: str) -> List[str]:
        basis_631G_star = []
        basis_631G_star.append("6-31G* 15 5")
        basis_631G_star.append("0 1   0.995621   0.028507   0.000000   0.000000   0.000000")
        basis_631G_star.append("     -0.013353   0.000000   0.000000  -0.000001  -0.002076")
        basis_631G_star.append("     -0.002076  -0.002076   0.000000   0.000000   0.000000")
        basis_631G_star.append("0 1   0.097077  -0.179454   0.178051   0.308393   0.125888")
        basis_631G_star.append("     -0.181305   0.093381   0.161741   0.066029   0.000102")
        basis_631G_star.append("     -0.020671   0.005300  -0.020773  -0.008476  -0.014681")
        basis_631G_star.append("0 1   0.097077  -0.179454  -0.356101  -0.000000   0.125888")
        basis_631G_star.append("     -0.181304  -0.186762  -0.000000   0.066030  -0.031058")
        basis_631G_star.append("0 1   0.097077  -0.179454   0.178051  -0.308393   0.125888")
        basis_631G_star.append("     -0.181304   0.093381  -0.161741   0.066029   0.000102")
        basis_631G_star.append("     -0.020672   0.005300   0.020773  -0.008476   0.014681")
        basis_631G_star.append("1 0  -0.097061   0.179413   0.000000   0.000000   0.377708")
        basis_631G_star.append("      0.181288   0.000000   0.000000   0.198114  -0.010489")
        basis_631G_star.append("     -0.010489   0.036248  -0.000000  -0.000000  -0.000000")

        basis_sto_3g = []
        basis_sto_3g.append("STO-3G 5 5")
        basis_sto_3g.append("0 1   0.991926   0.038355   0.000000   0.000000  -0.000000")
        basis_sto_3g.append("0 1  -0.110715   0.313793  -0.233504  -0.404441  -0.165113")
        basis_sto_3g.append("1 0  -0.110716   0.313795  -0.000000  -0.000000   0.495336")
        basis_sto_3g.append("0 1  -0.110715   0.313793   0.467008   0.000000  -0.165112")
        basis_sto_3g.append("0 1  -0.110715   0.313793  -0.233504   0.404441  -0.165112")

        if basis_set == "STO-3G":
            return basis_sto_3g
        else:
            return basis_631G_star


