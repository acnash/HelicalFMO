import os

import MDAnalysis as mda

from src.controllers.controller import Controller
from src.logger_config import get_logger


class FMOController(Controller):

    def __init__(self):
        super().__init__()

        self.logger = get_logger(__name__)
        self.pdb_files = []
        self.universe_list = []
        self.output_location = ""

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
        for universe in self.universe_list:
            chain_A_residues = universe.select_atoms("chainid A").residues
            chain_B_residues = universe.select_atoms("chainid B").residues

            self.__work_on_chain(chain_A_residues)
            self.__work_on_chain(chain_B_residues)

    def __work_on_chain(self, residues: mda.AtomGroup):

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

        # Save the reordered PDB file
        new_u.atoms.write("reordered_chainA.pdb")

        print("I'm here....")




    def __add_gamess_instructions(self):
        pass
