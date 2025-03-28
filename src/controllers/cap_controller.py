import MDAnalysis as mda
import numpy as np

from src.controllers.controller import Controller
from src.models import pdb_reader


class CapController(Controller):

    def __init__(self, input_file, input_directory):
        super().__init__()
        self.input_file = input_file
        self.input_directory = input_directory

        self.universe_list = []

    def validate_inputs(self):
        if self.input_file:
            universe = pdb_reader.read_pdb_file(self.input_file)
            self.universe_list.append(universe)
        else:
            self.universe_list = pdb_reader.read_pdb_folder(self.input_directory)

    def run_controller(self):
        for universe in self.universe_list:
            chain_A_residues = universe.select_atoms("chainid A").residues
            for i, residue in enumerate(chain_A_residues):
                resnum = residue.resid

                # if the residue had an immediate left neighbour (resid - 1) then we don't need to cap it.
                has_prev = i > 0 and chain_A_residues[i-1].resid == resnum - 1
                # if the residue had an immediate right neighbour (resid +1) then we don't need to cap it.
                has_next = i < len(chain_A_residues) - 1 and chain_A_residues[i + 1].resid == resnum + 1

                if not has_prev:
                    print("Process hydrogen cap on the N atom")
                    N = residue.atoms.select_atoms("name N")[0].position
                    C_alpha = residue.atoms.select_atoms("name CA")[0].position
                    bond_vector = C_alpha - N
                    bond_vector /= np.linalg.norm(bond_vector)
                    NH_bond_length = 1.0
                    H_position = N - NH_bond_length * bond_vector
                    resname = residue.resname
                    chain = residue.segid
                    resid = residue.resid

                    # The above is horribly wrong.
                if not has_next:
                    print("Process hydrogen cap on the C atom")