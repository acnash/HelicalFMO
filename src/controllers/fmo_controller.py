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
        chain_A = self.universe.select_atoms("chain A")
        chain_B = self.universe.select_atoms("chain B")

        self.__work_on_chain(chain_A)
        self.__work_on_chain(chain_B)

    def __work_on_chain(self, chain: mda.AtomGroup):
        residues = [(res.resname, res.resid) for res in chain.residues]

        print("bob")

    def __add_gamess_instructions(self):
        pass