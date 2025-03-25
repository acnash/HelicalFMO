import MDAnalysis as mda

from src.controllers.controller import Controller
from src.logger_config import get_logger


class FMOController(Controller):

    def __init__(self, pdb_file, output_gamess_file):
        super().__init__()

        self.logger = get_logger(__name__)
        self.pdb_file = pdb_file
        self.output_gamess_file = output_gamess_file
        self.universe = mda.Universe(self.pdb_file)

    def run_controller(self):
        chain_A = self.universe.select_atoms("chain A")
        chain_B = self.universe.select_atoms("chain B")

        self.__work_on_chain(chain_A)
        self.__work_on_chain(chain_B)

    def __work_on_chain(self, chain: mda.AtomGroup):
        pass

    def __add_gamess_instructions(self):
        pass