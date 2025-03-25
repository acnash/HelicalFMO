#
import os.path

from MDAnalysis.analysis import distances
import numpy as np

from src.controllers.controller import Controller
from src.models import pdb_cleaner, pdb_reader, pdb_writer
from src.logger_config import get_logger


class ContactController(Controller):
    def __init__(self, file_location, folder_location, output_folder, distance_cutoff):
        super().__init__()

        self.logger = get_logger(__name__)

        self.file_location = file_location
        self.folder_location = folder_location
        self.distance_cutoff = distance_cutoff
        self.output_folder = output_folder

        self.universe_list = []

    def validate_inputs(self, ignore_num_start_res, ignore_num_end_res):
        if self.file_location:
            universe = pdb_reader.read_pdb_file(self.file_location)
            if pdb_cleaner.check_chains_in_pdb(universe):
                universe = pdb_cleaner.collect_two_chains(universe, ignore_num_start_res, ignore_num_end_res)
                self.universe_list.append(universe)
                print("Read and processed two chains from one model")
                self.logger.info("Read and processed two chains from one model")
            else:
                return False

        else:
            universe_list = pdb_reader.read_pdb_folder(self.file_location)
            for universe in universe_list:
                if pdb_cleaner.check_chains_in_pdb(universe):
                    universe = pdb_cleaner.collect_two_chains(universe)
                    self.universe_list.append(universe)
                    print("Read and processed two chains from one model")
                    self.logger.info("Read and processed two chains from one model")
                else:
                    return False

        return True

    def run_controller(self):
        count = 0
        for universe in self.universe_list:
            close_residues = set()  # Set to store residues that are close enough

            chain_a_atoms = universe.select_atoms('chainid A')
            chain_b_atoms = universe.select_atoms('chainid B')

            # Compute pairwise distance matrix in one step
            distance_matrix = distances.distance_array(chain_a_atoms.positions, chain_b_atoms.positions)

            # Find pairs of atoms that are within the cutoff distance
            close_pairs = np.where(distance_matrix <= self.distance_cutoff)

            # Get the residues for the atoms in chain A
            close_residues.update(chain_a_atoms[close_pairs[0]].residues)
            # Get the residues for the atoms in chain B
            close_residues.update(chain_b_atoms[close_pairs[1]].residues)

            # Select atoms from the close residues (both chains A and B)
            close_atoms = universe.select_atoms('resid {} and (chainid A or chainid B)'.format(
                ' '.join(map(str, [res.resid for res in close_residues]))
            ))

            # Write the fragment atom records to a new PDB file
            pdb_file_name = "".join(["contacts_", str(count), ".pdb"])
            pdb_output_location = os.path.join(self.output_folder, pdb_file_name)
            pdb_writer.write_fragments_pdb(pdb_output_location, close_atoms)

            # Write to a text file the list of residues in contact distance
            fragment_u = pdb_reader.read_pdb_file(pdb_output_location)
            chain_a_atoms = fragment_u.select_atoms('chainid A')
            chain_b_atoms = fragment_u.select_atoms('chainid B')
            contact_data_file_name = "".join(["contacts_", str(count), ".txt"])
            contacts_output_location = os.path.join(self.output_folder, contact_data_file_name)
            pdb_writer.write_found_residue_contacts(contacts_output_location, [chain_a_atoms, chain_b_atoms], ["A", "B"])

            count = count + 1
