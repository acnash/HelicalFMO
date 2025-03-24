#
import os.path

from MDAnalysis.analysis import distances
import numpy as np

from controller import Controller
from src.models import pdb_cleaner, pdb_reader, pdb_writer


class ContactController(Controller):
    def __init__(self, file_location, folder_location, output_folder, distance_cutoff=8.0):
        super().__init__()

        self.file_location = file_location
        self.folder_location = folder_location
        self.distance_cutoff = distance_cutoff
        self.output_folder = output_folder

        self.universe_list = []

    def validate_inputs(self):
        if self.file_location:
            universe = pdb_reader.read_pdb_file(self.file_location)
            if pdb_cleaner.check_chains_in_pdb(universe):
                universe = pdb_cleaner.collect_two_chains(universe)
                self.universe_list.append(universe)
                print("Read and processed two chains from one model")
            else:
                return False

        else:
            universe_list = pdb_reader.read_pdb_folder(self.file_location)
            for universe in universe_list:
                if pdb_cleaner.check_chains_in_pdb(universe):
                    universe = pdb_cleaner.collect_two_chains(universe)
                    self.universe_list.append(universe)
                    print("Read and processed two chains from one model")
                else:
                    return False

        return True


    def run_controller(self):
        count = 0
        for universe in self.universe_list:
            close_residues = set()  # Set to store residues that are close enough

            chain_a_atoms = universe.select_atoms('chainid A')
            chain_b_atoms = universe.select_atoms('chainid B')

            # Loop through all atoms in chain A
            for atom_a in chain_a_atoms:
                # Calculate distances to all atoms in chain B
                distances_matrix = distances.distance_array(atom_a.position, chain_b_atoms.positions)

                # Find atoms in chain B within the cutoff distance
                close_atoms_b = chain_b_atoms[np.where(distances_matrix < self.distance_cutoff)[0]]

                # Add residues from chain B that are within the cutoff
                for atom_b in close_atoms_b:
                    close_residues.add(atom_b.residue)

            # Loop through all atoms in chain B
            for atom_b in chain_b_atoms:
                # Calculate distances to all atoms in chain A
                distances_matrix = distances.distance_array(atom_b.position, chain_a_atoms.positions)

                # Find atoms in chain A within the cutoff distance
                close_atoms_a = chain_a_atoms[np.where(distances_matrix <= self.distance_cutoff)[0]]

                # Add residues from chain A that are within the cutoff
                for atom_a in close_atoms_a:
                    close_residues.add(atom_a.residue)

            # Select atoms from the close residues (both chains A and B)
            close_atoms = universe.select_atoms('resid {} and (chainid A or chainid B)'.format(
                ' '.join([str(residue.resid) for residue in close_residues])
            ))

            file_name = "_".join("contacts", str(count), ".pdb")
            output_location = os.path.join(self.output_folder, file_name)
            pdb_writer.write_fragments_pdb(output_location, close_atoms)
            count = count + 1
