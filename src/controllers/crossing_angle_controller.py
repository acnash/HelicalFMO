import os
import glob
from typing import Dict, Union, List

import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import calc_bond_distance

from src.logger_config import get_logger
from src.controllers.controller import Controller


class CrossingAngleController(Controller):

    def __init__(self):
        super().__init__()
        self.logger = get_logger(__name__)

        self.output_folder = None
        self.crossing_angle = None
        self.angle_points = None

    def validate_controller(self, config_section: Dict[str, Union[str, int, float, bool]]):
        self.output_folder = config_section.get("output_folder")
        self.crossing_angle = config_section.get("crossing_angle")
        self.angle_points = config_section.get("angle_points")

        return True

    def run_controller(self):
        pdb_files = glob.glob(os.path.join(self.output_folder, "rot_*.pdb"))

        for pdb_file in pdb_files:
            print(f"Loading rotation file {pdb_file}")
            self.logger.info(f"Loading rotation file {pdb_file}")
            self.__set_crossing_angle(pdb_file)

    def __set_crossing_angle(self, pdb_file):
        helical_distances, helix_vector = self.__calculate_helix_pair_length(pdb_file)
        length_length = max(helical_distances)

        u = mda.Universe(pdb_file)

        # Select chains
        chainA = u.select_atoms("segid A")
        chainB = u.select_atoms("segid B")

        vector_A = chainA[-1].position - chainA[0].position
        vector_B = chainB[-1].position - chainB[0].position
        rotation_axis = np.cross(vector_A, vector_B)
        rotation_axis /= np.linalg.norm(rotation_axis)

        # Normalize the helix direction
        #helix_direction = helix_vector / helix_length

        #calculate the length of the helix (use the longest)
        for i in range(self.angle_points):
            u_copy = mda.Universe(pdb_file)
            chainB_copy = u_copy.select_atoms("segid B and backbone")

            # Find the position along the helix for this rotation center
            fraction = (i + 1) / (N + 1)  # Avoid endpoints exactly
            midpoint = start_pos + fraction * helix_vector

            # Define the rotation
            rot = R.from_rotvec(np.deg2rad(angle_deg) * rotation_axis)

            # Apply rotation to chain B
            for atom in chainB_copy:
                rel_pos = atom.position - midpoint
                rotated_pos = rot.apply(rel_pos)
                atom.position = rotated_pos + midpoint

            # Save the rotated structure
            output_filename = os.path.join(pdb_file.removesuffix(".pdb"), f"_{i}.pdb")
            with mda.Writer(output_filename, u_copy.atoms.n_atoms) as w:
                w.write(u_copy.atoms)

            print(f"Saved: {output_filename}")

    def __calculate_helix_pair_length(self, input_file: str) -> List:
        u = mda.Universe(input_file)

        chainA = u.select_atoms("segid A and backbone")
        chainB = u.select_atoms("segid B and backbone")

        first_atom_A = chainA.atoms[0]
        last_atom_A = chainA.atoms[-1]

        first_atom_B = chainB.atoms[0]
        last_atom_B = chainB.atoms[-1]

        distance_A = calc_bond_distance(first_atom_A.position, last_atom_A.position, box=None)
        distance_B = calc_bond_distance(first_atom_B.position, last_atom_B.position, box=None)

        if distance_A > distance_B:
            helix_vector = last_atom_A - first_atom_A
        else:
            helix_vector = last_atom_B - first_atom_B

        # return the helix vector of the longest helix of the pair
        return [distance_A, distance_B, helix_vector]
