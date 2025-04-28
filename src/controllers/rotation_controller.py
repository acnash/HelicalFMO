import MDAnalysis as mda
import numpy as np
from typing import Dict, Union
from MDAnalysis.lib.transformations import rotation_matrix

from src.controllers.controller import Controller
from src.controllers.crossing_angle_controller import CrossingAngleController
from src.logger_config import get_logger

from openmm import *

from src.models import pdb_writer


class RotationController(Controller):

    def __init__(self):
        super().__init__()
        self.logger = get_logger(__name__)

        self.input_file = None
        self.output_folder = None
        self.rotation_step = None
        self.both_helices = None
        self.cross_angle = None
        self.perform_crossing_angle = False
        self.crossing_angle_controller = None

    def validate_controller(self, config_section: Dict[str, Union[str, int, float, bool]]):
        self.input_file = config_section.get("file")
        self.output_folder = config_section.get("output_folder")
        self.rotation_step = config_section.get("rotation_angle")
        self.both_helices = config_section.get("both_helices")
        self.perform_crossing_angle = config_section.get("crossing_angle") is not None

        if not os.path.exists(self.input_file):
            print(f"Error: input_file {self.input_file} does not exist.")
            self.logger.error(f"Error: input_file {self.input_file} does not exist.")
            return False

        if self.perform_crossing_angle:
            self.crossing_angle_controller = CrossingAngleController()
            is_valid = self.crossing_angle_controller.validate_controller(config_section)
            if not is_valid:
                return False

        return True

    def run_controller(self):
        u = mda.Universe(self.input_file)
        chainA = u.select_atoms("segid A")
        chainB = u.select_atoms("segid B")

        # Get axes
        axis_A = self.__get_principal_axis(chainA)
        axis_B = self.__get_principal_axis(chainB)

        # Define rotation center (e.g. center of each helix)
        center_A = chainA.center_of_geometry()
        center_B = chainB.center_of_geometry()

        # Full nested loop: rotate B 360° per A 20° rotation
        for i in range(0, 360, self.rotation_step):
            # always reset to the original file before rotating further
            u = mda.Universe(self.input_file)
            chainA = u.select_atoms("segid A")
            chainB = u.select_atoms("segid B")

            step_rad = np.deg2rad(i)
            # Rotate Chain A
            self.__rotate_selection(chainA, step_rad, axis_A, center_A)

            for j in range(0, 360, self.rotation_step):
                # Rotate Chain B
                step_rad = np.deg2rad(j)
                self.__rotate_selection(chainB, step_rad, axis_B, center_B)

                # Save structure
                outname = os.path.join(self.output_folder, f"rot_{i}_{j}.pdb")
                with mda.Writer(outname, u.atoms.n_atoms) as w:
                    w.write(u.atoms)

                resolved = False
                while not resolved:
                    try:
                        print(f"Attempting to energy minimise {outname}.")
                        self.logger.info(f"Attempting to energy minimise {outname}.")
                        self._minimise_structure(outname, outname)
                        resolved = True
                    except:
                        print(f"Warning: Energy minimisation of {outname} failed. Adjusting the position of helix B by 0.1 A.")
                        self.logger.warning(f"Energy minimisation of {outname} failed. Adjusting the position of helix B by 0.1 A.")
                        self.__adjust_interhelical_distance(outname)

        print(f"Finished building rotated helices.")
        self.logger.info(f"Finished building rotated helices.")

        #now I adjust for crossing angle
        if self.perform_crossing_angle:
            print(f"Making adjustments to the crossing angle.")
            self.logger.info(f"Making adjustments to the crossing angle.")
            self.crossing_angle_controller.run_controller()


    def __adjust_interhelical_distance(self, input_file):
        u = mda.Universe(input_file)

        selA = u.select_atoms("segid A")
        selB = u.select_atoms("segid B")
        selB.translate([0.1, 0.0, 0.0])

        merged = mda.Merge(selA.atoms, selB.atoms)

        pdb_writer.write_fragments_pdb(input_file, merged)

        # reconstruct the unit cell information
        cryst1_line = "CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1           1\n"
        with open(input_file, "r") as pdb_file:
            pdb_content = pdb_file.readlines()

        # Find the CRYST1 line and replace it
        for i, line in enumerate(pdb_content):
            if line.startswith("CRYST1"):
                pdb_content[i] = cryst1_line
                break

        with open(input_file, "w") as modified_file:
            modified_file.writelines(pdb_content)

    def __get_principal_axis(self, sel):
        positions = sel.positions - sel.center_of_geometry()
        cov = np.cov(positions.T)
        eigvals, eigvecs = np.linalg.eigh(cov)
        axis = eigvecs[:, np.argmax(eigvals)]
        return axis / np.linalg.norm(axis)

    def __rotate_selection(self, selection, angle_rad, axis, center):
        selection.translate(-center)
        R = rotation_matrix(angle_rad, axis)[:3, :3]
        selection.rotate(R)
        selection.translate(center)
