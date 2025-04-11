import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.base import AnalysisFromFunction
from MDAnalysis.lib.transformations import rotation_matrix
import os

from src.controllers.controller import Controller


class RotationController(Controller):

    def __init__(self, input_file: str, output_folder: str, rotation_step: int):
        super().__init__()

        self.input_file = input_file
        self.output_folder = output_folder
        self.rotation_step = rotation_step

    def validate_inputs(self):
        pass

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

        # Rotation increment in degrees
        step_deg = self.rotation_step
        step_rad = np.deg2rad(step_deg)

        # Full nested loop: rotate B 360° per A 20° rotation
        structure_id = 0
        for i in range(0, 360, step_deg):
            # Rotate Chain A
            self.__rotate_selection(chainA, step_rad, axis_A, center_A)

            for j in range(0, 360, step_deg):
                # Rotate Chain B
                self.__rotate_selection(chainB, step_rad, axis_B, center_B)

                # Save structure
                outname = f"{self.output_folder}/rot_{structure_id:03d}.pdb"
                with mda.Writer(outname, u.atoms.n_atoms) as w:
                    w.write(u.atoms)

                # Undo Chain B rotation
                self.__rotate_selection(chainB, -step_rad, axis_B, center_B)

                structure_id += 1

            # Undo Chain A rotation before next A rotation
            self.__rotate_selection(chainA, -step_rad, axis_A, center_A)

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
