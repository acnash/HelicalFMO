import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.transformations import rotation_matrix
from sys import stdout

from src.controllers.controller import Controller

from openmm.app import *
from openmm import *
import openmm.unit as unit

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

                self.__minimise_structure(outname)

                # Undo Chain B rotation
                #self.__rotate_selection(chainB, -step_rad, axis_B, center_B)

            # Undo Chain A rotation before next A rotation
            #self.__rotate_selection(chainA, -step_rad, axis_A, center_A)


    def __minimise_structure(self, input_file: str):
        pdb = PDBFile(input_file)
        forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
        #modeller = Modeller(pdb.topology, pdb.positions)
        #modeller.addHydrogens(forcefield)  # ← ADD correct hydrogens
        #with open("../temp/dimer_h.pdb", "w") as file:
        #    PDBFile.writeFile(modeller.topology, modeller.positions, file)
        #pdb = PDBFile('../temp/dimer_h.pdb')
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1 * unit.nanometer,
                                         constraints=HBonds)
        integrator = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.001 * unit.picoseconds)
        simulation = Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy()


        #simulation.reporters.append(PDBReporter('output.pdb', 100))
        #simulation.reporters.append(StateDataReporter(stdout, 100, step=True, potentialEnergy=True, temperature=True))
        print(f"Resolving potential side chain clashes on {input_file}")
        try:
            simulation.step(500)
        except:


        minimized_positions = simulation.context.getState(getPositions=True).getPositions()
        with open(input_file, 'w') as f:
            PDBFile.writeFile(simulation.topology, minimized_positions, f)

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
