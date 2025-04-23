from sys import stdout
from typing import Dict, Union

from PeptideBuilder import Geometry
import PeptideBuilder
import Bio.PDB
import MDAnalysis as mda
import numpy as np
from openmm.app import *
from openmm import *
import openmm.unit as unit

from src.logger_config import get_logger
from src.controllers.controller import Controller
from src.models import pdb_writer


class GenHelixController(Controller):

    def __init__(self):
        super().__init__()
        self.logger = get_logger(__name__)

        self.sequence_A = None
        self.sequence_B = None
        self.output_file = None

    def validate_controller(self, config_section: Dict[str, Union[str, int, float, bool]]):
        self.sequence_A = config_section.get("seq_a")
        self.sequence_B = config_section.get("seq_b")
        self.output_file = config_section.get("output_file")

        if not self.sequence_A and not self.sequence_B:
            return False
        if not self.output_file:
            return False

        return True

    #def validate_inputs(self):
    #    if not self.sequence_A and not self.sequence_B:
    #        return False
    #    if not self.output_file:
    #        return False

    #    return True

    def run_controller(self):
        geo_A_list = []
        for aa in self.sequence_A:
            geo = Geometry.geometry(aa)
            geo.phi = -60
            geo.psi_im1 = -40
            geo_A_list.append(geo)

        structure_A = PeptideBuilder.initialize_res(geo_A_list[0])
        for geo in geo_A_list[1:]:
            PeptideBuilder.add_residue(structure_A, geo)
        PeptideBuilder.add_terminal_OXT(structure_A)

        if self.sequence_B:
            #make a heterodimer
            geo_B_list = []
            for aa in self.sequence_B:
                geo = Geometry.geometry(aa)
                geo.phi = -60
                geo.psi_im1 = -40
                geo_B_list.append(geo)

            structure_B = PeptideBuilder.initialize_res(geo_B_list[0])
            for geo in geo_B_list[1:]:
                PeptideBuilder.add_residue(structure_B, geo)
            PeptideBuilder.add_terminal_OXT(structure_B)
        else:
            #make a homodimer
            structure_B = PeptideBuilder.initialize_res(geo_A_list[0])
            for geo in geo_A_list[1:]:
                PeptideBuilder.add_residue(structure_B, geo)
            PeptideBuilder.add_terminal_OXT(structure_B)

        # save each structure to a temp file
        outA = Bio.PDB.PDBIO()
        outA.set_structure(structure_A)
        outA.save("../temp/helixA.pdb")

        outB = Bio.PDB.PDBIO()
        outB.set_structure(structure_B)
        outB.save("../temp/helixB.pdb")

        # calculate the average radius of each helix; keep the longest
        max_A, radius_A = self.__calculate_average_helix_radius("../temp/helixA.pdb")
        max_B, radius_B = self.__calculate_average_helix_radius("../temp/helixB.pdb")

        # load both structures MDAnalysis.
        uniA = mda.Universe("../temp/helixA.pdb")
        uniB = mda.Universe("../temp/helixB.pdb")

        # displace helix B by the longest maximum radius
        selA = uniA.select_atoms("all")
        selB = uniB.select_atoms("all")

        #this does not work
        self.__align_principal_axis_to_z(selA)
        self.__align_principal_axis_to_z(selB)

        # Translate helixB
        selB.translate([max(max_A*2, max_B*2)+13.8, 0.0, 0.0])

        # collect properties to rebuild the universe
        new_positions_A = []
        new_names_A = []
        new_resids_A = []
        new_resnames_A = []

        residuesA = selA.residues
        residue_atom_counts_A = []
        for i, residue in enumerate(residuesA):
            for atom in residue.atoms:

                new_names_A.append(atom.name)
                new_positions_A.append(atom.position)
            new_resnames_A.append(residue.resname)
            new_resids_A.append(i+1)
            residue_atom_counts_A.append(len(residue.atoms))

        new_positions_B = []
        new_names_B = []
        new_resids_B = []
        new_resnames_B = []

        residuesB = selB.residues
        residue_atom_counts_B = []
        for i, residue in enumerate(residuesB):
            for atom in residue.atoms:
                new_names_B.append(atom.name)
                new_positions_B.append(atom.position)
            new_resnames_B.append(residue.resname)
            new_resids_B.append(i+1)
            residue_atom_counts_B.append(len(residue.atoms))

        # add chain Label A


        atom_resindex_A = np.concatenate([
            np.full(count, i, dtype=int) for i, count in enumerate(residue_atom_counts_A)
        ])
        new_atoms_u_A = mda.Universe.empty(n_atoms=len(new_positions_A),
                                           n_residues=len(new_resnames_A),
                                           atom_resindex=atom_resindex_A,
                                           trajectory=True)
        new_atoms_u_A.add_TopologyAttr("name", new_names_A)
        new_atoms_u_A.add_TopologyAttr("resid", new_resids_A)
        new_atoms_u_A.add_TopologyAttr("resname", new_resnames_A)
        new_atoms_u_A.add_TopologyAttr("chainID")
        new_atoms_u_A.atoms.chainIDs = "A"
        new_atoms_u_A.atoms.positions = np.array(new_positions_A)

        # add chain Label B
        atom_resindex_B = np.concatenate([
            np.full(count, i, dtype=int) for i, count in enumerate(residue_atom_counts_B)
        ])
        new_atoms_u_B = mda.Universe.empty(n_atoms=len(new_positions_B),
                                           n_residues=len(new_resnames_B),
                                           atom_resindex=atom_resindex_B,
                                           trajectory=True)
        new_atoms_u_B.add_TopologyAttr("name", new_names_B)
        new_atoms_u_B.add_TopologyAttr("resid", new_resids_B)
        new_atoms_u_B.add_TopologyAttr("resname", new_resnames_B)
        new_atoms_u_B.add_TopologyAttr("chainID")
        new_atoms_u_B.atoms.chainIDs = "B"
        new_atoms_u_B.atoms.positions = np.array(new_positions_B)

        # combine and save to a file
        merged = mda.Merge(new_atoms_u_A.atoms, new_atoms_u_B.atoms)
        #merged.trajectory.unitcell = [100.0, 100.0, 100.0, 90.0, 90.0, 90.0]
        pdb_writer.write_fragments_pdb("../temp/initial_dimer.pdb", merged)

        # add the cryst1_line to the top of the file
        cryst1_line = "CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1           1\n"
        with open("../temp/initial_dimer.pdb", "r") as pdb_file:
            pdb_content = pdb_file.readlines()

        # Find the CRYST1 line and replace it
        for i, line in enumerate(pdb_content):
            if line.startswith("CRYST1"):
                pdb_content[i] = cryst1_line
                break

        with open("../temp/initial_dimer.pdb", "w") as modified_file:
            modified_file.writelines(pdb_content)

        self.__move_helix_to_center_of_unit_cell("../temp/initial_dimer.pdb")

        # load into OpenMM and run an energy minimisation
        self.__minimise_structure("../temp/initial_dimer.pdb")

        # when the program starts check whether ./temp exists. If it doesn't, make it

    def __move_helix_to_center_of_unit_cell(self, pdb_file):

        with open(pdb_file, "r") as file:
            for line in file:
                if line.startswith("CRYST1"):
                    # Extract the six cell dimensions from the CRYST1 line
                    cell_params = line[6:54].strip().split()
                    cell = np.array([float(val) for val in cell_params])
                    break
            else:
                raise ValueError("No CRYST1 line found in the PDB file")

        # Load the universe
        u = mda.Universe(pdb_file)

        # Select atoms of the helix (example: using Cα atoms)
        helix = u.select_atoms("all")

        # Parse the unit cell dimensions from the string
        #cell = np.array([float(val) for val in cell_string.split()])

        # Compute the center of the unit cell (only using the first three dimensions)
        unit_cell_center = cell[:3] / 2

        # Get the center of geometry of the helix (average position of all Cα atoms)
        helix_center = helix.center_of_geometry()

        # Calculate the shift needed to move the helix to the center of the unit cell
        shift_vector = unit_cell_center - helix_center

        # Apply the shift to the atom positions
        helix.positions += shift_vector

        # Optionally, save the new structure to a file
        with mda.Writer(pdb_file, u.atoms.n_atoms) as w:
            w.write(u)

        return helix.positions

    def __minimise_structure(self, input_file: str):
        pdb = PDBFile(input_file)
        forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
        modeller = Modeller(pdb.topology, pdb.positions)
        modeller.addHydrogens(forcefield)  # ← ADD correct hydrogens
        with open("../temp/dimer_h.pdb", "w") as file:
            PDBFile.writeFile(modeller.topology, modeller.positions, file)
        pdb = PDBFile('../temp/dimer_h.pdb')
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1 * unit.nanometer,
                                         constraints=HBonds)
        integrator = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
        simulation = Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy()

        minimized_positions = simulation.context.getState(getPositions=True).getPositions()
        simulation.reporters.append(PDBReporter('output.pdb', 1000))
        simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))
        simulation.step(10000)
        with open(self.output_file, 'w') as f:
            PDBFile.writeFile(simulation.topology, minimized_positions, f)

    def __align_principal_axis_to_z(self, selection):
        # Calculate the principal axis (longest moment of inertia axis)
        #principal_axis = selection.principal_axes()[0]
        #principal_axis /= np.linalg.norm(principal_axis)

        com = selection.center_of_mass()

        # Calculate the axis of the helix by fitting a line through the backbone atoms
        backbone = selection.select_atoms("backbone")  # Or use all atoms if needed
        backbone_positions = backbone.positions
        p0 = backbone_positions[0]
        p1 = backbone_positions[-1]

        # Helix axis is the vector between the first and last backbone atoms
        helix_axis = p1 - p0
        helix_axis /= np.linalg.norm(helix_axis)  # Normalize the helix axis

        # Define the z-axis
        z_axis = np.array([0.0, 0.0, 1.0])

        # Compute rotation axis (cross product) and angle
        rot_axis = np.cross(helix_axis, z_axis)
        rot_axis_norm = np.linalg.norm(rot_axis)

        if rot_axis_norm < 1e-6:
            # Already aligned with z-axis (or exactly opposite)
            if np.dot(principal_axis, z_axis) < 0:
                # 180° rotation around any perpendicular axis
                rot_axis = np.array([1.0, 0.0, 0.0])
                angle = np.pi
            else:
                return  # Already aligned
        else:
            rot_axis /= rot_axis_norm
            angle = np.arccos(np.clip(np.dot(helix_axis, z_axis), -1.0, 1.0))

        # Build the rotation matrix
        #R = rotation_matrix(angle, rot_axis)

        R = self.__rotation_matrix(angle, rot_axis)

        # Apply the rotation
        selection.rotate(R[:3, :3])

    def __rotation_matrix(self, angle, axis):
        """
        Creates a rotation matrix for a given angle and axis.
        This uses the Rodrigues' rotation formula.
        """
        axis = axis / np.linalg.norm(axis)  # Normalize axis

        cos_angle = np.cos(angle)
        sin_angle = np.sin(angle)
        x, y, z = axis

        # Rotation matrix formula
        R = np.array([
            [cos_angle + x ** 2 * (1 - cos_angle), x * y * (1 - cos_angle) - z * sin_angle,
             x * z * (1 - cos_angle) + y * sin_angle],
            [y * x * (1 - cos_angle) + z * sin_angle, cos_angle + y ** 2 * (1 - cos_angle),
             y * z * (1 - cos_angle) - x * sin_angle],
            [z * x * (1 - cos_angle) - y * sin_angle, z * y * (1 - cos_angle) + x * sin_angle,
             cos_angle + z ** 2 * (1 - cos_angle)]
        ])

        return R

    def __calculate_average_helix_radius(self, pdb_file) -> (float, float):
        u = mda.Universe(pdb_file)
        helix = u.select_atoms("name CA")

        # get positions and centre them
        positions = helix.positions
        center = helix.center_of_geometry()
        positions_centered = positions - center

        # compute covariance matrix of positions
        cov = np.cov(positions_centered.T)

        # get eigenvectors (principal axes) and eigenvalues
        eigvals, eigvecs = np.linalg.eigh(cov)

        # principal axis is the eigenvector corresponding to the largest eigenvalue
        principal_axis = eigvecs[:, np.argmax(eigvals)]
        principal_axis = principal_axis / np.linalg.norm(principal_axis)

        # compute radial distances from each point to the principal axis
        radii = []
        for pos in positions_centered:
            # project position onto axis
            projection_length = np.dot(pos, principal_axis)
            projection = projection_length * principal_axis
            radial_vector = pos - projection
            radius = np.linalg.norm(radial_vector)
            radii.append(radius)

        max_radius = np.max(radii)
        average_radius = np.median(radii)

        return max_radius, average_radius







#------------------------------------
#For an alpha helix

#from PeptideBuilder import Geometry
#import PeptideBuilder

# create a peptide consisting of 6 glycines
#geo = Geometry.geometry("G")
#geo.phi = -60
#geo.psi_im1 = -40
#structure = PeptideBuilder.initialize_res(geo)
#for i in range(5):
#    PeptideBuilder.add_residue(structure, geo)
# add terminal oxygen (OXT) to the final glycine
#PeptideBuilder.add_terminal_OXT(structure)

#import Bio.PDB

#out = Bio.PDB.PDBIO()
#out.set_structure(structure)
#out.save("example.pdb")