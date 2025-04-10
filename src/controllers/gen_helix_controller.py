import MDAnalysis
from sys import stdout
from src.controllers.controller import Controller
from PeptideBuilder import Geometry
import PeptideBuilder
import Bio.PDB
import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.transformations import rotation_matrix
from openmm.app import *
from openmm import *
import openmm.unit as unit

from src.models import pdb_writer


class GenHelixController(Controller):

    def __init__(self, sequence_A, sequence_B, output_file):
        super().__init__()

        self.sequence_A = sequence_A
        self.sequence_B = sequence_B
        self.output_file = output_file

    def validate_inputs(self):
        if not self.sequence_A and not self.sequence_B:
            return False
        if not self.output_file:
            return False

        return True

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
        _, radius_A = self.__calculate_average_helix_radius("../temp/helixA.pdb")
        _, radius_B = self.__calculate_average_helix_radius("../temp/helixB.pdb")

        # load both structures MDAnalysis.
        uniA = mda.Universe("../temp/helixA.pdb")
        uniB = mda.Universe("../temp/helixB.pdb")

        # displace helix B by the longest maximum radius
        selA = uniA.select_atoms("all")
        selB = uniB.select_atoms("all")

        self.__align_principal_axis_to_z(selA)
        self.__align_principal_axis_to_z(selB)

        # Translate helixB
        selB.translate([max(radius_A, radius_B), 0.0, 0.0])

        new_positions_A = []
        new_names_A = []
        new_resids_A = []
        new_resnames_A = []

        residuesA = uniA.select_atoms("all").residues
        for i, residue in enumerate(residuesA):
            for atom in residue:
                new_names_A.append(atom.name)
                new_positions_A.append(atom.positions)
            new_resnames_A.append(residue.resname)
            new_resids_A.append(i+1)

        new_positions_B = []
        new_names_B = []
        new_resids_B = []
        new_resnames_B = []

        residuesB = uniA.select_atoms("all").residues
        for i, residue in enumerate(residuesB):
            for atom in residue:
                new_names_B.append(atom.name)
                new_positions_B.append(atom.positions)
            new_resnames_B.append(residue.resname)
            new_resids_B.append(i+1)

        # add chain Label A
        new_atoms_u_A = mda.Universe.empty(n_atoms=len(new_positions_A),
                                           n_residues=len(new_positions_A),
                                           atom_resindex=np.arange(len(new_positions_A)),
                                           trajectory=True)
        new_atoms_u_A.add_TopologyAttr("name", new_names_A)
        new_atoms_u_A.add_TopologyAttr("resid", new_resids_A)
        new_atoms_u_A.add_TopologyAttr("resname", new_resnames_A)
        new_atoms_u_A.add_TopologyAttr("chainID")
        new_atoms_u_A.atoms.chainIDs = "A"
        new_atoms_u_A.atoms.positions = np.array(new_positions_A)

        # add chain Label B
        new_atoms_u_B = mda.Universe.empty(n_atoms=len(new_positions_B),
                                           n_residues=len(new_positions_B),
                                           atom_resindex=np.arange(len(new_positions_B)),
                                           trajectory=True)
        new_atoms_u_B.add_TopologyAttr("name", new_names_B)
        new_atoms_u_B.add_TopologyAttr("resid", new_resids_B)
        new_atoms_u_B.add_TopologyAttr("resname", new_resnames_B)
        new_atoms_u_B.add_TopologyAttr("chainID")
        new_atoms_u_B.atoms.chainIDs = "B"
        new_atoms_u_B.atoms.positions = np.array(new_positions_B)

        # combine and save to a file
        merged = mda.Merge(new_atoms_u_A.atoms, new_atoms_u_B.atoms)
        pdb_writer.write_fragments_pdb("../temp/initial_dimer.pdb", merged)

        # add the cryst1_line to the top of the file
        cryst1_line = "CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1           1\n"
        with open("../temp/initial_dimer.pdb", "r") as pdb_file:
            pdb_content = pdb_file.readlines()
        pdb_content.insert(0, cryst1_line)
        with open("../temp/initial_dimer.pdb", "w") as modified_file:
            modified_file.writelines(pdb_content)

        # load into OpenMM and run an energy minimisation
        self.__minimise_structure("../temp/initial_dimer.pdb")

        # when the program starts check whether ./temp exists. If it doesn't, make it


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
        with open(self.output_file, 'w') as f:
            PDBFile.writeFile(simulation.topology, minimized_positions, f)
        #simulation.reporters.append(PDBReporter('output.pdb', 1000))
        #simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))
        #simulation.step(10000)

    def __align_principal_axis_to_z(selection: MDAnalysis.AtomGroup):
        # Calculate the principal axes
        principal_axes = selection.principal_axes()
        # The first principal axis corresponds to the largest eigenvalue
        main_axis = principal_axes[0]
        # Calculate the rotation matrix to align this axis with the z-axis
        R = rotation_matrix(main_axis, [0, 0, 1])
        # Apply the rotation
        selection.rotate(R)

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
        average_radius = np.average(radii)

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