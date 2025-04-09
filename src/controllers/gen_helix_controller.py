from src.controllers.controller import Controller
from PeptideBuilder import Geometry
import PeptideBuilder
import Bio.PDB
import MDAnalysis as mda
import numpy as np
from pymol import cmd


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

        #geo_A = Geometry.geometry(self.sequence_A)
        #geo_A.phi = -60
        #geo_A.psi_im1 = -40
        #structure_A = PeptideBuilder.initialize_res(geo_A)
        #PeptideBuilder.add_residue(structure_A, geo_A)
        #PeptideBuilder.add_terminal_OXT(structure_A)
        if self.sequence_B:
            geo_B = Geometry.geometry(self.sequence_B)
            geo_B.phi = -60
            geo_B.psi_im1 = -40
            structure_B = PeptideBuilder.initialize_res(geo_B)
            PeptideBuilder.add_residue(structure_B, geo_B)
            PeptideBuilder.add_terminal_OXT(structure_B)

        # process as a homodimer
        if not self.sequence_B:
            out = Bio.PDB.PDBIO()
            out.set_structure(structure_A)
            out.save(self.output_file)
        else:
            # process as a heterodimer
            pass #both A and B

        # add hydrogens using open source Pymol
        cmd.reinitialize()
        cmd.load(self.output_file, "protein")
        cmd.h_add("protein")
        cmd.save(self.output_file, "protein")
        cmd.quit()



        self.__calculate_average_helix_radius(self.output_file)


    def __calculate_average_helix_radius(self, pdb_file) -> float:
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

        return max_radius







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