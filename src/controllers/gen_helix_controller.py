from src.controllers.controller import Controller
from PeptideBuilder import Geometry
import PeptideBuilder
import Bio.PDB

class GenHelixController(Controller):

    def __init__(self, sequence_A, sequence_B, output_file):
        super().__init__()

        self.sequence_A = sequence_A
        self.sequence_B = sequence_B
        self.output_file = output_file

    def run_controller(self):

        geo_A = Geometry.geometry(self.sequence_A)
        geo_A.phi = -60
        geo_A.psi_im1 = -40
        structure_A = PeptideBuilder.initialize_res(geo_A)
        PeptideBuilder.add_residue(structure_A, geo_A)
        PeptideBuilder.add_terminal_OXT(structure_A)
        if self.sequence_B:
            geo_B = Geometry.geometry(self.sequence_B)
            geo_B.phi = -60
            geo_B.psi_im1 = -40
            structure_B = PeptideBuilder.initialize_res(geo_B)
            PeptideBuilder.add_residue(structure_B, geo_B)
            PeptideBuilder.add_terminal_OXT(structure_B)

        if not self.sequence_B:
            out = Bio.PDB.PDBIO()
            out.set_structure(structure_A)
            out.save(self.output_file)
        else:
            pass #both A and B





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