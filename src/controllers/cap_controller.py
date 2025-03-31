import MDAnalysis as mda
import numpy as np

from src.controllers.controller import Controller
from src.models import pdb_reader, pdb_writer, pdb_cleaner


class CapController(Controller):

    def __init__(self, input_file, input_directory):
        super().__init__()
        self.input_file = input_file
        self.input_directory = input_directory

        self.universe_list = []

    def validate_inputs(self):
        if self.input_file:
            universe = pdb_reader.read_pdb_file(self.input_file)
            self.universe_list.append(universe)
        else:
            self.universe_list = pdb_reader.read_pdb_folder(self.input_directory)

    def run_controller(self):
        for universe in self.universe_list:

            new_positions_A = []
            new_names_A = []
            new_resids_A = []
            new_resnames_A = []

            new_positions_B = []
            new_names_B = []
            new_resids_B = []
            new_resnames_B = []

            chain_A_residues = universe.select_atoms("chainid A").residues
            for i, residue in enumerate(chain_A_residues):
                resnum = residue.resid

                # if the residue had an immediate left neighbour (resid - 1) then we don't need to cap it.
                has_prev = i > 0 and chain_A_residues[i-1].resid == resnum - 1
                # if the residue had an immediate right neighbour (resid +1) then we don't need to cap it.
                has_next = i < len(chain_A_residues) - 1 and chain_A_residues[i + 1].resid == resnum + 1

                if not has_prev:
                    print("Process hydrogen cap on the N atom")
                    N = residue.atoms.select_atoms("name N")[0].position
                    C_alpha = residue.atoms.select_atoms("name CA")[0].position
                    bond_vector = C_alpha - N
                    bond_vector /= np.linalg.norm(bond_vector)
                    NH_bond_length = 1.0
                    H_position = N - NH_bond_length * bond_vector

                    new_positions_A.append(H_position)
                    new_names_A.append("HNT")
                    new_resids_A.append(resnum)
                    new_resnames_A.append(residue.resname)

                if not has_next:
                    print("Process hydrogen cap on the C atom")
                    C = residue.atoms.select_atoms("name C")[0].position
                    C_alpha = residue.atoms.select_atoms("name CA")[0].position
                    bond_vector = C_alpha - C
                    bond_vector /= np.linalg.norm(bond_vector)
                    CH_bond_length = 1.0
                    H_position = C - CH_bond_length * bond_vector

                    new_positions_A.append(H_position)
                    new_names_A.append("HCT")
                    new_resids_A.append(resnum)
                    new_resnames_A.append(residue.resname)

            chain_B_residues = universe.select_atoms("chainid B").residues
            for i, residue in enumerate(chain_B_residues):
                resnum = residue.resid

                # if the residue had an immediate left neighbour (resid - 1) then we don't need to cap it.
                has_prev = i > 0 and chain_B_residues[i - 1].resid == resnum - 1
                # if the residue had an immediate right neighbour (resid +1) then we don't need to cap it.
                has_next = i < len(chain_B_residues) - 1 and chain_B_residues[i + 1].resid == resnum + 1

                if not has_prev:
                    print("Process hydrogen cap on the N atom")
                    N = residue.atoms.select_atoms("name N")[0].position
                    C_alpha = residue.atoms.select_atoms("name CA")[0].position
                    bond_vector = C_alpha - N
                    bond_vector /= np.linalg.norm(bond_vector)
                    NH_bond_length = 1.0
                    H_position = N - NH_bond_length * bond_vector

                    new_positions_B.append(H_position)
                    new_names_B.append("HNT")
                    new_resids_B.append(resnum)
                    new_resnames_B.append(residue.resname)

                if not has_next:
                    print("Process hydrogen cap on the C atom")
                    C = residue.atoms.select_atoms("name C")[0].position
                    C_alpha = residue.atoms.select_atoms("name CA")[0].position
                    bond_vector = C_alpha - C
                    bond_vector /= np.linalg.norm(bond_vector)
                    CH_bond_length = 1.0
                    H_position = C - CH_bond_length * bond_vector

                    new_positions_B.append(H_position)
                    new_names_B.append("HCT")
                    new_resids_B.append(resnum)
                    new_resnames_B.append(residue.resname)

            if new_positions_A:
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

            if new_positions_B:
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

            if new_positions_A and new_positions_B:
                # Merge old universe with new hydrogen atoms
                merged = mda.Merge(universe.atoms, new_atoms_u_A.atoms, new_atoms_u_B.atoms)
                pdb_writer.write_pdb_to_xyz(merged, "output_temp.xyz")
                pdb_writer.write_fragments_pdb("output_temp.pdb", merged)
                #u = mda.Universe("output_test.pdb")
                #merged = pdb_cleaner.sort_universe(u)
                #pdb_writer.write_fragments_pdb("output_test.pdb", merged)
            elif new_positions_A and not new_positions_B:
                merged = mda.Merge(universe.atoms, new_atoms_u_A.atoms)
                #merged.atoms.write("output_test.pdb")
                pdb_writer.write_fragments_pdb("output_test.pdb", merged)
            elif not new_positions_A and new_positions_B:
                merged = mda.Merge(universe.atoms, new_atoms_u_B.atoms)
                #merged.atoms.write("output_test.pdb")
                pdb_writer.write_fragments_pdb("output_test.pdb", merged)
            else:
                #universe.atoms.write("output_test.pdb")
                pdb_writer.write_fragments_pdb("output_test.pdb", universe)
