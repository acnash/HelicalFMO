from typing import Dict, Union

import MDAnalysis as mda
import numpy as np

from src.controllers.controller import Controller
from src.models import pdb_reader, pdb_writer
from src.logger_config import get_logger


class CapController(Controller):

    def __init__(self):
        super().__init__()
        logger = get_logger(__name__)

        self.input_file = None
        self.input_directory = None

        self.universe_list = []

    def validate_controller(self, config_section: Dict[str, Union[str, int, float, bool]]):
        self.input_file = config_section.get("file")
        self.input_directory = config_section.get("folder")

        if self.input_file:
            universe = pdb_reader.read_pdb_file(self.input_file)
            self.universe_list.append(universe)
        else:
            self.universe_list = pdb_reader.read_pdb_folder(self.input_directory)


    #def validate_inputs(self):
    #    if self.input_file:
    #        universe = pdb_reader.read_pdb_file(self.input_file)
    #        self.universe_list.append(universe)
    #    else:
    #        self.universe_list = pdb_reader.read_pdb_folder(self.input_directory)

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
                total_charge = self.__calculate_charge(merged)

                # make a reduced file for PSI4
                self.__reduce_structure(merged, "just_backbone.pdb")  # via MDAnalysis
                just_backbone_universe = mda.Universe("just_backbone.pdb")  # via plain text editing

                pdb_writer.write_pdb_to_psi4in(just_backbone_universe, "output_temp.psi4in", 0)

                # run PSI4
                # get the output
                # remap positions of hydrogens to the 'merged' universe object

                #pdb_writer.write_fragments_pdb("output_temp.pdb", just_backbone_universe)

                # pdb_writer.write_pdb_to_psi4in(merged, "output_temp.psi4in", total_charge)
                # pdb_writer.write_fragments_pdb("output_temp.pdb", merged)

            elif new_positions_A and not new_positions_B:
                merged = mda.Merge(universe.atoms, new_atoms_u_A.atoms)
                total_charge = self.__calculate_charge(merged)
                pdb_writer.write_pdb_to_psi4in(merged, "output_temp.psi4in", total_charge)
                pdb_writer.write_fragments_pdb("output_test.pdb", merged)
            elif not new_positions_A and new_positions_B:
                merged = mda.Merge(universe.atoms, new_atoms_u_B.atoms)
                total_charge = self.__calculate_charge(merged)
                pdb_writer.write_pdb_to_psi4in(merged, "output_temp.psi4in", total_charge)
                pdb_writer.write_fragments_pdb("output_test.pdb", merged)
            else:
                total_charge = self.__calculate_charge(universe)
                pdb_writer.write_pdb_to_psi4in(universe, "output_temp.psi4in", total_charge)
                pdb_writer.write_fragments_pdb("output_test.pdb", universe)

    def __calculate_charge(self, u: mda.Universe) -> int:
        residue_charge = {
            "GLY": 0,
            "VAL": 0,
            "ALA": 0,
            "LEU": 0,
            "ILE": 0,
            "TRP": 0,
            "TYR": 0,
            "PHE": 0,
            "MET": 0,
            "SER": 0,
            "HIS": 0,
            "HID": 0,
            "HIE": 0,
            "HIP": 1,
            "GLN": 0,
            "GLU": -1,
            "ARG": 1,
            "LYS": 1,
            "PRO": 0,
            "THR": 0,
            "CYS": 0,
            "ASP": -1,
            "ASN": 0,
            "HYP": 0
        }

        total_charge = 0
        for residue in u.residues:
            total_charge = total_charge + residue_charge[residue.resname]

        return total_charge

    def __reduce_structure(self, u: mda.Universe, output_file: str) -> None:

        writer = mda.Writer(output_file, multiframe=False)
        # Loop through each residue
        for res in u.residues:
            if res.resname in {"GLY", "PRO"}:
                # Keep all atoms for glycine and proline
                selected_atoms = res.atoms
            else:
                # Select backbone atoms and CB
                selected_atoms = res.atoms.select_atoms("name N CA C O CB H HA HCT HNT")

            if selected_atoms:
                writer.write(selected_atoms)

        writer.close()

        with open(output_file, "r") as file:
            content = file.read()

        content = content.replace("CB ", "HX ")

        with open(output_file, "w") as file:
            file.write(content)

