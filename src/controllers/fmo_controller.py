import os
from typing import List, Any, Tuple

import MDAnalysis as mda

from src.controllers.controller import Controller
from src.logger_config import get_logger
from src.models.fmo_object import FMOObject


class FMOController(Controller):

    def __init__(self, basis: str, theory: str):
        super().__init__()

        self.logger = get_logger(__name__)
        self.output_location = ""
        self.pdb_files = []
        self.universe_list = []

        self.basis = basis
        self.theory = theory

        self.fmo_objects = []

        #self.instruction_list = []
        #self.fragment_name_list = []
        #self.system_charge_list = []
        #self.indat_list = []
        #self.icharge_list = []
        #self.num_fragments_list = []
        #self.mult_list = []
        #self.frag_boundary_list = []
        #self.fmoxyz_list = []

    def validate_inputs(self, pdb_file, pdb_folder, output_location):
        if pdb_file:
            self.pdb_files.append(pdb_file)
            self.universe_list.append(mda.Universe(pdb_file))
        else:
            self.pdb_files = [os.path.join(pdb_folder, f) for f in os.listdir(pdb_folder) if f.endswith(".pdb")]
            for file in self.pdb_files:
                self.universe_list.append(mda.Universe(file))

        self.output_location = output_location

    def run_controller(self):
        # operate at scale over all universe objects
        # an iteration results in an FMO GAMESS input file
        for universe in self.universe_list:
            fmo_object = FMOObject()

            # define the GAMESS output file location and file name
            file_path = universe.filename
            _, base_name = os.path.split(file_path)
            name_without_ext, _ = os.path.splitext(base_name)
            fmo_object.output_file = os.path.join(self.output_location, f"{name_without_ext}.inp")

            chain_A_residues = universe.select_atoms("chainid A").residues
            chain_B_residues = universe.select_atoms("chainid B").residues

            chain_A_uni = self.__reorder_atoms_in_chain(chain_A_residues)
            chain_B_uni = self.__reorder_atoms_in_chain(chain_B_residues)

            fmo_object.indat = self.__build_indat(universe)
            fmo_object.fragment_names = self.__build_fragnames(universe)

            icharge, system_charge = self.__build_icharge(universe)
            fmo_object.icharge = icharge
            fmo_object.system_charge = system_charge

            fmo_object.instruction = self.__add_gamess_instructions(system_charge)
            fmo_object.num_fragments = self.__build_nfrag(universe)
            fmo_object.mult = self.__build_multiplicity(universe)
            fmo_object.fmohyb = self.__write_hybrid_orbitals()

            chain_A_frag_boundary_ids = self.__find_fragment_boundaries(chain_A_uni)
            chain_B_frag_boundary_ids = self.__find_fragment_boundaries(chain_B_uni)
            fmo_object.frag_boundary = "".join([chain_A_frag_boundary_ids, chain_B_frag_boundary_ids])

            chain_A_fmoxyz_str, chain_A_element_list = self.__build_FMOXYZ(chain_A_uni)
            chain_B_fmoxyz_str, chain_B_element_list = self.__build_FMOXYZ(chain_B_uni)
            fmo_object.atomic_numbers = "".join(set(chain_A_element_list + chain_B_element_list))
            fmo_object.fmoxyz = "".join([chain_A_fmoxyz_str, chain_B_fmoxyz_str])

            self.fmo_objects.append(fmo_object)

        self.__write_fmo_files()


    def __write_fmo_files(self):
        for fmo_object in self.fmo_objects:
            output_file = fmo_object.output_file
            output_str = fmo_object.write_instructions()

            with open(output_file, "w") as file:
                file.write(output_str)

    def __build_multiplicity(self, u: mda.Universe) -> str:
        mult = "mult(1)="
        for residue in u.residues:
            mult = "".join([mult, "1", ","])
        mult = "".join([mult, "\n"])

        return mult


    def __build_nfrag(self, u: mda.Universe) -> str:
        count = 0
        for residue in u.residues:
            count = count + 1
        return f"nfrag={str(count)}\n"


    def __build_icharge(self, u: mda.Universe) -> Tuple[str, int]:
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
        system_charge = 0
        icharge = "icharg(1)="
        for residue in u.residues:
            charge = residue_charge[residue.resname]
            icharge = "".join([icharge, str(charge), ","])
            system_charge = system_charge + charge

        icharge = "".join([icharge, "\n"])

        return icharge, system_charge

    def __build_fragnames(self, u: mda.Universe) -> str:
        fragnames = "frgnam(1)="
        for residue in u.residues:
            fragnames = ",".join([fragnames, residue.resname])
        fragnames = "".join([fragnames, "\n"])

        return fragnames

    def __build_indat(self, u: mda.Universe) -> str:
        indat_str = "indat(1)="
        count = 1
        for residue in u.residues:
            for atom in residue.atoms:
                indat_str = "".join([indat_str, str(count)])
                indat_str = "".join([indat_str, ","])
            indat_str = "".join([indat_str, "\n"])
            count = count + 1

        return indat_str

    def __build_FMOXYZ(self, u: mda.Universe) -> Tuple[str, List[str]]:
        element_atomic_number = {
            "H": 1,  # Hydrogen
            "He": 2,  # Helium
            "Li": 3,  # Lithium
            "Be": 4,  # Beryllium
            "B": 5,  # Boron
            "C": 6,  # Carbon
            "N": 7,  # Nitrogen
            "O": 8,  # Oxygen
            "F": 9,  # Fluorine
            "Ne": 10,  # Neon
            "Na": 11,  # Sodium
            "Mg": 12,  # Magnesium
            "Al": 13,  # Aluminum
            "Si": 14,  # Silicon
            "P": 15,  # Phosphorus
            "S": 16,  # Sulfur
            "Cl": 17  # Chlorine
        }

        xyz_list = []
        element_list = []

        # Loop through all atoms in the Universe
        for atom in u.atoms:
            element = atom.name[:2] if atom.name[:2] in element_atomic_number else atom.name[0]  # Extract element symbol from atom name
            atomic_number = element_atomic_number[element]  # Placeholder (MDAnalysis doesn't store atomic numbers)
            element_data = "".join([element, " ", str(atomic_number), "\n"])
            element_list.append(element_data)
            xyz_entry = f"{element} {atomic_number} {atom.position[0]:.3f} {atom.position[1]:.3f} {atom.position[2]:.3f} \n"
            xyz_list.append(xyz_entry)

        return "".join(xyz_list), element_list

    def __find_fragment_boundaries(self, u: mda.Universe) -> str: #List[tuple[Any, Any]]:
        residues = u.residues

        # List to store recorded atom indices
        #c_ca_indices = []
        boundary_str = ""

        # Loop through residues except the last one
        for i in range(len(residues) - 1):
            current_res = residues[i]
            next_res = residues[i + 1]

            # Check if the next residue follows sequentially
            if current_res.resid + 1 == next_res.resid:
                # Get "C" and "CA" atom indices from the current residue
                c_atom = current_res.atoms.select_atoms("name C")
                ca_atom = current_res.atoms.select_atoms("name CA")

                # Record atom indices if found
                if len(c_atom) > 0 and len(ca_atom) > 0:
                    boundary_str = "".join([boundary_str, "-", str(ca_atom[0].index), " ", str(c_atom[0].index), " ", self.basis, "\n"])
                    #c_ca_indices.append(( "".join(["-",str(ca_atom[0].index)]), str(c_atom[0].index)))

        # Print results
        #print("Recorded C and CA atom indices:", c_ca_indices)
        return boundary_str


    def __reorder_atoms_in_chain(self, residues: mda.AtomGroup) -> mda.Universe:

        # List to store reordered AtomGroups
        reordered_atom_groups = []

        # Iterate through each residue in Chain A
        for res in residues:
            atoms = res.atoms  # Get all atoms in the residue

            # Select CA and C atoms
            ca_atom = atoms.select_atoms("name CA")
            c_atom = atoms.select_atoms("name C")

            # Select all other atoms except CA and C
            other_atoms = atoms.select_atoms("not name CA C")

            # Merge atoms in new order (other atoms first, then CA and C)
            reordered_atom_group = other_atoms + ca_atom + c_atom

            # Store the reordered AtomGroup
            reordered_atom_groups.append(reordered_atom_group)

        # Merge all reordered AtomGroups into a new universe
        new_u = mda.Merge(*reordered_atom_groups)  # Unpacking AtomGroups

        return new_u

    def __add_gamess_instructions(self, system_charge: int) -> str:
        #instructions = []
        # The basis, theory and memory instructions are consistent across the PDB files in the one run
        if self.basis == "6-31G*":
            basis_str = " $BASIS GBASIS=N31 NGAUSS=6 NDFUNC=1 $END\n"
        else:
            basis_str = " $BASIS GBASIS=STO NGAUSS=3 $END\n"
        memory_str = " $SYSTEM MEMORY=1000000 $END\n"

        # The control line differs according to the system charge
        #for system_charge in self.system_charge_list:
        if self.theory == "HF":
            control_str = f" $CONTROL SCFTYP=RHF RUNTYP=ENERGY ICHARG={str(system_charge)} $END\n"
        else:
            control_str = f" $CONTROL SCFTYP=RHF RUNTYP=ENERGY MPLEVL=2 ICHARG={str(system_charge)} $END\n"

        return "".join([basis_str, control_str, memory_str])
            #instructions.append("".join([basis_str, control_str, memory_str]))


    def __write_hybrid_orbitals(self) -> str:
        basis_631G_star = []
        basis_631G_star.append("6-31G* 15 5\n")
        basis_631G_star.append("0 1   0.995621   0.028507   0.000000   0.000000   0.000000\n")
        basis_631G_star.append("     -0.013353   0.000000   0.000000  -0.000001  -0.002076\n")
        basis_631G_star.append("     -0.002076  -0.002076   0.000000   0.000000   0.000000\n")
        basis_631G_star.append("0 1   0.097077  -0.179454   0.178051   0.308393   0.125888\n")
        basis_631G_star.append("     -0.181305   0.093381   0.161741   0.066029   0.000102\n")
        basis_631G_star.append("     -0.020671   0.005300  -0.020773  -0.008476  -0.014681\n")
        basis_631G_star.append("0 1   0.097077  -0.179454  -0.356101  -0.000000   0.125888\n")
        basis_631G_star.append("     -0.181304  -0.186762  -0.000000   0.066030  -0.031058\n")
        basis_631G_star.append("0 1   0.097077  -0.179454   0.178051  -0.308393   0.125888\n")
        basis_631G_star.append("     -0.181304   0.093381  -0.161741   0.066029   0.000102\n")
        basis_631G_star.append("     -0.020672   0.005300   0.020773  -0.008476   0.014681\n")
        basis_631G_star.append("1 0  -0.097061   0.179413   0.000000   0.000000   0.377708\n")
        basis_631G_star.append("      0.181288   0.000000   0.000000   0.198114  -0.010489\n")
        basis_631G_star.append("     -0.010489   0.036248  -0.000000  -0.000000  -0.000000\n")

        basis_sto_3g = []
        basis_sto_3g.append("STO-3G 5 5\n")
        basis_sto_3g.append("0 1   0.991926   0.038355   0.000000   0.000000  -0.000000\n")
        basis_sto_3g.append("0 1  -0.110715   0.313793  -0.233504  -0.404441  -0.165113\n")
        basis_sto_3g.append("1 0  -0.110716   0.313795  -0.000000  -0.000000   0.495336\n")
        basis_sto_3g.append("0 1  -0.110715   0.313793   0.467008   0.000000  -0.165112\n")
        basis_sto_3g.append("0 1  -0.110715   0.313793  -0.233504   0.404441  -0.165112\n")

        if self.basis == "STO-3G":
            return "".join(basis_sto_3g)
        else:
            return "".join(basis_631G_star)


