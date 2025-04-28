from sys import stdout
from abc import ABC, abstractmethod
from typing import Dict, Union
import numpy as np
from openmm.app import *
from openmm import *
import openmm.unit as unit
from openmm import CustomExternalForce


class Controller(ABC):

    def __init__(self):
        pass

    @abstractmethod
    def run_controller(self):
        pass

    @abstractmethod
    def validate_controller(self, config_section: Dict[str, Union[str, int, float, bool]]) -> bool:
        pass

    def _minimise_structure(self, input_file: str, output_file: str, update_hydrogens=False):
        #get the unit cell dimensions
        with open(input_file, "r") as file:
            for line in file:
                if line.startswith("CRYST1"):
                    # Extract the six cell dimensions from the CRYST1 line
                    cell_params = line[6:54].strip().split()
                    cell = np.array([float(val) for val in cell_params])
                    break
            else:
                raise ValueError("No CRYST1 line found in the PDB file")

        z_length = float(cell[2])/10 # converted to nanometers

        z0 = (z_length/2) * unit.nanometer  # center of bilayer
        core_half_thickness = 1.5 * unit.nanometer  # thickness of hydrophobic core
        k_base_penalty = 5.0  # base penalty (kJ/mol/nm^2)

        # Define hydrophobicity scores (simple: hydrophobic=0, hydrophilic=1)
        hydrophobic_residues = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'TYR'}
        hydrophilic_residues = {'ASP', 'GLU', 'ASN', 'GLN', 'ARG', 'LYS', 'HIS', 'SER', 'THR'}

        # Default hydrophobicity map: residue name → penalty multiplier
        residue_penalty = {}

        for resname in hydrophobic_residues:
            residue_penalty[resname] = 0.2  # low penalty

        for resname in hydrophilic_residues:
            residue_penalty[resname] = 1.0  # full penalty

        # Other residues (like GLY, PRO, CYS) → medium penalty
        default_penalty = 0.5

        # Create a force
        # "Only penalize the atom if it leaves the membrane core, and make the penalty quadratic in how far outside it is, scaled by its hydrophobicity."
        # Each part explained:
        # z: the z-position of the particle.
        # z0: the center of the bilayer (usually z0 = 0 nm).
        # core_half_thickness: half the thickness of the membrane hydrophobic core (say, ±1.5 nm).
        # abs(z - z0): the distance of the atom from the membrane center in the z direction.
        # abs(z-z0) - core_half_thickness:
        # If this is positive, it means the atom is outside the membrane core.
        # If this is negative, it means the atom is inside the membrane core.
        # step(abs(z-z0) - core_half_thickness):
        # This is 0 if inside the core (safe inside).
        # This is 1 if outside the core (apply penalty).
        # step(x) is a function that is 0 when x < 0, and 1 when x > 0.
        # (abs(z-z0) - core_half_thickness)^2:
        # Once you're outside, the penalty grows quadratically with distance from the membrane boundary.
        # Small violations give small penalties; large violations are punished strongly.
        # hydro_penalty:
        # A strength parameter that depends on the hydrophobicity of the residue:
        # Small for hydrophobic residues → they can survive near the edge without huge penalty.
        # Large for hydrophilic residues → they are punished heavily for being outside.

        bilayer_force = CustomExternalForce(
            "hydro_penalty * step(abs(z-z0)-core_half_thickness) * (abs(z-z0) - core_half_thickness)^2"
        )

        # Global parameters
        bilayer_force.addGlobalParameter("z0", z0.value_in_unit(unit.nanometer))
        bilayer_force.addGlobalParameter("core_half_thickness", core_half_thickness.value_in_unit(unit.nanometer))

        # Per-particle parameters
        bilayer_force.addPerParticleParameter("hydro_penalty")  # penalty strength for this atom

        pdb = PDBFile(input_file)

        # Loop over atoms and assign hydrophobicity
        for atom in pdb.topology.atoms():
            resname = atom.residue.name
            penalty_multiplier = residue_penalty.get(resname, default_penalty)
            atom_penalty = k_base_penalty * penalty_multiplier  # scale by base penalty
            bilayer_force.addParticle(atom.index, [atom_penalty])

        forcefield = ForceField('amber99sb.xml', 'tip3p.xml')

        if update_hydrogens:
            modeller = Modeller(pdb.topology, pdb.positions)
            modeller.addHydrogens(forcefield)  # ← ADD correct hydrogens
            with open("../temp/dimer_h.pdb", "w") as file:
                PDBFile.writeFile(modeller.topology, modeller.positions, file)
            pdb = PDBFile('../temp/dimer_h.pdb')

        system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1 * unit.nanometer,
                                         constraints=HBonds)

        # Add the force to the system
        system.addForce(bilayer_force)

        integrator = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
        simulation = Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy()

        minimized_positions = simulation.context.getState(getPositions=True).getPositions()
        simulation.reporters.append(PDBReporter('output.pdb', 1000))
        simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))
        simulation.step(1)
        with open(output_file, 'w') as f:
            PDBFile.writeFile(simulation.topology, minimized_positions, f)

