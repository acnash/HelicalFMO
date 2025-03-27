from src.controllers.controller import Controller


class GromacsController(Controller):

    def __init__(self, index_file, topology_file, output_label):
        self.index_file = index_file
        self.topology_file = topology_file
        self.output_label = output_label

    def build_tpr(self, mdp_file, structure_file, sim_type):
        if sim_type == "energy_min":
            output_file = "".join([self.output_label, "em.tpr"])
            command = ["gmx", "grompp", ",-f", mdp_file, "-c", structure_file, "-r", structure_file, "-p", self.topology_file, "-o", output_file]
        elif sim_type == "first_eq":
            output_file = "".join([self.output_label, "eq.tpr"])
            command = ["gmx", "grompp", "-f", mdp_file, "-c", structure_file, "-r", structure_file, "-n", self.index_file, "-p", self.topology_file, "-o", output_file]
        elif sim_type == "next_eq":
            output_file = "".join([self.output_label, "eq.tpr"])
            command = f"gmx grompp -f {mdp_file} -c {structure_file} -r {structure_file} -n {self.index_file} -p {self.topology_file} -o {output_file}"
        else:
            pass

    def run_tpr(self):
        pass
