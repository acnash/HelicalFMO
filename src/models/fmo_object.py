

class FMOObject:

    def __init__(self):
        self.output_file = ""
        self.instruction = ""
        self.fragment_names = ""
        self.system_charge = 0
        self.indat = ""
        self.icharge = ""
        self.num_fragments = ""
        self.mult = ""
        self.frag_boundary = ""
        self.fmoxyz = ""
        self.title = ""
        self.atomic_numbers = ""
        self.fmohyb = ""

    def write_instructions(self) -> str:
        output_str = self.instruction
        output_str = "".join([output_str, " $FMO\n"])
        output_str = "".join([output_str, "  ", self.num_fragments])
        output_str = "".join([output_str, self.fragment_names])
        output_str = "".join([output_str, self.icharge])
        output_str = "".join([output_str, self.mult])
        output_str = "".join([output_str, self.indat])
        output_str = "".join([output_str, " $END\n"])
        output_str = "".join([output_str, " $fmoprp nprint=9 naodir=220 modorb=1 ipieda=1 $end\n"])
        output_str = "".join([output_str, " $DATA\n"])
        output_str = "".join([output_str, self.title])
        output_str = "".join([output_str, "c1\n"])
        output_str = "".join([output_str, self.atomic_numbers])
        output_str = "".join([output_str, " $END\n"])
        output_str = "".join([output_str, " $FMOHYB\n"])
        output_str = "".join([output_str, self.fmohyb])
        output_str = "".join([output_str, " $END\n"])
        output_str = "".join([output_str, " $FMOBND\n"])
        output_str = "".join([output_str, self.frag_boundary])
        output_str = "".join([output_str, " $END\n"])
        output_str = "".join([output_str, " $FMOXYZ\n"])
        output_str = "".join([output_str, self.fmoxyz])
        output_str = "".join([output_str, " $END\n"])

        return output_str
