

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

        split_fragment_names = self.__wrap_on_commas(self.fragment_names)
        output_str = "".join([output_str, split_fragment_names, "\n"])

        split_icharge = self.__wrap_on_commas(self.icharge)
        output_str = "".join([output_str, split_icharge, "\n"])

        split_mult = self.__wrap_on_commas(self.mult)
        output_str = "".join([output_str, split_mult, "\n"])

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

    def __wrap_on_commas(self, s: str, max_length=79) -> str:
        parts = s.split(",")
        lines = []
        current_line = ""

        for i, part in enumerate(parts):
            segment = part + ("," if i < len(parts) - 1 else "")
            if len(current_line) + len(segment) <= max_length:
                current_line += segment
            else:
                lines.append(current_line.strip())
                current_line = segment
        if current_line:
            lines.append(current_line.strip())

        return "\n".join(lines)
