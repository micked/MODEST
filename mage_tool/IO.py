
class Mutation:
    def __init__(self, mut_format, mut, pos=0, mut_type="", ref_genome=""):
        """Mutation
        
        mut_format can be:
        arrow: position required
            point mutation: A->T, pos, mut_change="Point_Mutation"
            insertion:      A, pos, mut_change="Insertion"
            deletion        3, pos, mut_change="Deletion"
        eq: position required
            point mutation: A=T, pos
            insertion:      =T or A=AT, pos
            deletion:       A=, pos
        genome: search genome for mutation (eq format)
            AATGATA[ATG=GT]ATGATA
            
        """
        
        if mut_format.lower() == "arrow":
            self._parse_arrow(mut, pos, mut_type)
        elif mut_format.lower() == "eq":
            self._parse_eq(mut, pos)
        else:
            raise Exception("Format: \"{}\" unknown".format(mut_format))
    
    #
    # Printing
    #
    
    def __str__(self):
        return "Mutation: [{}={}] at pos {}".format(
            self.before,
            self.after,
            self.position)
    
    #
    # Parsers
    #
    
    def _parse_arrow(self, mut, pos, mut_type):
        """Parse arrow format"""
        self.position = pos
        if mut_type.lower() == "point_mutation":  
            self.before = mut.split("->")[0]
            self.after = mut.split("->")[1]
        elif mut_type.lower() == "insertion":
            self.before = ""
            self.after = mut
        elif mut_type == "deletion" or "large_deletion":
            self.before = ref_genome[mut_pos:mut_pos+int(mut)]
            self.after = ""
        else:
            raise Exception("Unknown mut_type: " + mut_type)
            
    def _parse_eq(self, mut, pos):
        self.position = pos
        mut = mut.strip("[]")
        self.before, self.after = mut.split("=")
