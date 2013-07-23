amino_acids = [ "ARG",
                "HIE",
                "HID",
                "HIP",
                "LYS",
                "ASP",
                "GLU",
                "SER",
                "THR",
                "ASN",
                "GLN",
                "CYS",
                "GLY",
                "PRO",
                "ALA",
                "ILE",
                "LEU",
                "MET",
                "PHE",
                "TRP",
                "TYR",
                "VAL" ]

def_parameters = {}
# C, alpha-C bonds
for amino_acid in amino_acids:
    def_parameters[(amino_acid, ("C", "CA"))] = (0.1, 0.5)
# alpha-C, beta-C bonds
def_parameters[("ARG", ("CA", "CB"))] = (0.2, 0.5)
def_parameters[("HIE", ("CA", "CB"))] = (0.3, 0.5)
def_parameters[("HID", ("CA", "CB"))] = (0.3, 0.5)
def_parameters[("HIP", ("CA", "CB"))] = (0.3, 0.5)
def_parameters[("LYS", ("CA", "CB"))] = (0.2, 0.5)
def_parameters[("ASP", ("CA", "CB"))] = (0.5, 0.5)
def_parameters[("GLU", ("CA", "CB"))] = (0.3, 0.5)
def_parameters[("SER", ("CA", "CB"))] = (1.0, 0.5)
def_parameters[("THR", ("CA", "CB"))] = (1.0, 0.5)
def_parameters[("ASN", ("CA", "CB"))] = (0.5, 0.5)
def_parameters[("GLN", ("CA", "CB"))] = (0.3, 0.5)
def_parameters[("CYS", ("CA", "CB"))] = (1.0, 0.5)
def_parameters[("ALA", ("CA", "CB"))] = (1.0, 0.5)
def_parameters[("VAL", ("CA", "CB"))] = (1.0, 0.5)
def_parameters[("ILE", ("CA", "CB"))] = (0.5, 0.5)
def_parameters[("LEU", ("CA", "CB"))] = (0.5, 0.5)
def_parameters[("MET", ("CA", "CB"))] = (0.5, 0.5)
def_parameters[("PHE", ("CA", "CB"))] = (0.3, 0.5)
def_parameters[("TYR", ("CA", "CB"))] = (0.3, 0.5)
def_parameters[("TRP", ("CA", "CB"))] = (0.3, 0.5)
# beta-C, gamma-C bonds
def_parameters[("ARG", ("CB", "CG"))] = (0.3, 0.5)
def_parameters[("HIE", ("CB", "CG"))] = (0.5, 0.5)
def_parameters[("HID", ("CB", "CG"))] = (0.5, 0.5)
def_parameters[("HIP", ("CB", "CG"))] = (0.5, 0.5)
def_parameters[("LYS", ("CB", "CG"))] = (0.3, 0.5)
def_parameters[("ASP", ("CB", "CG"))] = (1.0, 0.5)
def_parameters[("GLU", ("CB", "CG"))] = (0.5, 0.5)
def_parameters[("ASN", ("CB", "CG"))] = (1.0, 0.5)
def_parameters[("GLN", ("CB", "CG"))] = (0.5, 0.5)
def_parameters[("ILE", ("CB", "CG1"))] = (1.0, 0.5)
def_parameters[("LEU", ("CB", "CG"))] = (1.0, 0.5)
def_parameters[("MET", ("CB", "CG"))] = (0.7, 0.5)
def_parameters[("PHE", ("CB", "CG"))] = (0.5, 0.5)
def_parameters[("TYR", ("CB", "CG"))] = (0.5, 0.5)
def_parameters[("TRP", ("CB", "CG"))] = (0.4, 0.5)
# gamma-C, delta-C bonds
def_parameters[("ARG", ("CG", "CD"))] = (0.5, 0.5)
def_parameters[("LYS", ("CG", "CD"))] = (0.5, 0.5)
def_parameters[("GLU", ("CG", "CD"))] = (1.0, 0.5)
def_parameters[("GLN", ("CG", "CD"))] = (1.0, 0.5)
def_parameters[("MET", ("CG", "SD"))] = (0.5, 0.5)
# delta, epsilon bonds
def_parameters[("ARG", ("CD", "NE"))] = (0.5, 0.5)
def_parameters[("LYS", ("CD", "CE"))] = (0.5, 0.5)