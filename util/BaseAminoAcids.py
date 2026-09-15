# ------------------------------------------------------------
# Functions related to base amino acid data handling
# ------------------------------------------------------------

# ------------------------------------------------------------
# Base amino acid handlers
# ------------------------------------------------------------
#canonical 3-letter residue names for standard amino acids plus common variants
_BASE_AA_RESNAMES = {
    #standard 20
    "ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE",
    "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
    #variants seen in some pipelines (present in legacy BaseAminoAcids.py)
    "HSD","HID",            #histidine delta-protonated naming in some FFs
    "HSE","HIE",            #histidine epsilon-protonated naming in some FFs
    "HIP",                  #doubly protonated histidine
    "CYX",                  #disulfide-bonded cysteine naming in some FFs
    "CYM"
}

def IsBaseAminoAcid(residue):
    """Return True if `residue` is recognised as a base amino-acid residue name."""
    return str(residue).strip().upper() in _BASE_AA_RESNAMES

# ------------------------------------------------------------
# Extended amino acid handlers
# ------------------------------------------------------------
#canonical 3-letter residue names for standard amino acids plus common variants
_EXTENDED_AA_RESNAMES = {"SEC","MSE"}

def IsExtendedAminoAcid(residue):
    """Return True if `residue` is recognised as an extended amino-acid residue name."""
    return str(residue).strip().upper() in _EXTENDED_AA_RESNAMES

# ------------------------------------------------------------
# Post-translational modification handlers
# ------------------------------------------------------------
#canonical 3-letter residue names for PTMs
_BASE_PTM_RESNAMES = {
    #phosphorylation
    "SEP","TPO","PTR","NEP",
    #acylation
    "ALY","AME","FME",
    #methylation
    "MLY","MLZ","M3L","AGM","MMO","MGN","HIC","MHS","A1I9G","SMC",
    #hydroxylation
    "TRX","TOQ","4HT","4PQ","HRP","HTR","UX8","HYP","BHD","BH2","CSO","TYQ","LYZ","ILX",
    #carboxylation
    "KCX",
    #cyclisation
    "PCA","7YO",
    #other
    "SCH","2CO","CSD","0CS","CSX","CSS","TSY","TQZ"
}

def IsPostTranslationalModification(residue):
    """Return True if `residue` is recognised as a PTM residue name."""
    return str(residue).strip().upper() in _BASE_PTM_RESNAMES

# ------------------------------------------------------------
# External/artificial caps (Schrödinger/Maestro-like, peptide cappers)
# ------------------------------------------------------------
def ExternalCapResnames():
    """Return a set of residue names that are considered external/artificial caps."""
    #keep these uppercase
    return {
        "ACE",        #acetyl cap
        "NME",        #N-methyl cap (common)
        "NMA",        #N-methylacetamide (common)
        "FOR",        #formyl cap
        "FMT",        #formate/formyl variants in some pipelines
        "ACET",       #alternate acetyl naming
        "METH"        #alternate methyl naming
    }

def IsExternalCap(residue):
    """True if residue name matches a known external/artificial cap."""
    return residue.strip().upper() in ExternalCapResnames()

# ------------------------------------------------------------
# Amino acid and derivativatives connectivity functions
# ------------------------------------------------------------
def BaseAminoAcidConnectivity(residue_name, terminal):
    #function returning the connectivity of a given residue
    aa_connectivity = []
    #basic amino acid connectivity data
    if residue_name.upper() == "ALA":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","HB1","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False]]
    elif residue_name.upper() == "ARG":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","CD","sing",False],
        ["CG","HG2","sing",False],
        ["CG","HG3","sing",False],
        ["CD","NE","sing",False],
        ["CD","HD2","sing",False],
        ["CD","HD3","sing",False],
        ["NE","CZ","sing",False],
        ["NE","HE","sing",False],
        ["CZ","NH1","sing",False],
        ["CZ","NH2","doub",False],
        ["NH1","HH11","sing",False],
        ["NH1","HH12","sing",False],
        ["NH2","HH21","sing",False],
        ["NH2","HH22","sing",False]]
    elif residue_name.upper() == "ASN":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","OD1","doub",False],
        ["CG","ND2","sing",False],
        ["ND2","HD21","sing",False],
        ["ND2","HD22","sing",False]]
    elif residue_name.upper() == "ASP":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","OD1","doub",False],
        ["CG","OD2","sing",False],
        ["OD2","HD2","sing",False]]
    elif residue_name.upper() == "CYS":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","SG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["SG","HG","sing",False]]
    elif residue_name.upper() == "GLN":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","CD","sing",False],
        ["CG","HG2","sing",False],
        ["CG","HG3","sing",False],
        ["CD","OE1","doub",False],
        ["CD","NE2","sing",False],
        ["NE2","HE21","sing",False],
        ["NE2","HE22","sing",False]]
    elif residue_name.upper() == "GLU":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","CD","sing",False],
        ["CG","HG2","sing",False],
        ["CG","HG3","sing",False],
        ["CD","OE1","doub",False],
        ["CD","OE2","sing",False],
        ["OE2","HE2","sing",False]]
    elif residue_name.upper() == "GLY":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","HA2","sing",False],
        ["CA","HA3","sing",False],
        ["C","O","doub",False]]
    elif residue_name.upper() == "HIS":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","ND1","sing",True],
        ["CG","CD2","doub",True],
        ["ND1","CE1","doub",True],
        ["ND1","HD1","sing",False],
        ["CD2","NE2","sing",True],
        ["CD2","HD2","sing",False],
        ["CE1","NE2","sing",True],
        ["CE1","HE1","sing",False],
        ["NE2","HE2","sing",False]]
    elif residue_name.upper() == "ILE":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG1","sing",False],
        ["CB","CG2","sing",False],
        ["CB","HB","sing",False],
        ["CG1","CD1","sing",False],
        ["CG1","HG12","sing",False],
        ["CG1","HG13","sing",False],
        ["CG2","HG21","sing",False],
        ["CG2","HG22","sing",False],
        ["CG2","HG23","sing",False],
        ["CD1","HD11","sing",False],
        ["CD1","HD12","sing",False],
        ["CD1","HD13","sing",False],
        #
        ["CG1","CD","sing",False],
        ["CD1","HD11","sing",False],
        ["CD1","HD12","sing",False],
        ["CD1","HD13","sing",False]
        ]
    elif residue_name.upper() == "LEU":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","CD1","sing",False],
        ["CG","CD2","sing",False],
        ["CG","HG","sing",False],
        ["CD1","HD11","sing",False],
        ["CD1","HD12","sing",False],
        ["CD1","HD13","sing",False],
        ["CD2","HD21","sing",False],
        ["CD2","HD22","sing",False],
        ["CD2","HD23","sing",False]]
    elif residue_name.upper() == "LYS":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","CD","sing",False],
        ["CG","HG2","sing",False],
        ["CG","HG3","sing",False],
        ["CD","CE","sing",False],
        ["CD","HD2","sing",False],
        ["CD","HD3","sing",False],
        ["CE","NZ","sing",False],
        ["CE","HE2","sing",False],
        ["CE","HE3","sing",False],
        ["NZ","HZ1","sing",False],
        ["NZ","HZ2","sing",False],
        ["NZ","HZ3","sing",False]]
    elif residue_name.upper() == "MET":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","SD","sing",False],
        ["CG","HG2","sing",False],
        ["CG","HG3","sing",False],
        ["SD","CE","sing",False],
        ["CE","HE1","sing",False],
        ["CE","HE2","sing",False],
        ["CE","HE3","sing",False]]
    elif residue_name.upper() == "PHE":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","CD1","doub",True],
        ["CG","CD2","sing",True],
        ["CD1","CE1","sing",True],
        ["CD1","HD1","sing",False],
        ["CD2","CE2","doub",True],
        ["CD2","HD2","sing",False],
        ["CE1","CZ","doub",True],
        ["CE1","HE1","sing",False],
        ["CE2","CZ","sing",True],
        ["CE2","HE2","sing",False],
        ["CZ","HZ","sing",False]]
    elif residue_name.upper() == "PRO":
        aa_connectivity = [["N","CA","sing",False],
        ["N","CD","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","CD","sing",False],
        ["CG","HG2","sing",False],
        ["CG","HG3","sing",False],
        ["CD","HD2","sing",False],
        ["CD","HD3","sing",False]]
    elif residue_name.upper() == "SER":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","OG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["OG","HG","sing",False]]
    elif residue_name.upper() == "THR":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","OG1","sing",False],
        ["CB","CG2","sing",False],
        ["CB","HB","sing",False],
        ["OG1","HG1","sing",False],
        ["CG2","HG21","sing",False],
        ["CG2","HG22","sing",False],
        ["CG2","HG23","sing",False]]
    elif residue_name.upper() == "TRP":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","CD1","doub",True],
        ["CG","CD2","sing",True],
        ["CD1","NE1","sing",True],
        ["CD1","HD1","sing",False],
        ["CD2","CE2","doub",True],
        ["CD2","CE3","sing",True],
        ["NE1","CE2","sing",True],
        ["NE1","HE1","sing",False],
        ["CE2","CZ2","sing",True],
        ["CE3","CZ3","doub",True],
        ["CE3","HE3","sing",False],
        ["CZ2","CH2","doub",True],
        ["CZ2","HZ2","sing",False],
        ["CZ3","CH2","sing",True],
        ["CZ3","HZ3","sing",False],
        ["CH2","HH2","sing",False]]
    elif residue_name.upper() == "TYR":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","CD1","doub",True],
        ["CG","CD2","sing",True],
        ["CD1","CE1","sing",True],
        ["CD1","HD1","sing",False],
        ["CD2","CE2","doub",True],
        ["CD2","HD2","sing",False],
        ["CE1","CZ","doub",True],
        ["CE1","HE1","sing",False],
        ["CE2","CZ","sing",True],
        ["CE2","HE2","sing",False],
        ["CZ","OH","sing",False],
        ["OH","HH","sing",False]]
    elif residue_name.upper() == "VAL":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG1","sing",False],
        ["CB","CG2","sing",False],
        ["CB","HB","sing",False],
        ["CG1","HG11","sing",False],
        ["CG1","HG12","sing",False],
        ["CG1","HG13","sing",False],
        ["CG2","HG21","sing",False],
        ["CG2","HG22","sing",False],
        ["CG2","HG23","sing",False]]
    #termini
    print("*******terminal*******",terminal,residue_name)
    if "C-term" in terminal:
        aa_connectivity.append(["C","OXT","sing",False])
        aa_connectivity.append(["OXT","HXT","sing",False])
        aa_connectivity.append(["C","OT1","doub",False])
        aa_connectivity.append(["C","OT2","sing",False])
        aa_connectivity.append(["OT2","HT2","sing",False])
    elif "C-cap" in terminal:
        aa_connectivity.append(["C","H","sing",False])
    if "N-term" in terminal:
        aa_connectivity.append(["N","H2","sing",False])
        if residue_name.upper() != "PRO": aa_connectivity.append(["N","H3","sing",False])
    elif "N-cap" in terminal:
        aa_connectivity.append(["N","C*","sing",False])
        aa_connectivity.append(["C*","O*","doub",False])
        aa_connectivity.append(["C*","H","sing",False])
    elif "ipa_termini" in terminal:
        aa_connectivity.append(["N","C#","sing",False])
        aa_connectivity.append(["C#","O#","doub",False])
    return aa_connectivity

def ExtendedAminoAcidConnectivity(residue_name, terminal):
    #function returning the connectivity of a given residue
    aa_connectivity = []
    #extended amino acid connectivity data
    if residue_name.upper() == "SEC":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","SE","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["SE","HG","sing",False]]
    elif residue_name.upper() == "MSE":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","SE","sing",False],
        ["CG","HG2","sing",False],
        ["CG","HG3","sing",False],
        ["SE","CE","sing",False],
        ["CE","HE1","sing",False],
        ["CE","HE2","sing",False],
        ["CE","HE3","sing",False]]
    #termini
    if "C-term" in terminal:
        aa_connectivity.append(["C","OXT","sing",False])
        aa_connectivity.append(["OXT","HXT","sing",False])
        aa_connectivity.append(["C","OT1","doub",False])
        aa_connectivity.append(["C","OT2","sing",False])
        aa_connectivity.append(["OT2","HT2","sing",False])
    elif "C-cap" in terminal:
        aa_connectivity.append(["C","H","sing",False])
    if "N-term" in terminal:
        aa_connectivity.append(["N","H2","sing",False])
        if residue_name.upper() != "PRO": aa_connectivity.append(["N","H3","sing",False])
    elif "N-cap" in terminal:
        aa_connectivity.append(["N","C*","sing",False])
        aa_connectivity.append(["C*","O*","doub",False])
        aa_connectivity.append(["C*","H","sing",False])
    elif "ipa_termini" in terminal:
        aa_connectivity.append(["N","C#","sing",False])
        aa_connectivity.append(["C#","O#","doub",False])
    return aa_connectivity

def PTMConnectivity(residue_name, terminal):
    #function returning the connectivity of a given residue, PTM version
    aa_connectivity = []
    #basic PTM connectivity data
    #phosphorylation
    if residue_name.upper() == "SEP":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","OG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["OG","P","sing",False],
        ["P","O1P","doub",False],
        ["P","O2P","sing",False],
        ["P","O3P","sing",False]]
    elif residue_name.upper() == "TPO":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","OG1","sing",False],
        ["CB","CG2","sing",False],
        ["CB","HB","sing",False],
        ["CG2","HG21","sing",False],
        ["CG2","HG22","sing",False],
        ["CG2","HG23","sing",False],
        ["OG1","P","sing",False],
        ["P","O1P","doub",False],
        ["P","O2P","sing",False],
        ["P","O3P","sing",False]]
    elif residue_name.upper() == "PTR":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","CD1","doub",True],
        ["CG","CD2","sing",True],
        ["CD1","CE1","sing",True],
        ["CD1","HD1","sing",False],
        ["CD2","CE2","doub",True],
        ["CD2","HD2","sing",False],
        ["CE1","CZ","doub",True],
        ["CE1","HE1","sing",False],
        ["CE2","CZ","sing",True],
        ["CE2","HE2","sing",False],
        ["CZ","OH","sing",False],
        ["OH","P","sing",False],
        ["P","O1P","doub",False],
        ["P","O2P","sing",False],
        ["P","O3P","sing",False]]
    elif residue_name.upper() == "NEP":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","ND1","sing",True],
        ["CG","CD2","doub",True],
        ["ND1","CE1","doub",True],
        ["ND1","HD1","sing",False],
        ["CD2","NE2","sing",True],
        ["CD2","HD2","sing",False],
        ["CE1","NE2","sing",True],
        ["CE1","HE1","sing",False],
        ["NE2","P","sing",False],
        ["P","O1P","doub",False],
        ["P","O2P","sing",False],
        ["P","O3P","sing",False]]
    #acylation
    elif residue_name.upper() == "ALY":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","CD","sing",False],
        ["CG","HG2","sing",False],
        ["CG","HG3","sing",False],
        ["CD","CE","sing",False],
        ["CD","HD2","sing",False],
        ["CD","HD3","sing",False],
        ["CE","NZ","sing",False],
        ["CE","HE2","sing",False],
        ["CE","HE3","sing",False],
        ["NZ","CH","sing",False],
        ["CH","OH","doub",False],
        ["CH","CH3","sing",False]
        ]
    elif residue_name.upper() == "AME":
        aa_connectivity = [["N","CA","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","SD","sing",False],
        ["CG","HG2","sing",False],
        ["CG","HG3","sing",False],
        ["SD","CE","sing",False],
        ["CE","HE1","sing",False],
        ["CE","HE2","sing",False],
        ["CE","HE3","sing",False],
        ["N","CT1","sing",False],
        ["CT1","CT2","sing",False],
        ["CT1","OT","doub",False]]
    elif residue_name.upper() == "FME":
        aa_connectivity = [["N","CA","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","SD","sing",False],
        ["CG","HG2","sing",False],
        ["CG","HG3","sing",False],
        ["SD","CE","sing",False],
        ["CE","HE1","sing",False],
        ["CE","HE2","sing",False],
        ["CE","HE3","sing",False],
        ["N","CN","sing",False],
        ["CN","O1","doub",False]]
    #methylation
    elif residue_name.upper() == "MLY":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","CD","sing",False],
        ["CG","HG2","sing",False],
        ["CG","HG3","sing",False],
        ["CD","CE","sing",False],
        ["CD","HD2","sing",False],
        ["CD","HD3","sing",False],
        ["CE","NZ","sing",False],
        ["CE","HE2","sing",False],
        ["CE","HE3","sing",False],
        ["NZ","CH1","sing",False],
        ["NZ","CH2","sing",False]]
    elif residue_name.upper() == "MLZ":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","CD","sing",False],
        ["CG","HG2","sing",False],
        ["CG","HG3","sing",False],
        ["CD","CE","sing",False],
        ["CD","HD2","sing",False],
        ["CD","HD3","sing",False],
        ["CE","NZ","sing",False],
        ["CE","HE2","sing",False],
        ["CE","HE3","sing",False],
        ["NZ","CM","sing",False]]
    elif residue_name.upper() == "M3L":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","CD","sing",False],
        ["CG","HG2","sing",False],
        ["CG","HG3","sing",False],
        ["CD","CE","sing",False],
        ["CD","HD2","sing",False],
        ["CD","HD3","sing",False],
        ["CE","NZ","sing",False],
        ["CE","HE2","sing",False],
        ["CE","HE3","sing",False],
        ["NZ","CM1","sing",False],
        ["NZ","CM2","sing",False],
        ["NZ","CM3","sing",False]]
    elif residue_name.upper() == "AGM":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","CD","sing",False],
        ["CG","HG2","sing",False],
        ["CG","HG3","sing",False],
        ["CD","NE1","sing",False],
        ["CD","CE2","sing",False],
        ["CD","HD3","sing",False],
        ["NE1","CZ","sing",False],
        ["NE1","HE","sing",False],
        ["CZ","NH1","sing",False],
        ["CZ","NH2","doub",False],
        ["NH1","HH11","sing",False],
        ["NH1","HH12","sing",False],
        ["NH2","HH21","sing",False],
        ["NH2","HH22","sing",False]]
    elif residue_name.upper() == "MGN":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB1","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB1","CG","sing",False],
        ["CA","CB2","sing",False],
        ["CG","CD","sing",False],
        ["CG","HG2","sing",False],
        ["CG","HG3","sing",False],
        ["CD","OE1","doub",False],
        ["CD","NE2","sing",False],
        ["NE2","HE21","sing",False],
        ["NE2","HE22","sing",False]]
    elif residue_name.upper() == "HIC":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","ND1","sing",True],
        ["CG","CD2","doub",True],
        ["ND1","CE1","doub",True],
        ["ND1","HD1","sing",False],
        ["CD2","NE2","sing",True],
        ["CD2","HD2","sing",False],
        ["CE1","NE2","sing",True],
        ["CE1","HE1","sing",False],
        ["NE2","CZ","sing",False]]
    elif residue_name.upper() == "MHS":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","ND1","sing",True],
        ["CG","CD2","doub",True],
        ["ND1","CE1","sing",True],
        ["ND1","CM","sing",False],
        ["CD2","NE2","sing",True],
        ["CD2","HD2","sing",False],
        ["CE1","NE2","doub",True],
        ["CE1","HE1","sing",False],
        ["NE2","HE2","sing",False]]
    elif residue_name.upper() == "A1I9G":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","CF","sing",False],
        ["CG","ND1","sing",True],
        ["CG","CD2","doub",True],
        ["ND1","CE1","doub",True],
        ["ND1","HD1","sing",False],
        ["CD2","NE2","sing",True],
        ["CD2","HD2","sing",False],
        ["CE1","NE2","sing",True],
        ["CE1","HE1","sing",False],
        ["NE2","HE2","sing",False]]
    elif residue_name.upper() == "SMC":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","SG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["SG","CS","sing",False]]
    elif residue_name.upper() == "MMO":
        aa_connectivity = [["N","CA","sing",False],
        ["N","CN","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","CD","sing",False],
        ["CG","HG2","sing",False],
        ["CG","HG3","sing",False],
        ["CD","NE","sing",False],
        ["CD","HD2","sing",False],
        ["CD","HD3","sing",False],
        ["NE","CZ","sing",False],
        ["NE","HE","sing",False],
        ["CZ","NH1","sing",False],
        ["CZ","NH2","doub",False],
        ["NH1","HH11","sing",False],
        ["NH1","HH12","sing",False],
        ["NH2","HH21","sing",False],
        ["NH2","HH22","sing",False]]
    #hydroxylation
    elif residue_name.upper() == "TRX":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","CD1","doub",True],
        ["CG","CD2","sing",True],
        ["CD1","NE1","sing",True],
        ["CD1","HD1","sing",False],
        ["CD2","CE2","doub",True],
        ["CD2","CE3","sing",True],
        ["NE1","CE2","sing",True],
        ["NE1","HE1","sing",False],
        ["CE2","CZ2","sing",True],
        ["CE3","CZ3","doub",True],
        ["CE3","HE3","sing",False],
        ["CZ2","CH2","doub",True],
        ["CZ2","HZ2","sing",False],
        ["CZ3","CH2","sing",True],
        ["CZ3","HZ3","sing",False],
        ["CH2","OH2","sing",False]]
    elif residue_name.upper() == "TOQ":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","CD1","doub",True],
        ["CG","CD2","sing",True],
        ["CD1","NE1","sing",True],
        ["CD1","HD1","sing",False],
        ["CD2","CE2","doub",True],
        ["CD2","CE3","sing",True],
        ["NE1","CE2","sing",True],
        ["NE1","HE1","sing",False],
        ["CE2","CZ2","sing",True],
        ["CE3","CZ3","doub",True],
        ["CE3","HE3","sing",False],
        ["CZ2","CH2","doub",True],
        ["CZ2","O7","sing",False],
        ["CZ3","CH2","sing",True],
        ["CZ3","HZ3","sing",False],
        ["CH2","O6","sing",False]]
    elif residue_name.upper() == "4HT":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","CD1","doub",True],
        ["CG","CD2","sing",True],
        ["CD1","NE1","sing",True],
        ["CD1","HD1","sing",False],
        ["CD2","CE2","doub",True],
        ["CD2","CE3","sing",True],
        ["NE1","CE2","sing",True],
        ["NE1","HE1","sing",False],
        ["CE2","CZ2","sing",True],
        ["CE3","CZ3","doub",True],
        ["CE3","OE3","sing",False],
        ["CZ2","CH2","doub",True],
        ["CZ2","HZ2","sing",False],
        ["CZ3","CH2","sing",True],
        ["CZ3","HZ3","sing",False],
        ["CH2","HH2","sing",False]]
    elif residue_name.upper() == "4PQ":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","CD1","doub",True],
        ["CG","CD2","sing",True],
        ["CD1","NE1","sing",True],
        ["CD1","HD1","sing",False],
        ["CD2","CE2","doub",True],
        ["CD2","CE3","sing",True],
        ["NE1","CE2","sing",True],
        ["NE1","HE1","sing",False],
        ["CE2","CZ2","sing",True],
        ["CE3","CZ3","doub",True],
        ["CE3","HE3","sing",False],
        ["CZ2","CH2","doub",True],
        ["CZ2","HZ2","sing",False],
        ["CZ3","CH2","sing",True],
        ["CZ3","OAD","sing",False],
        ["CH2","HH2","sing",False]]
    elif residue_name.upper() == "HRP":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","CD1","doub",True],
        ["CG","CD2","sing",True],
        ["CD1","NE1","sing",True],
        ["CD1","HD1","sing",False],
        ["CD2","CE2","doub",True],
        ["CD2","CE3","sing",True],
        ["NE1","CE2","sing",True],
        ["NE1","HE1","sing",False],
        ["CE2","CZ2","sing",True],
        ["CE3","CZ3","doub",True],
        ["CE3","HE3","sing",False],
        ["CZ2","CH2","doub",True],
        ["CZ2","HZ2","sing",False],
        ["CZ3","CH2","sing",True],
        ["CZ3","OX","sing",False],
        ["CH2","HH2","sing",False]]
    elif residue_name.upper() == "HTR":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","OH","sing",False],
        ["CB","HB3","sing",False],
        ["CG","CD1","doub",True],
        ["CG","CD2","sing",True],
        ["CD1","NE1","sing",True],
        ["CD1","HD1","sing",False],
        ["CD2","CE2","doub",True],
        ["CD2","CE3","sing",True],
        ["NE1","CE2","sing",True],
        ["NE1","HE1","sing",False],
        ["CE2","CZ2","sing",True],
        ["CE3","CZ3","doub",True],
        ["CE3","HE3","sing",False],
        ["CZ2","CH2","doub",True],
        ["CZ2","HZ2","sing",False],
        ["CZ3","CH2","sing",True],
        ["CZ3","HZ3","sing",False],
        ["CH2","HH2","sing",False]]
    elif residue_name.upper() == "UX8":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","O9","sing",False],
        ["CB","HB3","sing",False],
        ["CG","CD1","doub",True],
        ["CG","CD2","sing",True],
        ["CD1","NE1","sing",True],
        ["CD1","HD1","sing",False],
        ["CD2","CE2","doub",True],
        ["CD2","CE3","sing",True],
        ["NE1","CE2","sing",True],
        ["NE1","HE1","sing",False],
        ["CE2","CZ2","sing",True],
        ["CE3","CZ3","doub",True],
        ["CE3","HE3","sing",False],
        ["CZ2","CH2","doub",True],
        ["CZ2","HZ2","sing",False],
        ["CZ3","CH2","sing",True],
        ["CZ3","HZ3","sing",False],
        ["CH2","HH2","sing",False]]
    elif residue_name.upper() == "HYP":
        aa_connectivity = [["N","CA","sing",False],
        ["N","CD","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","CD","sing",False],
        ["CG","OD1","sing",False],
        ["CG","HG3","sing",False],
        ["CD","HD2","sing",False],
        ["CD","HD3","sing",False]]
    elif residue_name.upper() == "BHD":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","OB","sing",False],
        ["CG","OD1","doub",False],
        ["CG","OD2","sing",False],
        ["OD2","HD2","sing",False]]
    elif residue_name.upper() == "BH2":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","OB","sing",False],
        ["CG","OD1","doub",False],
        ["CG","OD2","sing",False],
        ["OD2","HD2","sing",False]]
    elif residue_name.upper() == "CSO":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","SG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["SG","OD","sing",False]]
    elif residue_name.upper() == "TYQ":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","CD1","doub",True],
        ["CG","CD2","sing",True],
        ["CD1","CE1","sing",True],
        ["CD1","OZ","sing",False],
        ["CD2","CE2","doub",True],
        ["CD2","HD2","sing",False],
        ["CE1","CZ","doub",True],
        ["CE1","HE1","sing",False],
        ["CE2","CZ","sing",True],
        ["CE2","N5","sing",False],
        ["CZ","OH","sing",False],
        ["OH","HH","sing",False]]
    elif residue_name.upper() == "LYZ":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","CD","sing",False],
        ["CG","HG2","sing",False],
        ["CG","HG3","sing",False],
        ["CD","CE","sing",False],
        ["CD","HD2","sing",False],
        ["CD","OH","sing",False],
        ["CE","NZ","sing",False],
        ["CE","HE2","sing",False],
        ["CE","HE3","sing",False],
        ["NZ","HZ1","sing",False],
        ["NZ","HZ2","sing",False],
        ["NZ","HZ3","sing",False]]
    elif residue_name.upper() == "ILX":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG1","sing",False],
        ["CB","CG2","sing",False],
        ["CB","HB","sing",False],
        ["CG1","CD1","sing",False],
        ["CG1","OG1","sing",False],
        ["CG1","HG13","sing",False],
        ["CG2","HG21","sing",False],
        ["CG2","HG22","sing",False],
        ["CG2","HG23","sing",False],
        ["CD1","HD11","sing",False],
        ["CD1","HD12","sing",False],
        ["CD1","OD1","sing",False]]
    #carboxylation
    elif residue_name.upper() == "KCX":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","CD","sing",False],
        ["CG","HG2","sing",False],
        ["CG","HG3","sing",False],
        ["CD","CE","sing",False],
        ["CD","HD2","sing",False],
        ["CD","HD3","sing",False],
        ["CE","NZ","sing",False],
        ["CE","HE2","sing",False],
        ["CE","HE3","sing",False],
        ["NZ","HZ1","sing",False],
        ["NZ","CX","sing",False],
        ["CX","OQ1","doub",False],
        ["CX","OQ2","sing",False]]
    #cyclisation
    elif residue_name.upper() == "PCA":
        aa_connectivity = [["N","CA","sing",False],
        ["N","CD","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","CD","sing",False],
        ["CG","HG2","sing",False],
        ["CG","HG3","sing",False],
        ["CD","OE","doub",False]]
    elif residue_name.upper() == "7YO":
        aa_connectivity = [["N","CA","sing",False],
        ["N","CD","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","CG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["CG","CD","sing",False],
        ["CG","HG2","sing",False],
        ["CG","HG3","sing",False],
        ["CD","O1","doub",False]]
    #other
    elif residue_name.upper() == "SCH":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","SG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["SG","SD","sing",False],
        ["SD","CE","sing",False]]
    elif residue_name.upper() == "2CO":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","SG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["SG","OD","sing",False],
        ["OD","OE","sing",False]]
    elif residue_name.upper() == "CSD":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","SG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["SG","OD1","sing",False],
        ["SG","OD2","doub",False]]
    elif residue_name.upper() == "0CS":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","SG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["SG","OD1","doub",False],
        ["SG","OD2","sing",False],
        ["OD2","OD3","sing",False]]
    elif residue_name.upper() == "CSX":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","SG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["SG","OD","doub",False]]
    elif residue_name.upper() == "CSS":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","SG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["SG","SD","sing",False]]
    elif residue_name.upper() == "TSY":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","SG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["SG","S1","sing",False],
        ["S1","S2","sing",False]]
    elif residue_name.upper() == "TQZ":
        aa_connectivity = [["N","CA","sing",False],
        ["N","H","sing",False],
        ["CA","C","sing",False],
        ["CA","CB","sing",False],
        ["CA","HA","sing",False],
        ["C","O","doub",False],
        ["CB","SG","sing",False],
        ["CB","HB2","sing",False],
        ["CB","HB3","sing",False],
        ["SG","S1","sing",False],
        ["S1","S2","sing",False],
        ["S2","S3","sing",False],
        ["S3","S4","sing",False]]
    #termini
    if "C-term" in terminal:
        aa_connectivity.append(["C","OXT","sing",False])
        aa_connectivity.append(["OXT","HXT","sing",False])
        aa_connectivity.append(["C","OT1","doub",False])
        aa_connectivity.append(["C","OT2","sing",False])
        aa_connectivity.append(["OT2","HT2","sing",False])
    elif "C-cap" in terminal:
        aa_connectivity.append(["C","H","sing",False])
    if "N-term" in terminal:
        aa_connectivity.append(["N","H2","sing",False])
        if residue_name.upper() != "PRO": aa_connectivity.append(["N","H3","sing",False])
    elif "N-cap" in terminal:
        aa_connectivity.append(["N","C*","sing",False])
        aa_connectivity.append(["C*","O*","doub",False])
        aa_connectivity.append(["C*","H","sing",False])
    elif "ipa_termini" in terminal:
        aa_connectivity.append(["N","C#","sing",False])
        aa_connectivity.append(["C#","O#","doub",False])
    return aa_connectivity

def AminoAcidConnectivity(residue_name, terminal):
    #function returning the connectivity of a given residue
    aa_connectivity = []
    #basic connectivity
    if str(residue_name).strip().upper() in _BASE_AA_RESNAMES: aa_connectivity = BaseAminoAcidConnectivity(residue_name,terminal)
    elif str(residue_name).strip().upper() in _EXTENDED_AA_RESNAMES: aa_connectivity = ExtendedAminoAcidConnectivity(residue_name,terminal)
    elif str(residue_name).strip().upper() in _BASE_PTM_RESNAMES: aa_connectivity = PTMConnectivity(residue_name,terminal)
    #extra connectivity, such as connectivity to the previous residue
    aa_connectivity.extend([
        ["Cp","N","sing",False],
        ["Cp","Op","doub",False],
        ["Cp","CAp","sing",False],
        ["Nf","C","sing",False],
        ["Nf","CAf","sing",False]])
    return aa_connectivity

# ------------------------------------------------------------
# Charge-related functions
# ------------------------------------------------------------
def GetLocalCharges(residue_name, atom_types):
    #function that returns local charges from specific residues
    local_charges = [0]*(len(atom_types) + 1)
    if residue_name == "M3L":
        for iatm in range(len(atom_types)):
            if atom_types[iatm] == "NZ":
                local_charges[iatm + 1] = 1
                local_charges[0] = 1
                break
    elif residue_name == "NAP":
        for iatm in range(len(atom_types)):
            if atom_types[iatm] == "N1N":
                local_charges[iatm + 1] = 1
                local_charges[0] = 1
                break
    elif residue_name == "F43":
        for iatm in range(len(atom_types)):
            if atom_types[iatm] == "NI":
                local_charges[iatm + 1] = 1
                local_charges[0] += 1
            elif atom_types[iatm] == "NC":
                local_charges[iatm + 1] = -1
                local_charges[0] += 1
            if local_charges[0] == 2: break
    return local_charges

# ------------------------------------------------------------
# Other functions
# ------------------------------------------------------------
def main_chain_atom(atom_name):
    #function checking whether an atom is part of a main chain
    main_chain = False
    if (atom_name == "C") or (atom_name == "CA") or (atom_name == "N") or (atom_name == "O"): main_chain = True
    return main_chain

__all__ = [
    "IsBaseAminoAcid",
    "IsExtendedAminoAcid",
    "IsPostTranslationalModification",
    "ExternalCapResnames",
    "IsExternalCap",
    "BaseAminoAcidConnectivity",
    "ExtendedAminoAcidConnectivity",
    "PTMConnectivity",
    "AminoAcidConnectivity",
    "GetLocalCharges",
    "main_chain_atom"
]
