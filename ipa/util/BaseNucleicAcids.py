# ------------------------------------------------------------
# Functions related to base nucleic acid data handling
# ------------------------------------------------------------

import util.Geometry as GEOM
import numpy as np

# ------------------------------------------------------------
# Base nucleic acid handlers
# ------------------------------------------------------------
#canonical letter residue names for standard nucleic acids plus common variants
_BASE_NA_RESNAMES = {
    #base RNA
    "A","C","G","U",
    #base DNA
    "DA","DC","DG","DT",
    #other, related residues
    "DU","DI","I"
}

_EXTENDED_NA_RESNAMES = {
    "PSU","H2U","5MC","OMC","OMG",
    "7MG","1MA","M2G","YYG","5MU",
    "UR3","5CM","8OG","CBR","BRU",
    "2MG","PYO","LCC","FHU","TFO"
}

def IsBaseNucleicAcid(residue):
    """Return True if `residue` is recognised as a base amino-acid residue name."""
    return str(residue).strip().upper() in _BASE_NA_RESNAMES

def IsExtendedNucleicAcid(residue):
    """Return True if `residue` is recognised as a base amino-acid residue name."""
    return str(residue).strip().upper() in _EXTENDED_NA_RESNAMES

# ------------------------------------------------------------
# Nucleic acid and derivativatives connectivity functions
# ------------------------------------------------------------
def BaseNucleicAcidConnectivity(residue_name, terminal):
    #function returning the connectivity of a given residue
    na_connectivity = []
    #basic amino acid connectivity data
    if residue_name.upper() == "A":
        na_connectivity = [["OP3","P","sing",False],
                           ["P","OP1","doub",False],
                           ["P","OP2","sing",False],
                           ["P","O5'","sing",False],
                           ["O5'","C5'","sing",False],
                           ["C5'","C4'","sing",False],
                           ["C4'","O4'","sing",False],
                           ["C4'","C3'","sing",False],
                           ["O4'","C1'","sing",False],
                           ["C3'","O3'","sing",False],
                           ["C3'","C2'","sing",False],
                           ["C2'","O2'","sing",False],
                           ["C2'","C1'","sing",False],
                           ["C1'","N9","sing",False],
                           ["N9","C8","sing",False],
                           ["N9","C4","sing",False],
                           ["C8","N7","doub",False],
                           ["N7","C5","sing",False],
                           ["C5","C6","sing",False],
                           ["C5","C4","doub",False],
                           ["C6","N6","sing",False],
                           ["C6","N1","doub",False],
                           ["N1","C2","sing",False],
                           ["C2","N3","doub",False],
                           ["N3","C4","sing",False]]
    elif residue_name.upper() == "DA":
        na_connectivity = [["OP3","P","sing",False],
                           ["P","OP1","doub",False],
                           ["P","OP2","sing",False],
                           ["P","O5'","sing",False],
                           ["O5'","C5'","sing",False],
                           ["C5'","C4'","sing",False],
                           ["C4'","O4'","sing",False],
                           ["C4'","C3'","sing",False],
                           ["O4'","C1'","sing",False],
                           ["C3'","O3'","sing",False],
                           ["C3'","C2'","sing",False],
                           ["C2'","C1'","sing",False],
                           ["C1'","N9","sing",False],
                           ["N9","C8","sing",False],
                           ["N9","C4","sing",False],
                           ["C8","N7","doub",False],
                           ["N7","C5","sing",False],
                           ["C5","C6","sing",False],
                           ["C5","C4","doub",False],
                           ["C6","N6","sing",False],
                           ["C6","N1","doub",False],
                           ["N1","C2","sing",False],
                           ["C2","N3","doub",False],
                           ["N3","C4","sing",False]]
    elif residue_name.upper() == "C":
        na_connectivity = [["OP3","P","sing",False],
                          ["P","OP1","doub",False],
                          ["P","OP2","sing",False],
                          ["P","O5'","sing",False],
                          ["O5'","C5'","sing",False],
                          ["C5'","C4'","sing",False],
                          ["C4'","O4'","sing",False],
                          ["C4'","C3'","sing",False],
                          ["O4'","C1'","sing",False],
                          ["C3'","O3'","sing",False],
                          ["C3'","C2'","sing",False],
                          ["C2'","O2'","sing",False],
                          ["C2'","C1'","sing",False],
                          ["C1'","N1","sing",False],
                          ["N1","C2","sing",False],
                          ["N1","C6","sing",False],
                          ["C2","O2","doub",False],
                          ["C2","N3","sing",False],
                          ["N3","C4","doub",False],
                          ["C4","N4","sing",False],
                          ["C4","C5","sing",False],
                          ["C5","C6","doub",False]]
    elif residue_name.upper() == "DC":
        na_connectivity = [["OP3","P","sing",False],
                           ["P","OP1","doub",False],
                           ["P","OP2","sing",False],
                           ["P","O5'","sing",False],
                           ["O5'","C5'","sing",False],
                           ["C5'","C4'","sing",False],
                           ["C4'","O4'","sing",False],
                           ["C4'","C3'","sing",False],
                           ["O4'","C1'","sing",False],
                           ["C3'","O3'","sing",False],
                           ["C3'","C2'","sing",False],
                           ["C2'","C1'","sing",False],
                           ["C1'","N1","sing",False],
                           ["N1","C2","sing",False],
                           ["N1","C6","sing",False],
                           ["C2","O2","doub",False],
                           ["C2","N3","sing",False],
                           ["N3","C4","doub",False],
                           ["C4","N4","sing",False],
                           ["C4","C5","sing",False],
                           ["C5","C6","doub",False]]
    elif residue_name.upper() == "G":
        na_connectivity = [["OP3","P","sing",False],
                           ["P","OP1","doub",False],
                           ["P","OP2","sing",False],
                           ["P","O5'","sing",False],
                           ["O5'","C5'","sing",False],
                           ["C5'","C4'","sing",False],
                           ["C4'","O4'","sing",False],
                           ["C4'","C3'","sing",False],
                           ["O4'","C1'","sing",False],
                           ["C3'","O3'","sing",False],
                           ["C3'","C2'","sing",False],
                           ["C2'","O2'","sing",False],
                           ["C2'","C1'","sing",False],
                           ["C1'","N9","sing",False],
                           ["N9","C8","sing",False],
                           ["N9","C4","sing",False],
                           ["C8","N7","doub",False],
                           ["N7","C5","sing",False],
                           ["C5","C6","sing",False],
                           ["C5","C4","doub",False],
                           ["C6","O6","doub",False],
                           ["C6","N1","sing",False],
                           ["N1","C2","sing",False],
                           ["C2","N2","sing",False],
                           ["C2","N3","doub",False],
                           ["N3","C4","sing",False]]
    elif residue_name.upper() == "DG":
        na_connectivity = [["OP3","P","sing",False],
                           ["P","OP1","doub",False],
                           ["P","OP2","sing",False],
                           ["P","O5'","sing",False],
                           ["O5'","C5'","sing",False],
                           ["C5'","C4'","sing",False],
                           ["C4'","O4'","sing",False],
                           ["C4'","C3'","sing",False],
                           ["O4'","C1'","sing",False],
                           ["C3'","O3'","sing",False],
                           ["C3'","C2'","sing",False],
                           ["C2'","C1'","sing",False],
                           ["C1'","N9","sing",False],
                           ["N9","C8","sing",False],
                           ["N9","C4","sing",False],
                           ["C8","N7","doub",False],
                           ["N7","C5","sing",False],
                           ["C5","C6","sing",False],
                           ["C5","C4","doub",False],
                           ["C6","O6","doub",False],
                           ["C6","N1","sing",False],
                           ["N1","C2","sing",False],
                           ["C2","N2","sing",False],
                           ["C2","N3","doub",False],
                           ["N3","C4","sing",False]]
    elif residue_name.upper() == "I":
        na_connectivity = [["OP3","P","sing",False],
                           ["P","OP1","doub",False],
                           ["P","OP2","sing",False],
                           ["P","O5'","sing",False],
                           ["O5'","C5'","sing",False],
                           ["C5'","C4'","sing",False],
                           ["C4'","O4'","sing",False],
                           ["C4'","C3'","sing",False],
                           ["O4'","C1'","sing",False],
                           ["C3'","O3'","sing",False],
                           ["C3'","C2'","sing",False],
                           ["C2'","O2'","sing",False],
                           ["C2'","C1'","sing",False],
                           ["C1'","N9","sing",False],
                           ["N9","C8","sing",False],
                           ["N9","C4","sing",False],
                           ["C8","N7","doub",False],
                           ["N7","C5","sing",False],
                           ["C5","C6","sing",False],
                           ["C5","C4","doub",False],
                           ["C6","O6","doub",False],
                           ["C6","N1","sing",False],
                           ["N1","C2","sing",False],
                           ["C2","N3","doub",False],
                           ["N3","C4","sing",False]]
    elif residue_name.upper() == "DI":
        na_connectivity = [["OP3","P","sing",False],
                           ["P","OP1","doub",False],
                           ["P","OP2","sing",False],
                           ["P","O5'","sing",False],
                           ["O5'","C5'","sing",False],
                           ["C5'","C4'","sing",False],
                           ["C4'","O4'","sing",False],
                           ["C4'","C3'","sing",False],
                           ["O4'","C1'","sing",False],
                           ["C3'","O3'","sing",False],
                           ["C3'","C2'","sing",False],
                           ["C2'","C1'","sing",False],
                           ["C1'","N9","sing",False],
                           ["N9","C8","sing",False],
                           ["N9","C4","sing",False],
                           ["C8","N7","doub",False],
                           ["N7","C5","sing",False],
                           ["C5","C6","sing",False],
                           ["C5","C4","doub",False],
                           ["C6","O6","doub",False],
                           ["C6","N1","sing",False],
                           ["N1","C2","sing",False],
                           ["C2","N3","doub",False],
                           ["N3","C4","sing",False]]
    elif residue_name.upper() == "U":
        na_connectivity = [["OP3","P","sing",False],
                           ["P","OP1","doub",False],
                           ["P","OP2","sing",False],
                           ["P","O5'","sing",False],
                           ["O5'","C5'","sing",False],
                           ["C5'","C4'","sing",False],
                           ["C4'","O4'","sing",False],
                           ["C4'","C3'","sing",False],
                           ["O4'","C1'","sing",False],
                           ["C3'","O3'","sing",False],
                           ["C3'","C2'","sing",False],
                           ["C2'","O2'","sing",False],
                           ["C2'","C1'","sing",False],
                           ["C1'","N1","sing",False],
                           ["N1","C2","sing",False],
                           ["N1","C6","sing",False],
                           ["C2","O2","doub",False],
                           ["C2","N3","sing",False],
                           ["N3","C4","sing",False],
                           ["C4","O4","doub",False],
                           ["C4","C5","sing",False],
                           ["C5","C6","doub",False]]
    elif residue_name.upper() == "DU":
        na_connectivity = [["OP3","P","sing",False],
                           ["P","OP1","doub",False],
                           ["P","OP2","sing",False],
                           ["P","O5'","sing",False],
                           ["O5'","C5'","sing",False],
                           ["C5'","C4'","sing",False],
                           ["C4'","O4'","sing",False],
                           ["C4'","C3'","sing",False],
                           ["O4'","C1'","sing",False],
                           ["C3'","O3'","sing",False],
                           ["C3'","C2'","sing",False],
                           ["C2'","C1'","sing",False],
                           ["C1'","N1","sing",False],
                           ["N1","C2","sing",False],
                           ["N1","C6","sing",False],
                           ["C2","O2","doub",False],
                           ["C2","N3","sing",False],
                           ["N3","C4","sing",False],
                           ["C4","O4","doub",False],
                           ["C4","C5","sing",False],
                           ["C5","C6","doub",False]]
    elif residue_name.upper() == "DT":
        na_connectivity = [["OP3","P","sing",False],
                           ["P","OP1","doub",False],
                           ["P","OP2","sing",False],
                           ["P","O5'","sing",False],
                           ["O5'","C5'","sing",False],
                           ["C5'","C4'","sing",False],
                           ["C4'","O4'","sing",False],
                           ["C4'","C3'","sing",False],
                           ["O4'","C1'","sing",False],
                           ["C3'","O3'","sing",False],
                           ["C3'","C2'","sing",False],
                           ["C2'","C1'","sing",False],
                           ["C1'","N1","sing",False],
                           ["N1","C2","sing",False],
                           ["N1","C6","sing",False],
                           ["C2","O2","doub",False],
                           ["C2","N3","sing",False],
                           ["N3","C4","sing",False],
                           ["C4","O4","doub",False],
                           ["C4","C5","sing",False],
                           ["C5","C7","sing",False],
                           ["C5","C6","doub",False]]
    #termini
    if "ipa_termini" in terminal:
        na_connectivity.append(["P","O#","sing",False])
        na_connectivity.append(["C#","O#","sing",False])
    return na_connectivity

def ExtendedNucleicAcidConnectivity(residue_name, terminal):
    #function returning the connectivity of a given residue
    na_connectivity = []#termini
    if "ipa_termini" in terminal:
        na_connectivity.append(["P","O#","sing",False])
        na_connectivity.append(["C#","O#","sing",False])
    return na_connectivity

# ------------------------------------------------------------
# Functions to remove bad protonation states
# ------------------------------------------------------------
def get_phosphate_charge(path2sdf):
    """A function that uses geometry considerations to identify the phosphate charge"""
    p_charge = 0
    n_oxygen_neighbour = 0
    #get geometry and connectivity
    connectivity,atomtp,atoms,_,_ = GEOM.ReadSDFWithTopology(path2sdf)
    #find phosphorous
    ipx = -1
    p_type = ""
    for iatm in range(len(atoms)):
        if atoms[iatm] == "P":
            ipx = iatm
            p_type = atomtp[iatm]
            break
    if ipx >= 0:
        #now look for its neighbouring oxygens
        for icnn in connectivity:
            if (icnn[0] == p_type) or (icnn[1] == p_type):
                o_type = icnn[1] if icnn[0] == p_type else icnn[0]
                for iatm in range(len(atomtp)):
                    if atomtp[iatm] == o_type and atoms[iatm] == "O":
                        n_oxygen_neighbour += 1
                        if icnn[2] == "sing":
                            nneighbours = 0
                            for jcnn in connectivity:
                                if (jcnn[0] == o_type) or (jcnn[1] == o_type): nneighbours += 1
                            if nneighbours == 1: p_charge -= 1
                        break
    else: p_charge = None
    return p_charge,n_oxygen_neighbour

def fix_phosphate_incomplete_octet(path2sdf, target_distance = 1.488):
    """Function that adds missing oxygen atom to a terminal phosphate"""
    connectivity,atomtp,atoms,geom,_ = GEOM.ReadSDFWithTopology(path2sdf)
    #find the phosphorous
    p_idx = -1
    p_atom_tp = ""
    coordinates_O = [0.0,0.0,0.0]
    idx_replace = -1
    old_element = ""
    new_element = ""
    for iatm in range(len(atoms)):
        if atoms[iatm] == "P":
            p_idx = iatm
            p_atom_tp = atomtp[iatm]
            break
    if p_idx >= 0:
        #find the P-X, X!=O bond
        for icnn in connectivity:
            if (icnn[0] == p_atom_tp) or (icnn[1] == p_atom_tp):
                the_other_type = icnn[1] if icnn[0] == p_atom_tp else icnn[0]
                for iatm in range(len(atomtp)):
                    if atomtp[iatm] == the_other_type:
                        idx_replace = iatm
                        break
                if (idx_replace >= 0) and (atoms[idx_replace] == "H"):
                    #fix the bond
                    Pcoord = np.asarray(geom[p_idx])
                    Hcoord = np.asarray(geom[idx_replace])
                    direction = Pcoord - Hcoord
                    norm = np.linalg.norm(direction)
                    direction /= norm
                    Ocoord = Pcoord - target_distance*direction
                    coordinates_O = list(Ocoord)
                    old_element = "H"
                    new_element = "O"
    return idx_replace,coordinates_O,old_element,new_element

# ------------------------------------------------------------
# Functions to add protons to phosphate
# ------------------------------------------------------------
def place_H_from_dihedral(Oref, P, Otarget, OH = 0.96, angle_deg = 108.5, dihedral_deg = 180.0):
    """Place H using internal coordinates:
           Oref - P - Otarget - H
       Parameters
       ----------
       Oref : xyz of reference oxygen
       P : xyz of phosphorus
       Otarget : xyz of oxygen receiving proton
       OH : O-H distance
       angle_deg : P-O-H angle
       dihedral_deg : Oref-P-O-H dihedral
       Returns
       -------
       H : xyz coordinate of proton"""

    def unit(v):
        n = np.linalg.norm(v)
        if n < 1e-8:
            raise ValueError("zero-length vector")
        return v/n
    
    Oref = np.asarray(Oref,dtype = float)
    P = np.asarray(P,dtype = float)
    Otarget = np.asarray(Otarget,dtype = float)
    theta = np.deg2rad(angle_deg)
    phi = np.deg2rad(dihedral_deg)
    #local frame at Otarget
    e1 = unit(P - Otarget)  #Otarget -> P
    #normal to plane Oref-P-Otarget
    n = unit(np.cross(Oref - P,Otarget - P))
    #perpendicular direction in the Otarget-centered frame
    e2 = unit(np.cross(n,e1))
    e3 = unit(np.cross(e1,e2))
    #build H direction
    hdir = (np.cos(theta)*e1 + np.sin(theta)*(np.cos(phi)*e2 + np.sin(phi)*e3))
    H_coord = Otarget + OH*unit(hdir)
    return H_coord

__all__ = [
    "IsBaseNucleicAcid",
    "IsExtendedNucleicAcid",
    "BaseNucleicAcidConnectivity",
    "get_phosphate_charge",
    "fix_phosphate_incomplete_octet",
    "place_H_from_dihedral"
]
