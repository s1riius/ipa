# ------------------------------------------------------------
# Functions to support reesidue-flipping
# ------------------------------------------------------------

def FlipASN(atoms_comb, atp_comb, index):
    #function to flip ASNs
    #swap OD1 with ND2 or OD2 with ND1
    for iatm in range(len(atoms_comb[index])):
        if atp_comb[index][iatm] == "OD1":
            atp_comb[index][iatm] = "ND2"
            atoms_comb[index][iatm] = "N"
        elif atp_comb[index][iatm] == "OD2":
            atp_comb[index][iatm] = "ND1"
            atoms_comb[index][iatm] = "N"
        elif atp_comb[index][iatm] == "ND1":
            atp_comb[index][iatm] = "OD2"
            atoms_comb[index][iatm] = "O"
        elif atp_comb[index][iatm] == "ND2":
            atp_comb[index][iatm] = "OD1"
            atoms_comb[index][iatm] = "O"
    return 0

def FlipGLN(atoms_comb, atp_comb, index):
    #function flipping GLNs
    #swap OE1 with NE2 or OE2 with NE1
    for iatm in range(len(atoms_comb[index])):
        if atp_comb[index][iatm] == "OE1":
            atp_comb[index][iatm] = "NE2"
            atoms_comb[index][iatm] = "N"
        elif atp_comb[index][iatm] == "OE2":
            atp_comb[index][iatm] = "NE1"
            atoms_comb[index][iatm] = "N"
        elif atp_comb[index][iatm] == "NE1":
            atp_comb[index][iatm] = "OE2"
            atoms_comb[index][iatm] = "O"
        elif atp_comb[index][iatm] == "NE2":
            atp_comb[index][iatm] = "OE1"
            atoms_comb[index][iatm] = "O"
    return 0

def FlipHIS(atoms_comb, atp_comb, index):
    #function flipping HISs
    #swap ND1 with CD2 or CD1 with ND2
    #swap NE2 with CE1 or CE2 with NE1
    for iatm in range(len(atoms_comb[index])):
        if atp_comb[index][iatm] == "ND1":
            atp_comb[index][iatm] = "CD2"
            atoms_comb[index][iatm] = "C"
        elif atp_comb[index][iatm] == "ND2":
            atp_comb[index][iatm] = "CD1"
            atoms_comb[index][iatm] = "C"
        elif atp_comb[index][iatm] == "CD1":
            atp_comb[index][iatm] = "ND2"
            atoms_comb[index][iatm] = "N"
        elif atp_comb[index][iatm] == "CD2":
            atp_comb[index][iatm] = "ND1"
            atoms_comb[index][iatm] = "N"
        elif atp_comb[index][iatm] == "NE1":
            atp_comb[index][iatm] = "CE2"
            atoms_comb[index][iatm] = "C"
        elif atp_comb[index][iatm] == "NE2":
            atp_comb[index][iatm] = "CE1"
            atoms_comb[index][iatm] = "C"
        elif atp_comb[index][iatm] == "CE1":
            atp_comb[index][iatm] = "NE2"
            atoms_comb[index][iatm] = "N"
        elif atp_comb[index][iatm] == "CE2":
            atp_comb[index][iatm] = "NE1"
            atoms_comb[index][iatm] = "N"
    return 0
