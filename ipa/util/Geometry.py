# ------------------------------------------------------------
# Geometry-related utilitarians
# ------------------------------------------------------------

import gemmi
import json
import subprocess
import util.Formatting as FORMAT

# ------------------------------------------------------------
# General functions
# ------------------------------------------------------------
def AtomFix(atom):
    #function that fixes element names in PDB files
    element = atom
    if atom == "CL": element = "Cl"
    elif atom == "NA": element = "Na"
    elif atom == "MG": element = "Mg"
    elif atom == "FE": element = "Fe"
    return element

# ------------------------------------------------------------
# XYZ-related functions
# ------------------------------------------------------------
def Write2XYZ(atmlist, geometry, name):
    #function writing an xyz file given geometry and atom list
    if not name.endswith(".xyz"): name += ".xyz"
    wfile = open(name,"w")
    wfile.write(str(len(atmlist)) + "\n\n")
    for idatm in range(len(atmlist)):
        element = AtomFix(atmlist[idatm])
        wfile.write(element + "   " + str(geometry[idatm][0]) + "    " + str(geometry[idatm][1]) + "    " + str(geometry[idatm][2]) + "\n")
    wfile.close()
    return 0

def ReadXYZ(path2file):
    #function reading an xyz and writing to containers
    if not path2file.endswith(".xyz"): path2file += ".xyz"
    rfile = open(path2file,"r")
    rfilelines = rfile.readlines()
    rfile.close()
    natms = int(rfilelines[0])
    atoms = []
    geometry = []
    for iatm in range(natms):
        line = rfilelines[iatm + 2].strip().rstrip()
        cleanline = FORMAT.CleanTheLine(line).split(" ")
        atoms.append(cleanline[0])
        geometry.append([float(cleanline[1]),float(cleanline[2]),float(cleanline[3])])
    return atoms,geometry

def MergeXYZ(xyz1, xyz2, newname):
    #function that merges two geometry files into a new xyz file
    rfile = open(xyz1,"r")
    rfilelines = rfile.readlines()
    rfile.close()
    natm1 = int(rfilelines[0])
    atoms = []
    geometry = []
    for idatm in range(2,natm1 + 2):
        line = rfilelines[idatm]
        cleanline = FORMAT.CleanTheLine(line)
        data = cleanline.split(" ")
        atoms.append(data[0])
        geometry.append([float(data[1]),float(data[2]),float(data[3])])
    rfile = open(xyz2,"r")
    rfilelines = rfile.readlines()
    rfile.close()
    natm2 = int(rfilelines[0])
    for idatm in range(2,natm2 + 2):
        line = rfilelines[idatm]
        cleanline = FORMAT.CleanTheLine(line)
        data = cleanline.split(" ")
        atoms.append(data[0])
        geometry.append([float(data[1]),float(data[2]),float(data[3])])
    Write2XYZ(atoms,geometry,newname)

# ------------------------------------------------------------
# SDF-related functions
# ------------------------------------------------------------
def ReadSDF(path2file):
    #function reading an sdf and writing to containers
    if not path2file.endswith(".sdf"): path2file += ".sdf"
    rfile = open(path2file,"r")
    rfilelines = rfile.readlines()
    natms = int(rfilelines[3][0:3])
    atoms = []
    geometry = []
    for iatm in range(natms):
        line = rfilelines[iatm + 4].strip().rstrip()
        cleanline = FORMAT.CleanTheLine(line).split(" ")
        atoms.append(cleanline[3])
        geometry.append([float(cleanline[0]),float(cleanline[1]),float(cleanline[2])])
    return atoms,geometry

def ReadMOL(path2file):
    #function reading a mol and writing to containers
    if not path2file.endswith(".mol"): path2file += ".mol"
    rfile = open(path2file,"r")
    rfilelines = rfile.readlines()
    natms = int(rfilelines[3][0:3])
    atoms = []
    geometry = []
    for iatm in range(natms):
        line = rfilelines[iatm + 4].strip().rstrip()
        cleanline = FORMAT.CleanTheLine(line).split(" ")
        atoms.append(cleanline[3])
        geometry.append([float(cleanline[0]),float(cleanline[1]),float(cleanline[2])])
    return atoms,geometry

def ReadSDFWithTopology(path2file):
    #function reading an sdf and writing connectivity, atomtype, atoms, geometry, and local charges to containers
    if not path2file.endswith(".sdf"): path2file += ".sdf"
    atoms,geometry = ReadSDF(path2file)
    rfile = open(path2file,"r")
    rfilelines = rfile.readlines()
    natms = int(rfilelines[3][0:3])
    nbnds = int(rfilelines[3][3:6])
    #build atomtype: C1, C2, O1, ...
    atomtype = []
    element_count = {}
    for iatm in range(natms):
        el = atoms[iatm]
        if el not in element_count:
            element_count[el] = 0
        element_count[el] += 1
        atomtype.append(el + str(element_count[el]))
    connectivity = []
    for ibnd in range(nbnds):
        line = rfilelines[4 + natms + ibnd].strip().rstrip()
        cleanline = FORMAT.CleanTheLine(line).split(" ")
        atom1 = int(cleanline[0])
        atom2 = int(cleanline[1])
        bond_order = cleanline[2]
        atomA_name = atomtype[atom1 - 1]
        atomB_name = atomtype[atom2 - 1]
        if bond_order == "2": bond_type = "doub"
        elif bond_order == "3": bond_type = "trip"
        elif bond_order == "4": bond_type = "arom"
        else: bond_type = "sing"
        connectivity.append((atomA_name,atomB_name,bond_type))
    #read local charges
    local_charges = [0 for _ in range(natms + 1)]
    for line in rfilelines:
        sline = line.strip().rstrip()
        if sline.startswith("M  CHG"):
            cleanline = FORMAT.CleanTheLine(sline).split(" ")
            nchg = int(cleanline[2])
            for ichg in range(nchg):
                atom_index = int(cleanline[3 + 2*ichg])
                atom_charge = int(cleanline[4 + 2*ichg])
                local_charges[atom_index] = atom_charge
                local_charges[0] += 1
    return connectivity,atomtype,atoms,geometry,local_charges

def ReadMOLWithTopology(path2file):
    #function reading a mol and writing connectivity, atomtype, atoms, geometry, and local charges to containers
    if not path2file.endswith(".mol"): path2file += ".mol"
    atoms,geometry = ReadMOL(path2file)
    rfile = open(path2file,"r")
    rfilelines = rfile.readlines()
    natms = int(rfilelines[3][0:3])
    nbnds = int(rfilelines[3][3:6])
    #build atomtype: C1, C2, O1, ...
    atomtype = []
    element_count = {}
    for iatm in range(natms):
        el = atoms[iatm]
        if el not in element_count:
            element_count[el] = 0
        element_count[el] += 1
        atomtype.append(el + str(element_count[el]))
    connectivity = []
    for ibnd in range(nbnds):
        line = rfilelines[4 + natms + ibnd].strip().rstrip()
        cleanline = FORMAT.CleanTheLine(line).split(" ")
        atom1 = int(cleanline[0])
        atom2 = int(cleanline[1])
        bond_order = cleanline[2]
        atomA_name = atomtype[atom1 - 1]
        atomB_name = atomtype[atom2 - 1]
        if bond_order == "2": bond_type = "doub"
        elif bond_order == "3": bond_type = "trip"
        elif bond_order == "4": bond_type = "arom"
        else: bond_type = "sing"
        connectivity.append((atomA_name,atomB_name,bond_type))
    #read local charges
    local_charges = [0 for _ in range(natms + 1)]
    for line in rfilelines:
        sline = line.strip().rstrip()
        if sline.startswith("M  CHG"):
            cleanline = FORMAT.CleanTheLine(sline).split(" ")
            nchg = int(cleanline[2])
            for ichg in range(nchg):
                atom_index = int(cleanline[3 + 2*ichg])
                atom_charge = int(cleanline[4 + 2*ichg])
                local_charges[atom_index] = atom_charge
                local_charges[0] += 1
    return connectivity,atomtype,atoms,geometry,local_charges

def ObabelXYZ2SDF(ifile,ofile):
    #function that converts an xyz to sdf using obabel
    subprocess.run(["obabel",ifile,"-O",ofile])
    return 0

def CleanSDFfromProtons(sdffile):
    #function that removes protons from an sdf file
    #FMCM: Potential Bug here
    #a potential bug comes from the fact that the indices are not down-shifted after removing protons
    #FMCM: Potential Bug here
    rfile = open(sdffile,"r")
    rfilelines = rfile.readlines()
    rfile.close()
    natoms = int(rfilelines[3][0:3])
    nbonds = int(rfilelines[3][3:6])
    natoms_new = natoms
    nbonds_new = nbonds
    Hindices = []
    lines2remove = []
    #loop over atoms to find the protons
    for idatm in range(natoms):
        line = rfilelines[idatm + 4].strip().rstrip()
        cleanline = FORMAT.CleanTheLine(line)
        linedata = cleanline.split(" ")
        if linedata[3] == "H": 
            Hindices.append(idatm + 1)
            lines2remove.append(idatm + 4)
            natoms_new -= 1
    #loop over bonds to find bonds to protons
    for ibond in range(nbonds):
        line = rfilelines[ibond + 4 + natoms].strip().rstrip()
        cleanline = FORMAT.CleanTheLine(line)
        linedata = cleanline.split(" ")
        atomA = int(linedata[0])
        atomB = int(linedata[1])
        if (atomA in Hindices) or (atomB in Hindices):
            #remove the line
            lines2remove.append(ibond + 4 + natoms)
            nbonds_new -= 1
    #now write the new sdf file
    wfile = open(sdffile,"w")
    wfile.write(rfilelines[0])
    wfile.write(rfilelines[1])
    wfile.write(rfilelines[2])
    atombondheader = rfilelines[3][6:]
    stringatoms = ""
    if natoms_new < 10: stringatoms = "  " + str(natoms_new)
    elif natoms_new < 100: stringatoms = " " + str(natoms_new)
    else: stringatoms = str(natoms_new)
    stringbonds = ""
    if nbonds_new < 10: stringbonds = "  " + str(nbonds_new)
    elif nbonds_new < 100: stringbonds = " " + str(nbonds_new)
    else: stringbonds = str(nbonds_new)
    wfile.write(stringatoms + stringbonds + atombondheader)
    for iline in range(4,len(rfilelines)):
        if iline in lines2remove: continue
        wfile.write(rfilelines[iline])
    wfile.close()

def Spaces_Number(number):
    #normalizing function for SDF files
    spaces = "    "
    if (number > -1000.0) and (number <= -100.0): spaces = " "
    elif (number > -100.0) and (number <= -10.0): spaces = "  "
    elif (number > -10.0) and (number < 0.0): spaces = "   "
    elif (number >= 0.0) and (number < 10.0): spaces = "    "
    elif (number >= 10.0) and (number < 100.0): spaces = "   "
    elif (number >= 100.0) and (number < 1000.0): spaces = "  "
    return spaces

def Spaces_Element(element):
    #function dictating how many spaces are added at the end of the element
    spaces = " "
    if len(element) == 2: spaces = ""
    return spaces

def WriteSDF(atomtype, atoms, geometry, connectivity, local_charges, sdfname, skip_proton = True):
    """Write an SDF (V2000) using mmCIF/PDB-derived connectivity, with robust name canonicalization.
       Changelog:
         - Canonicalize atom names consistently (strip quotes/whitespace, normalize unicode primes, ' -> *)
         - Build bonds_to_write first, then set nbonds = len(bonds_to_write) (counts line cannot drift)
         - Use dict mapping name->index (no fragile list.index lookups)
         - Hydrogen detection based on element list (`atoms`) rather than atom-name heuristic
         - Safer headers (non-empty title line)
       Input/output signature is unchanged."""
    
    #helper functions (local, no external dependencies)
    def canon_atom_name(name: str) -> str:
        """Strip quotes/whitespace and normalize primes: (’, ′, ʼ, etc.) -> ' then ' -> *."""
        if name is None:
            return ""
        s = str(name).strip()
        #remove surrounding single/double quotes repeatedly (handles "\"C1'\"" cases)
        while len(s) >= 2 and ((s[0] == s[-1]) and s[0] in ("'",'"')):
            s = s[1:-1].strip()
        #normalize a bunch of "prime-like" unicode characters to ASCII apostrophe
        #(common in copied text / CIF exports): right quote, modifier letter apostrophe, prime, double prime, etc.
        prime_like = {
            "\u2018", "\u2019",  #left/right single quotation mark
            "\u2032", "\u2033", "\u2034", "\u2035", "\u2036", "\u2037",  # prime family
            "\u02B9", "\u02BA",  #modifier letter prime/double prime
            "\u02BC",            #modifier letter apostrophe
            "\uFF07",            #fullwidth apostrophe
            "\u00B4",            #acute accent (often misused as apostrophe)
            "\u0060"}            #grave accent (often misused as apostrophe)
        for ch in prime_like:
            s = s.replace(ch,"'")
        #now convert any number of ASCII apostrophes to '*' to keep labels simple and parser-safe
        #O''' -> O***, C1' -> C1*
        s = s.replace("'","*")
        return s

    def is_hydrogen_atom_by_element(atom_idx0: int) -> bool:
        """True if element for atom index is hydrogen (or deuterium/tritium)."""
        try:
            el = str(atoms[atom_idx0]).strip().upper()
        except Exception:
            return False
        return el in {"H","D","T"}

    #open file
    if not sdfname.endswith(".sdf"): sdfname += ".sdf"
    wfile = open(sdfname,"w")
    #SDF header: 3 lines (title, program, comment)
    wfile.write(f"{sdfname}\n")
    wfile.write("    In-Pocket-Analysis SDF writer\n")
    wfile.write("\n")
    #build mapping from canonical atom names to 1-based SDF indices
    natoms = len(atomtype)
    canon_names = [canon_atom_name(nm) for nm in atomtype]
    name_to_idx = {}
    duplicates = set()
    for i,nm in enumerate(canon_names,start = 1):
        if nm in name_to_idx:
            duplicates.add(nm)
        else:
            name_to_idx[nm] = i
    #if duplicates exist, we cannot reliably map bonds by name.
    #keep behavior compatible (raise hard error rather than silently writing wrong bonds).
    if duplicates:
        wfile.close()
        raise ValueError(
            "WriteSDF: non-unique atom names after canonicalization (prime->* / quote stripping). "
            f"Duplicates: {sorted(list(duplicates))[:10]}{'...' if len(duplicates) > 10 else ''}"
        )
    #build list of bonds to write (and count from that list)
    bonds_to_write = []
    missing_atoms = set()
    for b in connectivity:
        #expected: (atomA_name, atomB_name, bond_type, ...)
        if len(b) < 3:
            continue
        atomA_raw,atomB_raw,btype = b[0],b[1],b[2]
        atomA = canon_atom_name(atomA_raw)
        atomB = canon_atom_name(atomB_raw)
        if not atomA or not atomB:
            continue
        idxA = name_to_idx.get(atomA,None)
        idxB = name_to_idx.get(atomB,None)
        if idxA is None or idxB is None:
            #store for debugging; skip writing such bonds
            if idxA is None:
                missing_atoms.add(atomA)
            if idxB is None:
                missing_atoms.add(atomB)
            continue
        #hydrogen filtering (robust): skip bonds involving H/D/T atoms based on element list
        if skip_proton and (is_hydrogen_atom_by_element(idxA - 1) or is_hydrogen_atom_by_element(idxB - 1)):
            continue
        #bond order mapping
        border = 1
        if btype == "doub":
            border = 2
        elif btype == "trip":
            border = 3
        elif btype == "arom":
            border = 4
        #else keep single
        bonds_to_write.append((idxA,idxB,border))
    nbonds = len(bonds_to_write)
    #for compatibility, we keep writing an SDF (counts are correct for what we actually wrote),
    #but you still get a clear warning on stdout if desired. Here we raise only if everything breaks.
    if natoms == 0:
        wfile.close()
        raise ValueError("WriteSDF: no atoms provided")
    #counts line (fixed-width)
    #V2000 counts line: 3-digit natoms + 3-digit nbonds (right aligned)
    wfile.write(f"{natoms:>3}{nbonds:>3}  0  0  0  0  0  0  0  0999 V2000\n")
    #atom block
    for idatm in range(len(atoms)):
        spaces_x = Spaces_Number(geometry[idatm][0])
        x_round = "{:.4f}".format(round(geometry[idatm][0], 4))
        spaces_y = Spaces_Number(geometry[idatm][1])
        y_round = "{:.4f}".format(round(geometry[idatm][1], 4))
        spaces_z = Spaces_Number(geometry[idatm][2])
        z_round = "{:.4f}".format(round(geometry[idatm][2], 4))
        spaces_el = Spaces_Element(atoms[idatm])
        wfile.write(spaces_x + x_round + spaces_y + y_round + spaces_z + z_round)
        wfile.write(" " + str(atoms[idatm]) + spaces_el + "  0  0  0  0  0  0  0  0  0  0  0  0\n")
    #bond block
    for idxA,idxB,border in bonds_to_write:
        #indices are 1-based already
        spacesA = "  "
        if 10 <= idxA < 100:
            spacesA = " "
        elif 100 <= idxA < 1000:
            spacesA = ""
        spacesB = "  "
        if 10 <= idxB < 100:
            spacesB = " "
        elif 100 <= idxB < 1000:
            spacesB = ""
        wfile.write(f"{spacesA}{idxA}{spacesB}{idxB}  {border}  0  0  0  0\n")
    if local_charges:
        if local_charges[0] > 0:
            spaces = "  "
            atom_charges = str(local_charges[0])
            if len(atom_charges) == 2: spaces = " "
            elif len(atom_charges) == 3: spaces = ""
            wfile.write("M  CHG" + spaces + atom_charges)
            icount = 0
            for idatm in range(len(atoms)):
                if local_charges[idatm + 1] != 0:
                    #atom index
                    spaces = "   "
                    atom_charges = str(idatm + 1)
                    if len(atom_charges) == 2: spaces = "  "
                    elif len(atom_charges) == 3: spaces = " "
                    elif len(atom_charges) == 4: spaces = ""
                    wfile.write(spaces + atom_charges)
                    #atom charge
                    spaces = "   "
                    atom_charges = str(local_charges[idatm + 1])
                    if len(atom_charges) == 2: spaces = "  "
                    elif len(atom_charges) == 3: spaces = " "
                    elif len(atom_charges) == 4: spaces = ""
                    wfile.write(spaces + atom_charges)
                    icount += 1
                if icount == local_charges[0]: 
                    #done
                    wfile.write("\n")
                    break
    wfile.write("M  END\n$$$$\n\n")
    wfile.close()
    #optional: if there were missing atoms in connectivity, you may want visibility.
    #we avoid printing by default to keep behavior clean; uncomment if desired.
    output_string = ""
    if missing_atoms:
        output_string = f"WriteSDF: connectivity referenced {len(missing_atoms)} atom names not in atomtype. Examples: {sorted(list(missing_atoms))[:10]}"
    return output_string

# ------------------------------------------------------------
# PDB-related functions
# ------------------------------------------------------------
def Write2PDB(pdbfile, name):
    #function writing a pdb file based on a previous, larger one
    if not name.endswith(".pdb"): name += ".pdb"
    wfile = open(name,"w")
    for idatm in range(len(pdbfile)):
        wfile.write(pdbfile[idatm])
    wfile.write("END\n")
    wfile.close()
    return 0

def ReadPDB_base(pdbfile,altloc): 
    #function to read a PDB file and defines the altloc
    geometry = []            #contains the geometry declared in the pdb files
    atoms = []               #contains the list of atoms
    atmtp = []               #list of atom types
    residues = []            #contains the list of residues
    resnumbr = []            #contains the residue numbering
    chainid = []             #chain information
    restp = []               #residue type, i.e., whether this is protein or something else
    inscode = []             #seqid codes
    pdbline = []             #store the PDB lines for each atom
    #process pdbfile
    if not pdbfile.endswith(".pdb"): pdbfile += ".pdb"
    rfile = open(pdbfile,"r")
    rfilelines = rfile.readlines()
    for iline in range(len(rfilelines)):
        line = rfilelines[iline]
        if line.startswith("ATOM") or line.startswith("HETATM"):
            atmname = line[12:16].strip().rstrip()
            residue = line[16:20].strip().rstrip()
            if len(residue) == 4:
                #keep only one occupation, for now only A
                if not residue.startswith(altloc): continue
                residue = residue[1:]
            schainid = line[21].strip().rstrip()
            resnr = int(line[22:26])
            icode = line[26]
            if (icode == " ") or (icode == ""): icode = "0"
            inscode.append(icode)
            xcoord = float(line[30:38])
            ycoord = float(line[38:46])
            zcoord = float(line[46:54])
            element = line[76:78].strip().rstrip()
            geometry.append([xcoord,ycoord,zcoord])
            atoms.append(element)
            residues.append(residue)
            atmtp.append(atmname)
            chainid.append(schainid)
            resnumbr.append(resnr)
            if line.startswith("ATOM"): restp.append("ATOM")
            elif line.startswith("HETATM"): restp.append("HETATM")
            pdbline.append(line)
        if line.startswith("ENDMDL") or line.startswith("END"): break
    return geometry,atoms,atmtp,residues,resnumbr,chainid,restp,inscode,pdbline

def ReadPDBsimple(pdbfile):
    #simple function to read pdbs
    altloc = "A"
    geometry,atoms,_,_,_,_,_,_,_ = ReadPDB_base(pdbfile,altloc)
    return atoms,geometry

def ReadPDB(pdbfile): 
    #simple function to read pdbs
    altloc = "A"
    geometry,atoms,atmtp,residues,resnumbr,chainid,restp,inscode,pdbline = ReadPDB_base(pdbfile,altloc)
    return geometry,atoms,atmtp,residues,resnumbr,chainid,restp,inscode,pdbline

def PDB2XYZ(pdbfile):
    #function converting a PDB file to xyz
    geometry,atoms,_,_,_,_,_,_,_ = ReadPDB(pdbfile)
    Write2XYZ(atoms,geometry,pdbfile.replace(".pdb","H"))
    return pdbfile.replace(".pdb","H.xyz")

def WriteGemmiPDB(geometry,atoms,atmtp,residues,resnumbr,chainid,restp,pdbname):
    #function that builds a Gemmi structure and prints a pdb
    structure = gemmi.Structure()
    structure.add_model(gemmi.Model('1'))
    #get all chains
    chaindata = [[chainid[0],1]]
    ichaincnt = 0
    for idatm in range(1,len(chainid)):
        if chainid[idatm] != chainid[idatm - 1]: 
            chaindata.append([chainid[idatm],1])
            ichaincnt += 1
        else: 
            if resnumbr[idatm] != resnumbr[idatm - 1]: chaindata[ichaincnt][1] += 1
    startatm = 0
    for model in structure:
        for chain in chaindata:
            #create chain
            ichain = gemmi.Chain(chain[0])
            #populate chain
            for ires in range(chain[1]):
                #create residue
                residue = gemmi.Residue()
                #populate residue
                newstartatm = startatm
                for idatm in range(newstartatm,len(atoms)):
                    if (idatm > newstartatm) and ((residues[idatm] != residues[idatm - 1]) or (resnumbr[idatm] != resnumbr[idatm - 1])): break
                    atom = gemmi.Atom()
                    atom.name = atmtp[idatm].replace("#","")
                    atom.element = gemmi.Element(atoms[idatm])
                    atom.pos.x = geometry[idatm][0]
                    atom.pos.y = geometry[idatm][1]
                    atom.pos.z = geometry[idatm][2]
                    residue.add_atom(atom)
                    startatm += 1
                #add residue
                startatm -= 1
                residue.name = residues[startatm]
                residue.seqid.num = resnumbr[startatm]
                residue.het_flag = restp[startatm][0]
                ichain.add_residue(residue)
                startatm += 1
                if (startatm < len(atoms)) and (chainid[startatm] != chain[0]): break
            #add chain
            model.add_chain(ichain)
    if not pdbname.endswith(".pdb"): pdbname += ".pdb"
    structure.write_pdb(pdbname)
    return 0

# ------------------------------------------------------------
# MOL2-related functions
# ------------------------------------------------------------
def ReadMOL2(path2file):
    #function reading a mol2 and writing to containers
    if not path2file.endswith(".mol2"): path2file += ".mol2"
    rfile = open(path2file,"r")
    rfilelines = rfile.readlines()
    atoms = []
    geometry = []
    in_atom_block = False
    for line in rfilelines:
        sline = line.strip().rstrip()
        if sline.startswith("@<TRIPOS>ATOM"):
            in_atom_block = True
            continue
        if sline.startswith("@<TRIPOS>") and in_atom_block:
            break
        if in_atom_block:
            cleanline = FORMAT.CleanTheLine(sline).split(" ")
            #mol2 atom line format:
            #0       1         2 3 4 5
            #atom_id atom_name x y z atom_type ...
            atoms.append(cleanline[5].split(".")[0])
            geometry.append([float(cleanline[2]),float(cleanline[3]),float(cleanline[4])])
    return atoms,geometry

def ReadMOL2WithTopology(path2file):
    #ensure correct extension
    if not path2file.endswith(".mol2"): path2file += ".mol2"
    #read geometry data
    atoms,geometry = ReadMOL2(path2file)
    #now only read topology
    topology = []
    with open(path2file,"r") as rfile:
        rfilelines = rfile.readlines()
    in_bond_block = False
    for line in rfilelines:
        sline = line.strip().rstrip()
        if not sline:
            continue
        #start bond block
        if sline.startswith("@<TRIPOS>BOND"):
            in_bond_block = True
            continue
        #exit bond block when a new section starts
        if sline.startswith("@<TRIPOS>") and in_bond_block:
            break
        if in_bond_block:
            cleanline = FORMAT.CleanTheLine(sline).split(" ")
            #mol2 bond format:
            #0       1        2        3
            #bond_id origin   target   bond_type
            atom1 = int(cleanline[1])
            atom2 = int(cleanline[2])
            bond_type = cleanline[3]
            topology.append((atom1,atom2,bond_type))
    return atoms,geometry,topology

def MOL2toXYZ(mol2file, restp):
    #conversion between mol2 and xyz
    #a few fixes are used for older versions of obabel that do not consider the pH
    geometry = []
    atoms = []
    rfile = open(mol2file,"r")
    rfilelines = rfile.readlines()
    fetch = False
    CDCZ = []
    NE = []
    for line in rfilelines:
        if line.startswith("@<TRIPOS>ATOM"):
            fetch = True
            continue
        elif "@<TRIPOS>UNITY_ATOM_ATTR" in line: break
        elif line.startswith("@<TRIPOS>BOND"): break
        if fetch:
            cleanline = FORMAT.CleanTheLine(line)
            data = cleanline.strip().rstrip().split(" ")
            aux = [float(data[2]),float(data[3]),float(data[4])]
            element = data[5]
            posdot = element.rfind(".")
            if posdot >= 0: element = data[5][:posdot]
            if restp.lower() == "glu" and ((data[1] == "HE1") or (data[1] == "HE2")): continue
            elif restp.lower() == "asp" and ((data[1] == "HD1") or (data[1] == "HD2")): continue
            elif restp.lower() == "arg" and ((data[1] == "CZ") or (data[1] == "CD")): CDCZ.append(aux)
            elif restp.lower() == "arg" and data[1] == "NE": NE = aux
            atoms.append(element)
            geometry.append(aux)
    #add missing protons
    if restp.lower() == "arg": 
        atoms.append("H")
        aux = [-0.5*(CDCZ[0][0] + CDCZ[1][0]),-0.5*(CDCZ[0][1] + CDCZ[1][1]),-0.5*(CDCZ[0][2] + CDCZ[1][2])]
        for idcoord in range(3):
            aux[idcoord] *= 1.5
            aux[idcoord] += 2.5*NE[idcoord]
        geometry.append(aux)
    Write2XYZ(atoms,geometry,mol2file.replace(".mol2",""))
    subprocess.run(["rm",mol2file])

def WriteMOL2(path2file, topology, atoms, geometry, mol_name = "MOL"):
    #function writing topology, atoms, and geometry to a mol2 file
    if not path2file.endswith(".mol2"): path2file += ".mol2"
    if len(atoms) != len(geometry):
        raise ValueError("atoms and geometry must have the same length")
    natoms = len(atoms)
    nbonds = len(topology)
    with open(path2file,"w") as wfile:
        wfile.write("@<TRIPOS>MOLECULE\n")
        wfile.write(f"{mol_name}\n")
        wfile.write(f"{natoms} {nbonds} 0 0 0\n")
        wfile.write("SMALL\n")
        wfile.write("NO_CHARGES\n\n")
        #ATOM block
        wfile.write("@<TRIPOS>ATOM\n")
        for i,(atom,coords) in enumerate(zip(atoms,geometry),start = 1):
            x,y,z = coords
            atom_name = f"{atom}{i}"
            atom_type = atom
            subst_id = 1
            subst_name = "MOL"
            charge = 0.0
            wfile.write(f"{i:7d} {atom_name:<8s} "
                        f"{x:12.6f} {y:12.6f} {z:12.6f} "
                        f"{atom_type:<6s} {subst_id:3d} {subst_name:<8s} {charge:10.6f}\n")
        #BOND block
        wfile.write("@<TRIPOS>BOND\n")
        for i,bond in enumerate(topology,start = 1):
            atom1,atom2,bond_type = bond
            wfile.write(f"{i:6d} {atom1:4d} {atom2:4d} {str(bond_type)}\n")

# ------------------------------------------------------------
# CIF-related functions
# ------------------------------------------------------------
_ALTLOC_META_PREFIX = "__IPA_MMCIF_META__"

def _normalise_altloc_label(altloc):
    """Normalise altloc labels to uppercase single-token strings; blank-like values -> ''."""
    if altloc is None:
        return ""
    sval = str(altloc).strip()
    if sval in ("",".","?","\x00"):
        return ""
    return sval.upper()

def _safe_float(value, default = 0.0):
    try:
        return float(value)
    except Exception:
        return float(default)

def _encode_mmcif_atom_meta(altloc = "", occupancy = 1.0, selection = None):
    """Store mmCIF-only atom metadata inside pdblines while preserving the public API:
       pdblines remains a list[str]."""
    payload = {"src": "mmCIF",
               "altloc": _normalise_altloc_label(altloc),
               "occ": _safe_float(occupancy,1.0)}
    if isinstance(selection,dict) and selection:
        payload["selection"] = selection
    return _ALTLOC_META_PREFIX + json.dumps(payload,separators = (",",":"),sort_keys = True)

def _decode_mmcif_atom_meta(raw):
    """Decode mmCIF atom metadata stored in pdblines; return None if unavailable."""
    if not isinstance(raw,str):
        return None
    if not raw.startswith(_ALTLOC_META_PREFIX):
        return None
    try:
        payload = json.loads(raw[len(_ALTLOC_META_PREFIX):])
        if isinstance(payload,dict):
            return payload
    except Exception:
        pass
    return None

def _format_altloc_site_label(chain_name, residue_name, residue_number, seqidcode, atom_name):
    """Human-readable atom-site label for warnings/logging."""
    chain_name = str(chain_name).strip()
    residue_name = str(residue_name).strip()
    residue_number = str(residue_number).strip()
    atom_name = str(atom_name).strip()
    seqidcode = str(seqidcode).strip()
    if seqidcode in (""," ","0"):
        return f"{residue_name} {chain_name}{residue_number} atom {atom_name}"
    return f"{residue_name} {chain_name}{residue_number}{seqidcode} atom {atom_name}"

def GetAltLocSelectionWarnings(pdblines):
    """Extract unique altloc fallback warnings from mmCIF metadata embedded in pdblines.
       Returns a list of pre-formatted warning strings suitable for logger.warning(...)."""
    warnings = []
    seen = set()
    for raw in pdblines:
        meta = _decode_mmcif_atom_meta(raw)
        if not meta:
            continue
        selection = meta.get("selection",{})
        if not isinstance(selection,dict):
            continue
        if selection.get("status") != "fallback_highest_occupancy":
            continue
        requested = _normalise_altloc_label(selection.get("requested_altloc","")) or "<?>"
        selected = _normalise_altloc_label(selection.get("selected_altloc","")) or "."
        available = selection.get("available_altlocs",[])
        if not isinstance(available,list):
            available = []
        available = [lab for lab in [_normalise_altloc_label(x) for x in available] if lab]
        site = str(selection.get("site","<?>")).strip() or "<?>"
        selected_occ = selection.get("selected_occupancy",None)
        msg = (f"Requested altloc '{requested}' not present for {site}; "
               f"available altloc(s): {', '.join(available) if available else '<none>'}. "
               f"Falling back to highest occupancy")
        if selected != ".":
            msg += f" altloc '{selected}'"
        if selected_occ is not None:
            try:
                msg += f" (occupancy = {float(selected_occ):.3f})"
            except Exception:
                pass
        msg += "."
        if msg not in seen:
            seen.add(msg)
            warnings.append(msg)
    return warnings

def ReadCIFStructure(mmciffile):
    #wrapper around Gemmi to read structure from mmCIF
    structure = gemmi.read_structure(mmciffile)
    atoms = []
    atmtp = []
    residues = []
    resnumbr = []
    chainid = []
    geometry = []
    restp = []
    seqidcodes = []
    pdblines = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atoms.append(atom.element.name)
                    atmtp.append(atom.name)
                    residues.append(residue.name)
                    resnumbr.append(residue.seqid.num)
                    chainid.append(chain.name)
                    geometry.append([atom.pos.x,atom.pos.y,atom.pos.z])
                    residue_type = (residue.het_flag == "A")*("ATOM") + (residue.het_flag == "H")*"HETATM"
                    restp.append(residue_type)
                    icode = residue.seqid.icode
                    if (icode == "") or (icode == " "): icode = "0"
                    seqidcodes.append(icode)
                    #pdblines.append("mmCIF")
                    #preserve true mmCIF altloc/occupancy metadata while keeping
                    #the public output type unchanged (pdblines is still list[str]).
                    pdblines.append(_encode_mmcif_atom_meta(altloc = getattr(atom,"altloc",""),occupancy = getattr(atom,"occ",1.0)))
                    #print(atom.element.name,atom.name,residue.seqid.num,residue.name,chain.name,type(atom.pos.x),type(atom.pos.y),type(atom.pos.z),residue_type)
        break       #take only the first model
    return geometry,atoms,atmtp,residues,resnumbr,chainid,restp,seqidcodes,pdblines

def mmCIF2PDB(mmciffile):
    #function that takes an mmCIF file and converts it to PDB with the full connect cards
    #doc = gemmi.cif.read_file(mmciffile)
    structure = gemmi.read_structure(mmciffile)
    structure.write_pdb(mmciffile.replace("cif","pdb"),gemmi.PdbWriteOptions(preserve_serial = True,conect_records = True))

# ------------------------------------------------------------
# Geometry-reading wrapper functions
# ------------------------------------------------------------
def ReadGeometry(path2file):
    #general wrapper
    atoms = []
    geometry = []
    if path2file.endswith(".xyz"): atoms,geometry = ReadXYZ(path2file)
    elif path2file.endswith(".sdf"): atoms,geometry = ReadSDF(path2file)
    elif path2file.endswith(".mol"): atoms,geometry = ReadMOL(path2file)
    elif path2file.endswith(".mol2"): atoms,geometry = ReadMOL2(path2file)
    elif path2file.endswith(".pdb"): atoms,geometry = ReadPDBsimple(path2file.replace(".pdb",""))
    return atoms, geometry

# ------------------------------------------------------------
# Protonation functions
# ------------------------------------------------------------
def ObabelProtonation(path2file, pH):
    #function protonating the ligand using obabel
    fileextension = ".pdb"*path2file.endswith(".pdb") + ".xyz"*path2file.endswith(".xyz") + ".sdf"*path2file.endswith(".sdf")
    #modify here to use pybel
    #->
    intermedmol2 = path2file.replace(fileextension,".mol2")
    finalmol2 = intermedmol2.replace(".mol2","H.mol2")
    subprocess.run(["obabel",path2file,"-O",intermedmol2])
    subprocess.run(["obabel",intermedmol2,"-O",finalmol2,"-p",str(pH)])
    subprocess.run(["rm",intermedmol2])
    #<-
    #get the total charge of the molecule
    totalcharge = 0.0
    rfile = open(finalmol2,"r")
    rfilelines = rfile.readlines()
    fetch = False
    for line in rfilelines:
        if "@<TRIPOS>ATOM" in line: 
            fetch = True
            continue
        elif "@<TRIPOS>UNITY_ATOM_ATTR" in line: break
        elif "@<TRIPOS>BOND" in line: break
        if fetch:
            cleanline = FORMAT.CleanTheLine(line)
            data = cleanline.strip().rstrip().split(" ")
            totalcharge += float(data[8])
    return round(totalcharge)

def AddCapProton(cap_atoms, cap_geom, xyzfile):
    #function that adds proton caps in an xyz file that is closest to an atom to cap
    #find the C atom
    atmC = -1
    for iatm in range(len(cap_atoms)):
        if cap_atoms[iatm] == "C": 
            atmC = iatm
            break
    if atmC > -1:
        mindistH = 10000000000000000000000000000000000000000000000000000000000000000.0
        imindistH = -1
        mindistX = 10000000000000000000000000000000000000000000000000000000000000000.0
        imindistX = -1
        res_atoms,res_geom = ReadGeometry(xyzfile)
        for iatm in range(len(res_atoms)):
            raux = 0.0
            for idcoord in range(3):
                intermval = cap_geom[atmC][idcoord] - res_geom[iatm][idcoord]
                raux += intermval*intermval
            if (res_atoms[iatm] == "H"):
                if (raux < mindistH):
                    mindistH = raux
                    imindistH = iatm
            else:
                if (raux < mindistX) and ((raux > 0.1)):
                    overlap = False
                    for ibtm in range(len(cap_atoms)):
                        if (cap_geom[ibtm][0] == res_geom[iatm][0]) and (cap_geom[ibtm][1] == res_geom[iatm][1]) and (cap_geom[ibtm][2] == res_geom[iatm][2]):
                            overlap = True
                            break
                    if not overlap:
                        mindistX = raux
                        imindistX = iatm
        if (mindistH < 2.25) and (res_atoms[imindistH] == "H"):
            cap_atoms.append("H")
            cap_geom.append(res_geom[imindistH])
        if (mindistX < 3.24):
            cap_atoms.append("H")
            cap_geom.append(res_geom[imindistX])
    return 0

# ------------------------------------------------------------
# Geometry-related functions
# ------------------------------------------------------------
def ReCenterAndWrite(reference, modify, output):
    #function that takes a reference function, gets its center of mass; gets another, similar structure to recenter, and writes to file
    ref_atoms,ref_geometry = ReadGeometry(reference)
    mod_atoms,mod_geometry = ReadGeometry(modify)
    ref_com = [0.0,0.0,0.0]
    mod_com = [0.0,0.0,0.0]
    nheteroatm = 0
    #assume that the 2 structures have the same order of atoms
    for iatm in range(len(ref_atoms)):
        if mod_atoms[iatm] == "H": continue
        nheteroatm += 1
        for idcoord in range(3):
            ref_com[idcoord] += ref_geometry[iatm][idcoord]
            mod_com[idcoord] += mod_geometry[iatm][idcoord]
    for idcoord in range(3):
        ref_com[idcoord] /= nheteroatm
        mod_com[idcoord] /= nheteroatm
    for iatm in range(len(ref_atoms)):
        for idcoord in range(3):
            mod_geometry[iatm][idcoord] += ref_com[idcoord] - mod_com[idcoord]
    Write2XYZ(mod_atoms,mod_geometry,output)
    return 0

def GeneralisedReCenterAndWrite(reference, modify, output, skip):
    #function that takes a reference function, gets its center of mass; gets another, similar structure to recenter, and writes to file
    #the difference with respect to the previous one is that not all atoms enter in the calculation of the center of mass
    ref_atoms,ref_geometry = ReadGeometry(reference)
    mod_atoms,mod_geometry = ReadGeometry(modify)
    ref_com = [0.0,0.0,0.0]
    mod_com = [0.0,0.0,0.0]
    nheteroatm = 0
    #assume that the 2 structures have the same order of atoms
    for iatm in range(len(ref_atoms)):
        if mod_atoms[iatm] == "H": continue
        elif iatm + 1 in skip: continue
        nheteroatm += 1
        for idcoord in range(3):
            ref_com[idcoord] += ref_geometry[iatm][idcoord]
            mod_com[idcoord] += mod_geometry[iatm][idcoord]
    for idcoord in range(3):
        ref_com[idcoord] /= nheteroatm
        mod_com[idcoord] /= nheteroatm
    for iatm in range(len(ref_atoms)):
        for idcoord in range(3):
            mod_geometry[iatm][idcoord] += ref_com[idcoord] - mod_com[idcoord]
    Write2XYZ(mod_atoms,mod_geometry,output)
    return 0

def _SelectOneAltLocOnly_legacy(geometry, atoms, atmtp, residues, resnumbr, chainid, restp, seqidcodes, pdblines, altloc):
    #legacy heuristic kept as fallback for callers that do not provide mmCIF metadata
    new_geometry = []
    new_atoms = []
    new_atmtp = []
    new_residues = []
    new_resnumbr = []
    new_chainid = []
    new_restp = []
    new_seqidcodes = []
    new_pdblines = []
    skip_list = [False]*len(atmtp)
    for idatm in range(len(atmtp)):
        if skip_list[idatm]: continue
        if idatm != len(atmtp) - 1:
            if (atmtp[idatm] == atmtp[idatm + 1]) and (resnumbr[idatm] == resnumbr[idatm + 1]):
                #altloc identified
                altlocs = [idatm,idatm + 1]
                #get all altlocs
                for idbtm in range(idatm + 2,len(atmtp)):
                    if atmtp[idatm] == atmtp[idbtm]: altlocs.append(idbtm)
                    else: break
                idxsel = 0
                if altloc == "A": idxsel = 0
                elif altloc == "B": idxsel = 1
                elif altloc == "C": idxsel = 2
                elif altloc == "D": idxsel = 3
                elif altloc == "E": idxsel = 4
                new_geometry.append(geometry[altlocs[idxsel]])
                new_atoms.append(atoms[altlocs[idxsel]])
                new_atmtp.append(atmtp[altlocs[idxsel]])
                new_residues.append(residues[altlocs[idxsel]])
                new_resnumbr.append(resnumbr[altlocs[idxsel]])
                new_chainid.append(chainid[altlocs[idxsel]])
                new_restp.append(restp[altlocs[idxsel]])
                new_seqidcodes.append(seqidcodes[altlocs[idxsel]])
                new_pdblines.append(pdblines[altlocs[idxsel]])
                #atoms to skip
                for ialt in range(len(altlocs)):
                    skip_list[altlocs[ialt]] = True
            else: 
                new_geometry.append(geometry[idatm])
                new_atoms.append(atoms[idatm])
                new_atmtp.append(atmtp[idatm])
                new_residues.append(residues[idatm])
                new_resnumbr.append(resnumbr[idatm])
                new_chainid.append(chainid[idatm])
                new_restp.append(restp[idatm])
                new_seqidcodes.append(seqidcodes[idatm])
                new_pdblines.append(pdblines[idatm])
        else: 
            new_geometry.append(geometry[idatm])
            new_atoms.append(atoms[idatm])
            new_atmtp.append(atmtp[idatm])
            new_residues.append(residues[idatm])
            new_resnumbr.append(resnumbr[idatm])
            new_chainid.append(chainid[idatm])
            new_restp.append(restp[idatm])
            new_seqidcodes.append(seqidcodes[idatm])
            new_pdblines.append(pdblines[idatm])
    return new_geometry,new_atoms,new_atmtp,new_residues,new_resnumbr,new_chainid,new_restp,new_seqidcodes,new_pdblines

def SelectOneAltLocOnly(geometry,atoms,atmtp,residues,resnumbr,chainid,restp,seqidcodes,pdblines,altloc):
    """Select exactly one atom per atom-site using true mmCIF altloc metadata when available.
       Public API intentionally unchanged:
         input  -> 9 parallel arrays + altloc
         output -> same 9 parallel arrays
       Behaviour:
         - atom-sites are grouped by (chain, resnum, icode, residue name, atom name)
         - if requested altloc is present for a site, pick that altloc
         - if no explicit altloc exists for a site, keep the single/common site
         - if explicit altlocs exist but requested altloc is missing, fall back to the
           highest-occupancy candidate and annotate pdblines so callers can warn"""
    requested_altloc = _normalise_altloc_label(altloc)
    meta_by_idx = [_decode_mmcif_atom_meta(raw) for raw in pdblines]
    #no mmCIF metadata available -> preserve legacy behaviour exactly
    if not any(meta is not None for meta in meta_by_idx):
        return _SelectOneAltLocOnly_legacy(geometry,atoms,atmtp,residues,resnumbr,chainid,restp,seqidcodes,pdblines,altloc)
    new_geometry = []
    new_atoms = []
    new_atmtp = []
    new_residues = []
    new_resnumbr = []
    new_chainid = []
    new_restp = []
    new_seqidcodes = []
    new_pdblines = []

    def _append_selected(idxsel,selection_info = None):
        new_geometry.append(geometry[idxsel])
        new_atoms.append(atoms[idxsel])
        new_atmtp.append(atmtp[idxsel])
        new_residues.append(residues[idxsel])
        new_resnumbr.append(resnumbr[idxsel])
        new_chainid.append(chainid[idxsel])
        new_restp.append(restp[idxsel])
        new_seqidcodes.append(seqidcodes[idxsel])
        meta = dict(meta_by_idx[idxsel] or {})
        if selection_info:
            meta["selection"] = selection_info
        if meta:
            new_pdblines.append(_encode_mmcif_atom_meta(altloc = meta.get("altloc",""),occupancy = meta.get("occ",1.0),selection = meta.get("selection",None)))
        else:
            new_pdblines.append(pdblines[idxsel])

    #group by true atom-site identity
    atom_groups = {}
    for idx in range(len(atmtp)):
        key = (str(chainid[idx]).strip(),
               int(resnumbr[idx]),
               str(seqidcodes[idx]).strip() if str(seqidcodes[idx]).strip() else "0",
               str(residues[idx]).strip(),
               str(atmtp[idx]).strip())
        if key not in atom_groups:
            atom_groups[key] = []
        atom_groups[key].append(idx)
    for key,idx_list in atom_groups.items():
        chain_name,residue_number,icode,residue_name,atom_name = key
        site_label = _format_altloc_site_label(chain_name,residue_name,residue_number,icode,atom_name)
        candidates = []
        for idx in idx_list:
            meta = meta_by_idx[idx] or {}
            candidates.append((idx,_normalise_altloc_label(meta.get("altloc","")),_safe_float(meta.get("occ",0.0),0.0)))
        
        #deterministic ranking: highest occupancy, then first appearance
        def _rank_key(entry):
            idx0,_,occ = entry
            return (-occ,idx0)

        explicit = [entry for entry in candidates if entry[1] != ""]
        requested = [entry for entry in explicit if entry[1] == requested_altloc]
        selection_info = None
        if requested:
            idxsel,selected_altloc,selected_occ = sorted(requested,key = _rank_key)[0]
        elif explicit:
            #requested altloc absent -> fall back to highest occupancy
            idxsel,selected_altloc,selected_occ = sorted(candidates,key = _rank_key)[0]
            selection_info = {"status": "fallback_highest_occupancy",
                              "requested_altloc": requested_altloc,
                              "selected_altloc": selected_altloc,
                              "selected_occupancy": selected_occ,
                              "available_altlocs": sorted({entry[1] for entry in explicit}),
                              "site": site_label}
        else:
            #no explicit altloc at this site -> keep the common/single site
            idxsel,selected_altloc,selected_occ = sorted(candidates,key = _rank_key)[0]
        _append_selected(idxsel,selection_info = selection_info)
    return new_geometry,new_atoms,new_atmtp,new_residues,new_resnumbr,new_chainid,new_restp,new_seqidcodes,new_pdblines

#ipa_residues -> self.residue_states_from_IPA
#gpath -> self.path_to_IPA_geometries
def RebuildPocketFromResidues(ipa_residues, gpath):
    """Function that rebuilts a geometry and atom list from a list of selected residues, including capping, protons, and overlapping atoms"""
    pocket_atoms = []
    pocket_geometry = []
    pocket_charge = 0
    status = 1
    for ires in range(len(ipa_residues)):
        residue = ipa_residues[ires]
        if (ires > 0):
            prev_residue = ipa_residues[ires - 1]
            if (residue[0] == prev_residue[0]) and (residue[1] == prev_residue[1]) and (residue[2] == prev_residue[2]):
                #two residues are repeated, this means ambiguity in protonation states
                status = 0
                break
        residue_file = gpath + residue[4]
        residue_atoms,residue_geometry = ReadXYZ(residue_file)
        pocket_charge += (residue[3] == "+") - (residue[3] == "-")
        pocket_charge += (residue[3] == "+1") - (residue[3] == "-1")
        pocket_charge += 2*(residue[3] == "++") - 2*(residue[3] == "--")
        pocket_charge += 2*(residue[3] == "+2") - 2*(residue[3] == "-2")
        for iatm in range(len(residue_geometry)):
            atom_A = residue_geometry[iatm]
            clash = False
            if len(pocket_geometry) > 0:
                for ibtm in range(len(pocket_geometry)):
                    atom_B = pocket_geometry[ibtm]
                    distance = [atom_A[0] - atom_B[0],atom_A[1] - atom_B[1],atom_A[2] - atom_B[2]]
                    R2 = distance[0]*distance[0] + distance[1]*distance[1] + distance[2]*distance[2]
                    if R2 < 0.5: 
                        clash = True
                        if pocket_atoms[ibtm] == "H":
                            #replace proton
                            pocket_atoms[ibtm] = residue_atoms[iatm]
                            pocket_geometry[ibtm] = atom_A
                    if clash: break
            if not clash: 
                pocket_geometry.append(atom_A)
                pocket_atoms.append(residue_atoms[iatm])
    return status,pocket_atoms,pocket_geometry,pocket_charge
