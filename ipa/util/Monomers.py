# ------------------------------------------------------------
# Monomer-library-related functions
# ------------------------------------------------------------

from util.monomer_library_ccds import _ALL_CCDs_MONOMER_LIBRARY
import util.Formatting as FORMAT
import gemmi

# ------------------------------------------------------------
# General, monomer-related functions
# ------------------------------------------------------------
def MonomerExistsInLibrary(residue):
    """Return True if `residue` is found in the monomer library list."""
    return str(residue).strip().upper() in _ALL_CCDs_MONOMER_LIBRARY

def GetMonomerPath(residue):
    """Get the path to the monomer file."""
    #if we're here, then the monomer file exists
    base_folder = "util/monomers_library/"
    monomer_order_folder = residue[0].lower() + "/"
    monomer_file = residue.upper() + ".cif"
    return base_folder + monomer_order_folder + monomer_file

# ------------------------------------------------------------
# Connectivity extractors
# ------------------------------------------------------------
def ExtractCCDconnectivityGEMMI(cif_path: str):
    """Extract [atom1, atom2, bo, arom] from a CCP4/CCD-style monomer CIF.
       Designed for MonomerLibrary/monomers layout where files may contain multiple CIF blocks."""
    doc = gemmi.cif.read_file(cif_path)
    bond_blocks = []
    for b in doc:
        cat = b.find_mmcif_category("_chem_comp_bond")
        if cat is not None and len(cat) > 0:
            bond_blocks.append(b)
    if len(bond_blocks) == 0:
        atom_count = 0
        for b in doc:
            atom_cat = b.find_mmcif_category("_chem_comp_atom")
            if atom_cat is not None:
                atom_count += len(atom_cat)
        if atom_count == 1:
            return []
        else:
            raise ValueError(f"{cif_path}: no _chem_comp_bond table found")
    if len(bond_blocks) > 1:
        raise ValueError(f"{cif_path}: multiple blocks contain _chem_comp_bond ({len(bond_blocks)})")
    bonds = bond_blocks[0].find_mmcif_category("_chem_comp_bond")
    type_map = {
        "SINGLE": "sing", "SING": "sing",
        "DOUBLE": "doub", "DOUB": "doub",
        "TRIPLE": "trip", "TRIP": "trip",
        "AMIDE": "sing",
        "AROMATIC": "arom",
        "DELOCALIZED": "arom"
    }
    cols = list(getattr(bonds,"tags",[]))
    if not cols and hasattr(bonds,"get_tags"):
        cols = list(bonds.get_tags())
    if not cols:
        raise ValueError(f"{cif_path}: cannot determine _chem_comp_bond column names")
    def col_idx_any(*wanted_names: str) -> int:
        #wanted_names are full tag names like '_chem_comp_bond.type'
        for name in wanted_names:
            if name in cols:
                return cols.index(name)
        raise ValueError(f"{cif_path}: Column name not found: {wanted_names} (available: {cols})")
    i_a1 = col_idx_any("_chem_comp_bond.atom_id_1")
    i_a2 = col_idx_any("_chem_comp_bond.atom_id_2")
    #prefer 'type', else use 'value_order'
    i_type = col_idx_any("_chem_comp_bond.type","_chem_comp_bond.value_order")
    #prefer 'aromatic', else use 'pdbx_aromatic_flag'
    i_arom = col_idx_any("_chem_comp_bond.aromatic","_chem_comp_bond.pdbx_aromatic_flag")
    out = []
    for row in bonds:
        a1 = row[i_a1]
        a2 = row[i_a2]
        raw_type = row[i_type].strip().upper()
        bo = type_map.get(raw_type,raw_type.lower() if raw_type else "sing")
        arom_val = row[i_arom]
        if arom_val is None:
            arom = False
        else:
            arom = arom_val.strip().lower() in {"y","yes","1","true","t"}
        out.append([a1,a2,bo,arom])
    return out

def ExtractLigandConnectivity(mmciffile, chainID, resnumb):
    #function that extracts ligand connectivity data from the mmCIF
    #read the file
    rfile = open(mmciffile,"r")
    rfilelines = rfile.readlines()
    rfile.close()
    sresnumb = str(resnumb)
    #just to facilitate searches
    idx_chem_comp_bond_start = -1
    #look for ligand code
    found_scheme = False
    counter = 0
    mon_id = -1
    pdb_mon_id = -1
    auth_mon_id = -1
    pdb_strand_id = -1
    pdb_seq_num = -1
    auth_seq_num = -1
    smon_id = ""
    spdb_mon_id = ""
    sauth_mon_id = ""
    for idl in range(len(rfilelines)):
        line = rfilelines[idl]
        if line.startswith("_chem_comp_bond."): idx_chem_comp_bond_start = idl
        if line.startswith("_pdbx_nonpoly_scheme."):
            found_scheme = True
            if ".mon_id" in line: 
                mon_id = counter
                if line.replace("_pdbx_nonpoly_scheme.mon_id","").replace(" ","") + "blahblah" != "blahblah":
                    smon_id = line.replace("_pdbx_nonpoly_scheme.mon_id","").replace(" ","").replace("\n","")
            elif ".pdb_mon_id" in line: pdb_mon_id = counter
            elif ".auth_mon_id" in line: auth_mon_id = counter
            elif ".pdb_strand_id" in line: pdb_strand_id = counter
            elif ".pdb_seq_num" in line: pdb_seq_num = counter
            elif ".auth_seq_num" in line: auth_seq_num = counter
            counter += 1
            continue
        if found_scheme:
            if "#" in line: 
                if idx_chem_comp_bond_start < 0: idx_chem_comp_bond_start = idl
                break
            cleanline = FORMAT.CleanTheLine(line).split(" ")
            if (chainID == cleanline[pdb_strand_id]) and ((cleanline[auth_seq_num] == sresnumb) or (cleanline[pdb_seq_num] == sresnumb)):
                smon_id = cleanline[mon_id]
                spdb_mon_id = cleanline[pdb_mon_id]
                sauth_mon_id = cleanline[auth_mon_id]
                break
    #search for chem_comp_bond
    found_ccb = False
    counter = 0
    comp_id = -1
    atom_id_1 = -1
    atom_id_2 = -1
    value_order = -1
    pdbx_aromatic_flag = -1
    connectivity = []
    for idl in range(idx_chem_comp_bond_start,len(rfilelines)):
        line = rfilelines[idl]
        if line.startswith("_chem_comp_bond."):
            found_ccb = True
            if ".comp_id" in line: comp_id = counter
            elif ".atom_id_1" in line: atom_id_1 = counter
            elif ".atom_id_2" in line: atom_id_2 = counter
            elif ".value_order" in line: value_order = counter
            elif ".pdbx_aromatic_flag" in line: pdbx_aromatic_flag = counter
            counter += 1
            continue
        if found_ccb:
            if "#" in line: break
            cleanline = FORMAT.CleanTheLine(line).split(" ")
            if cleanline[comp_id] == smon_id:
                #good to store
                atomA = cleanline[atom_id_1]
                atomB = cleanline[atom_id_2]
                BO = cleanline[value_order]
                ar = cleanline[pdbx_aromatic_flag]
                aromaticity = False
                if ar.lower() == "y": aromaticity = True
                connectivity.append([atomA,atomB,BO,aromaticity])
    return connectivity

def GeneralisedExtractCCDconnectivityGEMMI(cif_path: str, monomer_code: str):
    """Extract [atom1, atom2, bo, arom] for a specific monomer_code from a CIF that
       may define multiple monomers (multiple blocks and/or multiple _chem_comp_bond tables).
       Supports both schemas seen in monomer libraries:
         - _chem_comp_bond.type + _chem_comp_bond.aromatic
         - _chem_comp_bond.value_order + _chem_comp_bond.pdbx_aromatic_flag"""
    mon = str(monomer_code).strip().upper()
    doc = gemmi.cif.read_file(str(cif_path))
    type_map = {
        #long forms
        "SINGLE": "sing", "DOUBLE": "doub", "TRIPLE": "trip",
        "AMIDE": "sing",
        "AROMATIC": "arom", "DELOCALIZED": "arom",
        #short forms (value_order)
        "SING": "sing", "DOUB": "doub", "TRIP": "trip"
    }

    def truthy(val: str) -> bool:
        return val.strip().lower() in {"y","yes","1","true","t"}
    
    def pick_cols(cols: list[str]):
        #resolve required atom columns
        
        def idx(name: str) -> int:
            if name in cols:
                return cols.index(name)
            #try suffix match, e.g., '_chem_comp_bond.atom_id_1'
            suf = "." + name
            for i,c in enumerate(cols):
                if c.endswith(suf):
                    return i
            raise ValueError(f"Column not found: {name} (available: {cols})")
        
        i_a1 = idx("atom_id_1")
        i_a2 = idx("atom_id_2")
        #bond order column: prefer 'type', else 'value_order'
        i_ord = None
        for candidate in ("_chem_comp_bond.type","type","_chem_comp_bond.value_order","value_order"):
            if candidate in cols:
                i_ord = cols.index(candidate)
                break
        if i_ord is None:
            #suffix fallback
            for i,c in enumerate(cols):
                if c.endswith(".type") or c.endswith(".value_order"):
                    i_ord = i
                    break
        if i_ord is None:
            raise ValueError(f"Bond order column not found (type/value_order). Available: {cols}")
        #aromatic flag column: prefer 'aromatic', else 'pdbx_aromatic_flag'
        i_arom = None
        for candidate in ("_chem_comp_bond.aromatic","aromatic","_chem_comp_bond.pdbx_aromatic_flag","pdbx_aromatic_flag"):
            if candidate in cols:
                i_arom = cols.index(candidate)
                break
        if i_arom is None:
            #suffix fallback
            for i,c in enumerate(cols):
                if c.endswith(".aromatic") or c.endswith(".pdbx_aromatic_flag"):
                    i_arom = i
                    break
        if i_arom is None:
            raise ValueError(f"Aromatic flag column not found (aromatic/pdbx_aromatic_flag). Available: {cols}")
        #comp_id column: required to select the monomer within a multi-component table
        i_comp = None
        for candidate in ("_chem_comp_bond.comp_id","comp_id"):
            if candidate in cols:
                i_comp = cols.index(candidate)
                break
        if i_comp is None:
            for i,c in enumerate(cols):
                if c.endswith(".comp_id"):
                    i_comp = i
                    break
        if i_comp is None:
            raise ValueError(f"comp_id column not found. Available: {cols}")
        return i_comp,i_a1,i_a2,i_ord,i_arom
    
    out = []
    found_any_table = False
    #scan all blocks; pull rows from any _chem_comp_bond tables whose comp_id matches monomer_code
    for b in doc:
        bonds = b.find_mmcif_category("_chem_comp_bond")
        if bonds is None or len(bonds) == 0:
            continue
        found_any_table = True
        #column names in your gemmi build are available as 'tags'
        cols = list(getattr(bonds,"tags",[]))
        if not cols and hasattr(bonds,"get_tags"):
            cols = list(bonds.get_tags())
        if not cols:
            raise ValueError(f"{cif_path}: cannot determine _chem_comp_bond column names")
        i_comp,i_a1,i_a2,i_ord,i_arom = pick_cols(cols)
        for row in bonds:
            if row[i_comp].strip().upper() != mon:
                continue
            a1 = row[i_a1]
            a2 = row[i_a2]
            raw_ord = row[i_ord].strip().upper()
            bo = type_map.get(raw_ord,raw_ord.lower() if raw_ord else "sing")
            arom_val = row[i_arom]
            arom = False if arom_val is None else truthy(arom_val)
            out.append([a1,a2,bo,arom])
    if not found_any_table:
        #if nothing was found, then lets search the monomer library
        if MonomerExistsInLibrary(monomer_code):
            path2monomer = GetMonomerPath(monomer_code)
            out = ExtractCCDconnectivityGEMMI(path2monomer)
        else: 
            raise ValueError(f"{cif_path}: no _chem_comp_bond table found")
    if not out:
        raise ValueError(f"{cif_path}: no bonds found for monomer code {mon}")
    return out
