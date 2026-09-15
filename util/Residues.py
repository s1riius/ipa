# ------------------------------------------------------------
# Residue-related utils
# ------------------------------------------------------------

from dataclasses import dataclass,field
from typing import List,Optional,Tuple
import numpy as np
import os
import util.BaseAminoAcids as AA
import util.Geometry as GEOM
import util.Monomers as MON
import util.VSEPR as VSEPR
import util.Ions as ION

@dataclass
class ResidueData:
    """Stores data for a single residue."""
    geometry: List[List[float]]
    atoms: List[str]
    atom_types: List[str]
    name: str
    number: int
    chain_id: str
    residue_type: str
    pdb_lines: List[str]
    cap_geometry: List[List[float]] = field(default_factory = list)
    cap_atoms: List[str] = field(default_factory = list)
    charge: int = 0
    terminal: str = "0"            #add terminal info
    residue_id: str = ""
    protonation_state: str = "0"
    seqidcode: str = "0"
    processed_atoms: List[str] = field(default_factory = list)
    processed_geometry: List[List[float]] = field(default_factory = list)
    protonated_file: Optional[str] = None
    ip_file: Optional[str] = None

    def __post_init__(self):
        """Set residue_id after initialization if not provided."""
        if not self.residue_id:
            self.residue_id = f"{self.name}_{self.chain_id}{self.number}"

    def is_standard_amino_acid(self) -> bool:
        """Check if this is a standard amino acid."""
        return AA.IsBaseAminoAcid(self.name[:3])

    def is_water(self) -> bool:
        """Check if this is a water molecule."""
        return self.name.upper() in ["HOH","SOL","WAT","H2O"]

def restructure_residues(geometry, atoms, atmtp, residues, resnumbr, chainid, restp, seqidcodes, pdblines, terminal) -> Tuple[List[ResidueData], List[Tuple[int, int, str]]]:
    """Transform a full list of residues into an organised dictionary."""
    residues_dict = {}
    non_protein_residue_info = {}
    index = 0
    for i,(geom,atom,atp,res,resnum,chain,rt,iseq,pdb,term) in enumerate(zip(geometry,atoms,atmtp,residues,resnumbr,chainid,restp,seqidcodes,pdblines,terminal)):
        if i > 0:
            if (resnumbr[i] != resnumbr[i - 1]) or (seqidcodes[i] != seqidcodes[i - 1]):
                index += 1
        res_key = index
        if res_key not in residues_dict:
            residues_dict[res_key] = ResidueData([],[],[],res,resnum,chain,rt,[],terminal = term,seqidcode = iseq)
            if rt == "HETATM": non_protein_residue_info[res_key] = [index,resnum,chain]
            elif rt == "ATOM": non_protein_residue_info[res_key] = [-1,-1000,"-1"]
        residues_dict[res_key].geometry.append(geom)
        residues_dict[res_key].atoms.append(atom)
        residues_dict[res_key].atom_types.append(atp)
        residues_dict[res_key].pdb_lines.append(pdb)
        residues_dict[res_key].terminal = (term)
    return list(residues_dict.values()),list(non_protein_residue_info.values())

def GetActiveResType(listactiveres, ires_other, totalnumberofresidues) -> List[str]:
    """Function that gets the type of active resides from another list"""
    activerestype = []
    searchidx = 0
    for ires in range(len(listactiveres)):
        activerestype.append("ATOM")
        for jres in range(searchidx,totalnumberofresidues):
            if (listactiveres[ires][0] == ires_other[jres][1]) and (listactiveres[ires][1] == ires_other[jres][2]):
                activerestype[ires] = "HETATM"
                searchidx = jres
                break
    return activerestype

def GetConnectors(listactiveres, activerestype) -> List[Tuple[int,str]]:
    """Function that extracts the connecting residues needed for a certain molecular fragment"""
    connectors = []
    for idres in range(1,len(listactiveres)):
        if (listactiveres[idres - 1][1] == listactiveres[idres][1]) and (activerestype[idres - 1] == "ATOM") and (activerestype[idres] == "ATOM"):
            for jdres in range(listactiveres[idres - 1][0] + 1,listactiveres[idres][0]):
                connectors.append([jdres,listactiveres[idres][1]])
    return connectors

def AssembleSystem(residues, protonate: bool = False, tempfolder: str = "") -> Tuple[List[Tuple[float, float, float]], List[str], List[str], List[str], List[int], List[str], List[str], float]:
    """Assembling a single geometry from a collection of residues"""
    cut_geom = []
    cut_atmlist = []
    cut_atomtypelist = []
    cut_residuelist = []
    cut_residuenumbers = []
    cut_chainidentifiers = []
    cut_residuetypes = []
    total_charge = 0
    if not tempfolder.endswith("/"): tempfolder += "/"
    pdbname = tempfolder + "temp.pdb"
    xyznameH = pdbname.replace(".pdb","H.xyz")
    sdfname = tempfolder + "temp.sdf"
    xyzname = tempfolder + "temp.xyz"
    for residue in residues:
        #add base atoms
        if not protonate:
            for idatm in range(len(residue.atoms)):
                cut_geom.append(residue.geometry[idatm])
                cut_atomtypelist.append(residue.atom_types[idatm])
                cut_atmlist.append(residue.atoms[idatm])
                cut_residuelist.append(residue.name)
                cut_residuenumbers.append(residue.number)
                cut_chainidentifiers.append(residue.chain_id)
                cut_residuetypes.append(residue.residue_type)
            #add cap
            for idatm in range(len(residue.cap_atoms)):
                atom_type2add = "O" if residue.cap_atoms[idatm] == "O" else "C"
                cut_geom.append(residue.cap_geometry[idatm])
                cut_atomtypelist.append(atom_type2add)
                cut_atmlist.append(residue.cap_atoms[idatm])
                cut_residuelist.append(residue.name)
                cut_residuenumbers.append(residue.number)
                cut_chainidentifiers.append(residue.chain_id)
                cut_residuetypes.append(residue.residue_type)
        else: 
            ismetal,local_charge = ION.IsMetal(residue.name)
            #temp containers with geometry, needed because of the cap
            temp_geometry = []
            temp_atoms = []
            temp_atom_type = []
            temp_residue = []
            temp_residue_number = []
            temp_chain_id = []
            temp_residue_type = []
            #populate them with the base geometry
            for idatm in range(len(residue.atoms)):
                #build temp geometry to protonate
                temp_geometry.append(residue.geometry[idatm])
                temp_atoms.append(residue.atoms[idatm])
                temp_atom_type.append(residue.atom_types[idatm])
                temp_residue.append(residue.name)
                temp_residue_number.append(residue.number)
                temp_chain_id.append(residue.chain_id)
                temp_residue_type.append(residue.residue_type)
                #add to the global geometry
                cut_geom.append(residue.geometry[idatm])
                cut_atmlist.append(residue.atoms[idatm])
                cut_atomtypelist.append(residue.atom_types[idatm])
                cut_residuelist.append(residue.name)
                cut_residuenumbers.append(residue.number)
                cut_chainidentifiers.append(residue.chain_id)
                cut_residuetypes.append(residue.residue_type)
            #add cap
            for idatm in range(len(residue.cap_atoms)):
                #same logic
                temp_geometry.append(residue.cap_geometry[idatm])
                temp_atoms.append(residue.cap_atoms[idatm])
                atom_type2add = "O" if residue.cap_atoms[idatm] == "O" else "C"
                temp_atom_type.append(atom_type2add)
                temp_residue.append("GLY")
                temp_residue_number.append(residue.number - 1)
                temp_chain_id.append(residue.chain_id)
                temp_residue_type.append(residue.residue_type)
                cut_geom.append(residue.cap_geometry[idatm])
                cut_atmlist.append(residue.cap_atoms[idatm])
                cut_atomtypelist.append(atom_type2add)
                cut_residuelist.append(residue.name)
                cut_residuenumbers.append(residue.number)
                cut_chainidentifiers.append(residue.chain_id)
                cut_residuetypes.append(residue.residue_type)
            if AA.IsBaseAminoAcid(residue.name) or AA.IsExtendedAminoAcid(residue.name) or AA.IsPostTranslationalModification(residue.name):
                #now we write a temp PDB file for protonating
                GEOM.WriteGemmiPDB(temp_geometry,temp_atoms,temp_atom_type,temp_residue,temp_residue_number,temp_chain_id,temp_residue_type,pdbname)
                local_charge = VSEPR.ProteinProtonation(pdbname,xyznameH.replace(".xyz",""))
                total_charge += local_charge
                #get the new geometry
                atomsH,geometryH = GEOM.ReadXYZ(xyznameH)
                #add the new atoms/protons
                atoms_added = len(atomsH) - len(temp_atoms)
                for idatm in range(atoms_added):
                    cut_geom.append(geometryH[idatm + len(temp_atoms)])
                    cut_atmlist.append(atomsH[idatm + len(temp_atoms)])
                    cut_atomtypelist.append("H")
                    cut_residuelist.append(residue.name)
                    cut_residuenumbers.append(residue.number)
                    cut_chainidentifiers.append(residue.chain_id)
                    cut_residuetypes.append(residue.residue_type)
            elif AA.IsExternalCap(residue.name):
                #to be implemented
                a = 1
            elif ismetal:
                #do nothing here, no cap to add either
                total_charge += local_charge
            else: 
                #use the default VSEPR protonation
                if MON.MonomerExistsInLibrary(residue.name): 
                    #protonate based on base info and topology from PDB
                    path2monomerfile = MON.GetMonomerPath(residue.name)
                    connectivity = MON.ExtractCCDconnectivityGEMMI(path2monomerfile)
                    local_charges = AA.GetLocalCharges(residue.name,residue.atom_types)
                    GEOM.WriteSDF(residue.atom_types,residue.atoms,residue.geometry,connectivity,local_charges,sdfname)
                    print("ADD VSEPR-SDF PROTONATION HERE")
                    local_charge = VSEPR.VSEPRProtonation(sdfname,xyznameH.replace(".xyz",""))
                    total_charge += local_charge
                else: 
                    #do the generalised protonation
                    GEOM.Write2XYZ(temp_atoms,temp_geometry,sdfname.replace(".sdf",".xyz"))
                    local_charge = VSEPR.VSEPRProtonation(xyzname,xyznameH.replace(".xyz",""))
                #get the new geometry
                atomsH,geometryH = GEOM.ReadXYZ(xyznameH)
                #add the new atoms/protons
                atoms_added = len(atomsH) - len(temp_atoms)
                for idatm in range(atoms_added):
                    cut_geom.append(geometryH[idatm + len(temp_atoms)])
                    cut_atmlist.append(atomsH[idatm + len(temp_atoms)])
                    cut_atomtypelist.append("H")
                    cut_residuelist.append(residue.name)
                    cut_residuenumbers.append(residue.number)
                    cut_chainidentifiers.append(residue.chain_id)
                    cut_residuetypes.append(residue.residue_type)
                    total_charge += local_charge
    #remove temp files
    if os.path.isfile(pdbname): os.remove(pdbname)
    if os.path.isfile(sdfname): os.remove(sdfname)
    if os.path.isfile(xyzname): os.remove(xyzname)
    if os.path.isfile(xyzname): os.remove(xyznameH)
    if protonate:
        #if we protonated, there might be clashes, so we need to address this
        import copy
        temp_geometry = []
        temp_atoms = []
        temp_atom_type = []
        temp_residue = []
        temp_residue_number = []
        temp_chain_id = []
        temp_residue_type = []
        for idatm in range(len(cut_atmlist)):
            if (cut_atmlist[idatm] != "H"):
                if (cut_atmlist[idatm] == "C") or (cut_atmlist[idatm] == "O"):
                    clash = False
                    for idbtm in range(idatm):
                        Rvec = [cut_geom[idatm][0] - cut_geom[idbtm][0],cut_geom[idatm][1] - cut_geom[idbtm][1],cut_geom[idatm][2] - cut_geom[idbtm][2]]
                        R2 = Rvec[0]*Rvec[0] + Rvec[1]*Rvec[1] + Rvec[2]*Rvec[2]
                        if R2 < 1.0e-1:
                            #repeated atom
                            clash = True
                            break
                    if clash: continue
                temp_geometry.append(cut_geom[idatm])
                temp_atoms.append(cut_atmlist[idatm])
                temp_atom_type.append(cut_atomtypelist[idatm])
                temp_residue.append(cut_residuelist[idatm])
                temp_residue_number.append(cut_residuenumbers[idatm])
                temp_chain_id.append(cut_chainidentifiers[idatm])
                temp_residue_type.append(cut_residuetypes[idatm])
            else: 
                #check for clashes
                clash = False
                for idbtm in range(len(cut_atmlist)):
                    if idatm == idbtm: continue
                    Rvec = [cut_geom[idatm][0] - cut_geom[idbtm][0],cut_geom[idatm][1] - cut_geom[idbtm][1],cut_geom[idatm][2] - cut_geom[idbtm][2]]
                    R2 = Rvec[0]*Rvec[0] + Rvec[1]*Rvec[1] + Rvec[2]*Rvec[2]
                    if R2 < 0.5:
                        #clash
                        clash = True
                        break
                if not clash:
                    temp_geometry.append(cut_geom[idatm])
                    temp_atoms.append(cut_atmlist[idatm])
                    temp_atom_type.append(cut_atomtypelist[idatm])
                    temp_residue.append(cut_residuelist[idatm])
                    temp_residue_number.append(cut_residuenumbers[idatm])
                    temp_chain_id.append(cut_chainidentifiers[idatm])
                    temp_residue_type.append(cut_residuetypes[idatm])
        #now copy back
        cut_geom = copy.deepcopy(temp_geometry)
        cut_atmlist = copy.deepcopy(temp_atoms)
        cut_atomtypelist = copy.deepcopy(temp_atom_type)
        cut_residuelist = copy.deepcopy(temp_residue)
        cut_residuenumbers = copy.deepcopy(temp_residue_number)
        cut_chainidentifiers = copy.deepcopy(temp_chain_id)
        cut_residuetypes = copy.deepcopy(temp_residue_type)
    return (cut_geom,cut_atmlist,cut_atomtypelist,cut_residuelist,cut_residuenumbers,cut_chainidentifiers,cut_residuetypes,total_charge)

def rotate_point_around_axis(point, axis_p1, axis_p2, angle_deg = 180.0):
    """Rotate `point` around the axis defined by axis_p1 -> axis_p2.
       Coordinates should be array-like: [x, y, z]."""
    p = np.asarray(point,dtype = float)
    a = np.asarray(axis_p1,dtype = float)
    b = np.asarray(axis_p2,dtype = float)
    axis = b - a
    norm = np.linalg.norm(axis)
    if norm == 0:
        raise ValueError("Rotation axis has zero length")
    u = axis/norm
    v = p - a
    theta = np.deg2rad(angle_deg)
    #Rodrigues rotation formula
    v_rot = (v*np.cos(theta)
             + np.cross(u,v)*np.sin(theta)
             + u*np.dot(u,v)*(1.0 - np.cos(theta)))
    return (a + v_rot).tolist()

def flip_hydroxyl_groups(residue) -> ResidueData:
    """Function that flips hydroxyl groups to find additional interaction patterns."""
    #find the oxygen and carbon atoms
    idx_oxygen = -1
    idx_carbon = -1
    if residue.name == "SER":
        for iatm in range(len(residue.atom_types)):
            atom = residue.atom_types[iatm]
            if atom == "OG":
                idx_oxygen = iatm
            elif atom == "CB":
                idx_carbon = iatm
            if (idx_oxygen != -1) and (idx_carbon != -1): break
    elif residue.name == "THR":
        for iatm in range(len(residue.atom_types)):
            atom = residue.atom_types[iatm]
            if atom == "OG1":
                idx_oxygen = iatm
            elif atom == "CB":
                idx_carbon = iatm
            if (idx_oxygen != -1) and (idx_carbon != -1): break
    elif residue.name == "TYR":
        for iatm in range(len(residue.atom_types)):
            atom = residue.atom_types[iatm]
            if atom == "OH":
                idx_oxygen = iatm
            elif atom == "CZ":
                idx_carbon = iatm
            if (idx_oxygen != -1) and (idx_carbon != -1): break
    if idx_oxygen >= 0 and idx_carbon >= 0:
        o_coord = residue.geometry[idx_oxygen]
        c_coord = residue.geometry[idx_carbon]
        #now find the closest proton
        idx_hydrogen = -1
        for iatm in range(len(residue.processed_atoms)):
            atom = residue.processed_atoms[iatm]
            if atom == "H":
                h_coord = residue.processed_geometry[iatm]
                dist_vec = [o_coord[0] - h_coord[0],o_coord[1] - h_coord[1],o_coord[2] - h_coord[2]]
                R2 = dist_vec[0]*dist_vec[0] + dist_vec[1]*dist_vec[1] + dist_vec[2]*dist_vec[2]
                if R2 < 1.21:          #1.1*1.1
                    idx_hydrogen = iatm
                    break
        #now rotate
        if idx_hydrogen != -1:
            h_coord = residue.processed_geometry[idx_hydrogen]
            new_H = rotate_point_around_axis(point = h_coord,axis_p1 = c_coord,axis_p2 = o_coord,angle_deg = 180.0)
            residue.processed_geometry[idx_hydrogen] = new_H
            residue.name += "_flip"                      #mark as flipped
            residue.residue_id += "f"                    #mark as flipped
            if str(residue.protonated_file).endswith(".xyz"):
                residue.protonated_file = str(residue.protonated_file).replace(".xyz","_flip.xyz")
                GEOM.Write2XYZ(residue.processed_atoms,residue.processed_geometry,str(residue.protonated_file))
            elif str(residue.protonated_file).endswith(".sdf"):
                connectivity,atomtype,_,_,local_charges = GEOM.ReadSDFWithTopology(str(residue.protonated_file))
                residue.protonated_file = str(residue.protonated_file).replace(".sdf","_flip.sdf")
                GEOM.WriteSDF(atomtype,residue.processed_atoms,residue.processed_geometry,connectivity,local_charges,str(residue.protonated_file),skip_proton = False)
            return residue
        else: 
            return None
    else: 
        return None