# ------------------------------------------------------------
# Capping util
# ------------------------------------------------------------

import os
import copy
import subprocess
import uuid
from pathlib import Path
import numpy as np
import util.Geometry as GEOM
import util.UnitConversion as UC

rad2deg = UC.rad2degree

def Eformaldehyde(solvent):
    #function that returns the energy and location of the geometry file for the free formaldehyde
    energy = 0.0
    path2file = ""
    if solvent == "gas":
        energy = -7.17564755863244930
        path2file = "util/formaldehyde/CH2O_gas"
    elif solvent == "water":
        energy = -7.17690971163739633
        path2file = "util/formaldehyde/CH2O_water"
    elif solvent == "octanol wet":
        energy = -7.17755624880021337
        path2file = "util/formaldehyde/CH2O_octoh-wet"
    elif solvent == "octanol":
        energy = -7.17792728879785802
        path2file = "util/formaldehyde/CH2O_octoh"
    elif solvent == "methanol":
        energy = -7.17842323339754618
        path2file = "util/formaldehyde/CH2O_meoh"
    elif solvent == "ethanol":
        energy = -7.17869038781591762
        path2file = "util/formaldehyde/CH2O_etoh"
    elif solvent == "dmso":
        energy = -7.17990281296664357
        path2file = "util/formaldehyde/CH2O_dmso"
    elif solvent == "dmf":
        energy = -7.18024518457347583
        path2file = "util/formaldehyde/CH2O_dmf"
    elif solvent == "acetone":
        energy = -7.17981494783519469
        path2file = "util/formaldehyde/CH2O_acetone"
    return energy,path2file

def _resolve_align_exe() -> str:
    """Resolve align.exe relative to the IPA repo, not the current working directory."""
    util_dir = Path(__file__).resolve().parent
    repo_dir = util_dir.parent                                 #folder containing in_pocket.py
    candidates = [util_dir / "align.exe",                      #util/align.exe
                  repo_dir / "align.exe",                      #align.exe next to in_pocket.py
                  repo_dir / "exe" / "auto_compiled" / "align.exe",
                  repo_dir / "core" / "cpp_util" / "align.exe",
                  repo_dir / "programs" / "align.exe"]
    for p in candidates:
        if p.is_file():
            return str(p)
    raise FileNotFoundError("align.exe not found. Tried:\n"
                            + "\n".join(f"  - {c}" for c in candidates)
                            + "\nFix: compile it or place it next to in_pocket.py or inside util/.")

def AddAtom2Cart(dist, angle, dihed, Acoord, Bcoord, Ccoord):
    #function to calculate cartesian coordinates of a new atom, A, given internal coordinates
    #dist - distance to a certain atom B
    #angle - the angle made between A-B and another atom C
    #dihed - the torsional angle made between A-B-C and another atom D
    sinT = np.sin(angle)
    cosT = np.cos(angle)
    sinP = np.sin(dihed)
    cosP = np.cos(dihed)
    xcoord = dist*cosT
    ycoord = dist*cosP*sinT
    zcoord = dist*sinP*sinT
    BCdist = np.array(Bcoord) - np.array(Ccoord)
    ACdist = np.array(Acoord) - np.array(Bcoord)
    norm = 1.0/np.linalg.norm(ACdist)
    ACdist *= norm
    BCAC = np.cross(BCdist,ACdist)
    norm = 1.0/np.linalg.norm(BCAC)
    BCAC *= norm
    BCACAC = np.cross(BCAC,ACdist)
    new_xcoord = Acoord[0] - ACdist[0]*xcoord + BCACAC[0]*ycoord + BCAC[0]*zcoord
    new_ycoord = Acoord[1] - ACdist[1]*xcoord + BCACAC[1]*ycoord + BCAC[1]*zcoord
    new_zcoord = Acoord[2] - ACdist[2]*xcoord + BCACAC[2]*ycoord + BCAC[2]*zcoord
    return [new_xcoord,new_ycoord,new_zcoord]

def AlignMolecules(structi1, structi2, structo):
    #wrapper around the program to align structures
    align_exe = _resolve_align_exe()
    proc = subprocess.Popen([align_exe,structi1,structi2,structo],stdout = subprocess.PIPE,stderr = subprocess.PIPE)
    out,err = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError("align.exe failed\n"
                           f"Command: {align_exe} {structi1} {structi2} {structo}\n"
                           f"stdout:\n{out.decode(errors = 'replace')}\n"
                           f"stderr:\n{err.decode(errors = 'replace')}\n")

def ReorientCap(OptCap, ProtCap):
    #function that does all the necessary math for reorienting the cap to match the protein environment
    #coordinates of extra-proton for the nasty math
    rOH = 1.1
    aHOC = 120.0*rad2deg
    dHOCH = 90.0*rad2deg
    opt_atoms,opt_geometry = GEOM.ReadXYZ(OptCap)
    prot_atoms,prot_geometry = GEOM.ReadXYZ(ProtCap)
    shift = copy.deepcopy(prot_geometry[0])
    Hcoord_opt = AddAtom2Cart(rOH,aHOC,dHOCH,opt_geometry[1],opt_geometry[0],opt_geometry[2])
    Hcoord_prot = AddAtom2Cart(rOH,aHOC,dHOCH,prot_geometry[1],prot_geometry[0],prot_geometry[2])
    newopt_geometry = copy.deepcopy(opt_geometry)
    newopt_geometry.append(Hcoord_opt)
    newprot_geometry = copy.deepcopy(prot_geometry)
    newprot_geometry.append(Hcoord_prot)
    newopt_atoms = copy.deepcopy(opt_atoms)
    newopt_atoms.append("H")
    newprot_atoms = copy.deepcopy(prot_atoms)
    newprot_atoms.append("H")
    temp_opt_filename = OptCap + f"_tmp_{uuid.uuid4().hex}"
    temp_prot_filename = ProtCap + f"_tmp_{uuid.uuid4().hex}"
    GEOM.Write2XYZ(newopt_atoms,newopt_geometry,temp_opt_filename)
    GEOM.Write2XYZ(newprot_atoms,newprot_geometry,temp_prot_filename)
    AlignMolecules(temp_prot_filename + ".xyz",temp_opt_filename + ".xyz",temp_prot_filename + "2")
    prot_atoms,prot_geometry = GEOM.ReadXYZ(temp_prot_filename + "2")
    new_prot_atoms = []
    new_prot_geometry = []
    shift[0] -= prot_geometry[0][0]
    shift[1] -= prot_geometry[0][1]
    shift[2] -= prot_geometry[0][2]
    for iatm in range(4):
        aux = [0.0,0.0,0.0]
        for idcoord in range(3):
            aux[idcoord] = prot_geometry[iatm][idcoord] + shift[idcoord]
        new_prot_atoms.append(prot_atoms[iatm])
        new_prot_geometry.append(aux)
    #at this point the geometries are aligned, though there might be a discrepance due to fine differences in the centers of mass
    #fix by shifting according to the carbon coordinates
    os.remove(temp_prot_filename + "2.xyz")
    os.remove(temp_opt_filename + ".xyz")
    os.remove(temp_prot_filename + ".xyz")
    newxyzfile = ProtCap                                 #temp_prot_filename + "2"
    GEOM.Write2XYZ(new_prot_atoms,new_prot_geometry,newxyzfile)
