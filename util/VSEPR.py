# ------------------------------------------------------------
# VSEPR-related utils
# ------------------------------------------------------------

import re
import subprocess
from pathlib import Path
from collections import Counter
from typing import Dict,Optional,Tuple

#import util.Formatting as FORMAT
import util.Geometry as GEOM

def _resolve_executable(exe_name: str) -> str:
    cwd_exe = Path(".") / exe_name
    if cwd_exe.is_file():
        return str(cwd_exe.resolve())
    repo_exe = Path(__file__).resolve().parents[1] / exe_name
    if repo_exe.is_file():
        return str(repo_exe)
    #preserve the old intended failure target, but make it explicit/absolute.
    return str(cwd_exe.resolve())

def CCDProtonation(inputfile, pdbfile, ccdfile, outputfile):
    #function to protonate using a CCD file, to get back a xyz and sdf file, along with molecular charge
    exefile = _resolve_executable("CCD_Protonate.exe")
    proc = subprocess.Popen([exefile,inputfile,pdbfile,ccdfile,outputfile],stdout = subprocess.PIPE,stderr = subprocess.PIPE)
    out,_ = proc.communicate()
    outputfile = out.decode("utf-8").split("\n")
    totalcharge = int(outputfile[0])
    return totalcharge

def _parse_last_int_from_stdout(stdout_text: str) -> int:
    """Return the last standalone integer found in stdout."""
    int_tokens = []
    for line in stdout_text.splitlines():
        s = line.strip()
        if re.fullmatch(r"[+-]?\d+", s):
            int_tokens.append(int(s))
    if not int_tokens:
        raise ValueError(f"Could not parse integer charge from VSEPR stdout: {stdout_text!r}")
    return int(int_tokens[-1])

def _xyz_atom_summary(xyz_path: str) -> Dict[str, object]:
    """Return a compact atom summary for an XYZ file."""
    atoms,_ = GEOM.ReadXYZ(xyz_path)
    elems = [str(a).strip().upper() for a in atoms]
    heavy = [el for el in elems if el != "H"]
    return {"n_atoms": len(elems),
            "n_h": sum(1 for el in elems if el == "H"),
            "n_heavy": len(heavy),
            "heavy_signature": dict(Counter(heavy)),}

def VSEPRProtonation(inputfile, outputfile):
    """connection to VSEPR protonation code"""
    exefile = _resolve_executable("VSEPR_Protonate.exe")
    proc = subprocess.Popen([exefile,inputfile,outputfile],stdout = subprocess.PIPE,stderr = subprocess.PIPE)
    out,_ = proc.communicate()
    #return totalcharge
    stdout_text = out.decode("utf-8",errors = "replace")
    totalcharge = _parse_last_int_from_stdout(stdout_text)
    return int(totalcharge)

def VSEPRChargeProbeXYZ(inputfile: str, outputfile: str) -> Tuple[Optional[int], Dict[str, object]]:
    """Use VSEPR_Protonate.exe as a charge probe only. The executable may add/remove hydrogens. We therefore:
         - read the raw charge printed for the VSEPR output geometry
         - compare input/output XYZ atom compositions
         - reject the probe if the heavy-atom framework changes
         - correct the raw charge back to the original input ligand charge using:
               q_input = q_vsepr - (H_out - H_in)"""
    exefile = _resolve_executable("VSEPR_Protonate.exe")
    proc = subprocess.Popen([exefile,inputfile,outputfile],stdout = subprocess.PIPE,stderr = subprocess.PIPE)
    out,err = proc.communicate()
    stdout_text = out.decode("utf-8",errors = "replace")
    stderr_text = err.decode("utf-8",errors = "replace")
    meta: Dict[str,object] = {"input_xyz": inputfile,
                              "output_base": outputfile,
                              "stdout": stdout_text,
                              "stderr": stderr_text,
                              "returncode": proc.returncode,}
    if proc.returncode not in (0,None):
        return None,meta
    raw_charge = _parse_last_int_from_stdout(stdout_text)
    output_xyz = outputfile if str(outputfile).lower().endswith(".xyz") else (str(outputfile) + ".xyz")
    input_summary = _xyz_atom_summary(inputfile)
    output_summary = _xyz_atom_summary(output_xyz)
    meta.update({"output_xyz": output_xyz,
                 "charge_vsepr": int(raw_charge),
                 "input_atoms": int(input_summary["n_atoms"]),
                 "output_atoms": int(output_summary["n_atoms"]),
                 "input_h": int(input_summary["n_h"]),
                 "output_h": int(output_summary["n_h"]),
                 "input_heavy_atoms": int(input_summary["n_heavy"]),
                 "output_heavy_atoms": int(output_summary["n_heavy"]),
                 "input_heavy_signature": input_summary["heavy_signature"],
                 "output_heavy_signature": output_summary["heavy_signature"],})
    heavy_ok = (input_summary["n_heavy"] == output_summary["n_heavy"] and input_summary["heavy_signature"] == output_summary["heavy_signature"])
    meta["heavy_atoms_unchanged"] = bool(heavy_ok)
    if not heavy_ok:
        return None,meta
    delta_h = int(output_summary["n_h"]) - int(input_summary["n_h"])
    adjusted_charge = int(raw_charge) - int(delta_h)
    meta["delta_h"] = int(delta_h)
    meta["charge_adjusted"] = int(adjusted_charge)
    return int(adjusted_charge),meta

def ProteinProtonation(inputfile, outputfile):
    #connection to VSEPR protonation code
    exefile = _resolve_executable("Protein_Protonate.exe")
    proc = subprocess.Popen([exefile,inputfile,outputfile],stdout = subprocess.PIPE,stderr = subprocess.PIPE)
    out,_ = proc.communicate()
    outputfile = out.decode("utf-8").split("\n")
    totalcharge = int(outputfile[0])
    return totalcharge

def VSEPRSDFMOL2Protonation(inputfile, outputfile):
    #connection to VSEPR protonation code
    exefile = _resolve_executable("SDFMOL2_Protonate.exe")
    proc = subprocess.Popen([exefile,inputfile,outputfile],stdout = subprocess.PIPE,stderr = subprocess.PIPE)
    out,_ = proc.communicate()
    outputfile = out.decode("utf-8").split("\n")
    totalcharge = int(outputfile[0])
    return totalcharge