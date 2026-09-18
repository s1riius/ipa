# ------------------------------------------------------------
# Functions related to ions and metals data
# ------------------------------------------------------------

from __future__ import annotations
import re
from typing import Tuple

#default charges commonly used in protein PDB contexts (extend as needed)
ION_CHARGES = {
    #halides
    "F": -1, "CL": -1, "BR": -1, "I": -1,
    #alkali metals
    "LI": 1, "NA": 1, "K": 1, "RB": 1, "CS": 1,
    #alkaline earth metals
    "MG": 2, "CA": 2, "SR": 2, "BA": 2,
    #common transition metals (typical defaults)
    "ZN": 2, "MN": 2, "FE": 2, "CO": 2, "NI": 2, "CU": 2, "V": 2,
    #other occasionally seen ions
    "AL": 3, "CD": 2, "HG": 2, "AG": 1
}

#explicit residue-name overrides (deterministic, do not “guess” too much)
EXPLICIT_OVERRIDES = {
    "FE3": 3, "FE2": 2,
    "CU1": 1, "CU2": 2,
    "MN2": 2, "ZN2": 2,
    "CA2": 2, "MG2": 2,
    "CL1": -1, "BR1": -1, "I1": -1
}

_RE_ELEM_SIGN = re.compile(r"^([A-Z]{1,2})([+-])$")                  # CL-, NA+
_RE_ELEM_SIGN_NUM = re.compile(r"^([A-Z]{1,2})([+-])(\d+)$")         # FE+3
_RE_ELEM_NUM_SIGN = re.compile(r"^([A-Z]{1,2})(\d+)([+-])$")         # FE3+
_RE_ELEM_NUM = re.compile(r"^([A-Z]{1,2})(\d+)$")                    # FE3, CA2, MN1

def _norm(name: str) -> str:
    return (name.strip().upper() or "")

# ------------------------------------------------------------
# Ion stuff
# ------------------------------------------------------------
def IsIon(residue: str) -> Tuple[bool,int]:
    """Determine whether a residue name corresponds to a simple ion and return (is_ion,charge).
       Supported:
         - Plain symbols: 'MN', 'ZN', 'CL', 'NA', ...
         - Annotated: 'FE3', 'FE3+', 'FE+3', 'CL-', 'CA2+', 'CU1', ..."""
    res = _norm(residue)
    if not res: return False,0
    #1) explicit overrides first
    if res in EXPLICIT_OVERRIDES:
        return True,EXPLICIT_OVERRIDES[res]
    #2) parse explicit sign + number
    m = _RE_ELEM_SIGN_NUM.match(res)  #FE+3
    if m:
        elem,sign,n = m.group(1),m.group(2),int(m.group(3))
        return True,(n if sign == "+" else -n)
    m = _RE_ELEM_NUM_SIGN.match(res)  #FE3+
    if m:
        elem,n,sign = m.group(1),int(m.group(2)),m.group(3)
        return True,(n if sign == "+" else -n)
    #3) parse simple +/- (assume magnitude 1)
    m = _RE_ELEM_SIGN.match(res)  # CL-
    if m:
        elem,sign = m.group(1),m.group(2)
        return True,(1 if sign == "+" else -1)
    #4) parse elem + number, but be conservative:
    #   accept only if the number matches the known default magnitude
    #   (prevents MN1 => +1; instead falls back to MN => +2)
    m = _RE_ELEM_NUM.match(res)  # FE3, CA2, MN1, ...
    if m:
        elem,n = m.group(1),int(m.group(2))
        if elem in ION_CHARGES and abs(ION_CHARGES[elem]) == n:
            return True,ION_CHARGES[elem]
        #otherwise: do NOT guess; fall back to plain elem lookup
    #5) plain lookup
    if res in ION_CHARGES:
        return True,ION_CHARGES[res]
    return False,0

# ------------------------------------------------------------
# Metals
# ------------------------------------------------------------
def IsMetal(residue: str) -> Tuple[bool,int]:
    """Detect whether `residue` is a metal; returns (is_metal: bool, charge: int)."""
    ismetal = False
    charge = 0
    res = residue.upper().strip()
    if (res == "LI") or (res == "NA") or (res == "K") or (res == "RB") or (res == "CS"):
        ismetal = True
        charge = 1
    elif (res == "MG") or (res == "CA") or (res == "SR") or (res == "BA") or (res == "ZN") or (res == "MN"):
        ismetal = True
        charge = 2
    elif (res == "FE") or (res == "CO") or (res == "NI") or (res == "CU") or (res == "V"):
        ismetal = True
        charge = 2
    elif (res == "AL"):
        ismetal = True
        charge = 3
    elif (res == "AG"):
        ismetal = True
        charge = 1
    return ismetal,charge

__all__ = [
    "IsIon",
    "IsMetal"
]
