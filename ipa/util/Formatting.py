# ------------------------------------------------------------
# functions related to general formatting of strings and so on
# ------------------------------------------------------------

import re
import numpy as np

def CleanTheLine(line):
    #function that removes all double spaces from string
    cleanline = line
    for idx in range(250):
        if "  " in cleanline: 
            line = cleanline.replace("  "," ")
            cleanline = line
        else: break
    return cleanline

def ReformatingString(string, expsz):
    #function that adds spaces to the beginning of a string to ensure it has a certain size
    new_string = string
    if string.startswith("-"): expsz += 1
    if len(string) != expsz:
        for iter in range(expsz):
            temp_string = " " + new_string
            new_string = temp_string
            if len(new_string) == expsz: break
    return new_string

# ------------------------------------------------------------------------------
# resolve util/monomers_library paths independent of current working dir.
# ------------------------------------------------------------------------------
def Format_Charge_State(charge: int) -> str:
    """Format an integer charge as a stable residue-state string ('0', '+1', '-2', ...)."""
    try:
        q = int(charge)
    except Exception:
        q = 0
    return "0" if q == 0 else f"{q:+d}"

def Parse_Charge_State(state) -> int:
    """Parse residue_state strings from ('+','--') or ('+1','-2','0','1') formats."""
    if state is None:
        return 0
    #numpy integer types (and plain ints)
    if isinstance(state,(int,np.integer)):
        return int(state)
    s = str(state).strip()
    if not s:
        return 0
    #legacy tokens like '+', '++', '---'
    if re.fullmatch(r"\++",s):
        return len(s)
    if re.fullmatch(r"-+",s):
        return -len(s)
    try:
        return int(s)
    except Exception:
        #last-resort fallbacks
        u = s.upper()
        if u in ("POS","PLUS","P"):
            return 1
        if u in ("NEG","MINUS","M"):
            return -1
        return 0

__all__ = [
    "CleanTheLine",
    "ReformatingString",
    "Format_Charge_State",
    "Parse_Charge_State"
]