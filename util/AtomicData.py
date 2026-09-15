# ------------------------------------------------------------
# Util with atom-related information
# ------------------------------------------------------------

from typing import Dict,Optional,Union
import numpy as np

#covalent radii + neighbor helpers (used to disambiguate "standard" states when --no-cap is active and backbone protons can compensate sidechain charge)
#reference radii (Å): B. Cordero et al., 2008, doi.org/10.1039/B801115J 
_COVALENT_RADII_ANGSTROM: Dict[int,float] = {
    1: 0.31,   #H
    2: 0.28,   #He
    3: 1.28,   #Li
    4: 0.96,   #Be
    5: 0.84,   #B
    6: 0.76,   #C
    7: 0.71,   #N
    8: 0.66,   #O
    9: 0.57,   #F
    10: 0.58,  #Ne
    11: 1.66,  #Na
    12: 1.41,  #Mg
    13: 1.21,  #Al
    14: 1.11,  #Si
    15: 1.07,  #P
    16: 1.05,  #S
    17: 1.02,  #Cl
    18: 1.06,  #Ar
    19: 2.03,  #K
    20: 1.76,  #Ca
    21: 1.70,  #Sc
    22: 1.60,  #Ti
    23: 1.53,  #V
    24: 1.39,  #Cr
    25: 1.61,  #Mn
    26: 1.52,  #Fe
    27: 1.50,  #Co
    28: 1.24,  #Ni
    29: 1.32,  #Cu
    30: 1.22,  #Zn
    31: 1.22,  #Ga
    32: 1.20,  #Ge
    33: 1.19,  #As
    34: 1.20,  #Se
    35: 1.20,  #Br
    36: 1.16,  #Kr
    37: 2.20,  #Rb
    38: 1.95,  #Sr
    39: 1.90,  #Y
    40: 1.75,  #Zr
    41: 1.64,  #Nb
    42: 1.54,  #Mo
    43: 1.47,  #Tc
    44: 1.46,  #Ru
    45: 1.42,  #Rh
    46: 1.39,  #Pd
    47: 1.45,  #Ag
    48: 1.44,  #Cd
    49: 1.42,  #In
    50: 1.39,  #Sn
    51: 1.39,  #Sb
    52: 1.38,  #Te
    53: 1.39,  #I
    54: 1.40,  #Xe
    55: 2.44,  #Cs
    56: 2.15,  #Ba
    57: 2.07,  #La
    58: 2.04,  #Ce
    59: 2.03,  #Pr
    60: 2.01,  #Nd
    61: 1.99,  #Pm
    62: 1.98,  #Sm
    63: 1.98,  #Eu
    64: 1.96,  #Gd
    65: 1.94,  #Tb
    66: 1.92,  #Dy
    67: 1.92,  #Ho
    68: 1.89,  #Er
    69: 1.90,  #Tm
    70: 1.87,  #Yb
    71: 1.87,  #Lu
    72: 1.75,  #Hf
    73: 1.70,  #Ta
    74: 1.62,  #W
    75: 1.51,  #Re
    76: 1.44,  #Os
    77: 1.41,  #Ir
    78: 1.36,  #Pt
    79: 1.36,  #Au
    80: 1.32,  #Hg
    81: 1.45,  #Tl
    82: 1.46,  #Pb
    83: 1.48,  #Bi
    84: 1.40,  #Po
    85: 1.50,  #At
    86: 1.50,  #Rn
    87: 2.60,  #Fr
    88: 2.21,  #Ra
    89: 2.15,  #Ac
    90: 2.06,  #Th
    91: 2.00,  #Pa
    92: 1.96,  #U
    93: 1.90,  #Np
    94: 1.87,  #Pu
    95: 1.80,  #Am
    96: 1.69,  #Cm
}

_ELEMENT_TO_Z: Dict[str,int] = {
    "H": 1,
    "He": 2,
    "Li": 3,
    "Be": 4,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "Ne": 10,
    "Na": 11,
    "Mg": 12,
    "Al": 13,
    "Si": 14,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "Ar": 18,
    "K": 19,
    "Ca": 20,
    "Sc": 21,
    "Ti": 22,
    "V": 23,
    "Cr": 24,
    "Mn": 25,
    "Fe": 26,
    "Co": 27,
    "Ni": 28,
    "Cu": 29,
    "Zn": 30,
    "Ga": 31,
    "Ge": 32,
    "As": 33,
    "Se": 34,
    "Br": 35,
    "Kr": 36,
    "Rb": 37,
    "Sr": 38,
    "Y": 39,
    "Zr": 40,
    "Nb": 41,
    "Mo": 42,
    "Tc": 43,
    "Ru": 44,
    "Rh": 45,
    "Pd": 46,
    "Ag": 47,
    "Cd": 48,
    "In": 49,
    "Sn": 50,
    "Sb": 51,
    "Te": 52,
    "I": 53,
    "Xe": 54,
    "Cs": 55,
    "Ba": 56,
    "La": 57,
    "Ce": 58,
    "Pr": 59,
    "Nd": 60,
    "Pm": 61,
    "Sm": 62,
    "Eu": 63,
    "Gd": 64,
    "Tb": 65,
    "Dy": 66,
    "Ho": 67,
    "Er": 68,
    "Tm": 69,
    "Yb": 70,
    "Lu": 71,
    "Hf": 72,
    "Ta": 73,
    "W": 74,
    "Re": 75,
    "Os": 76,
    "Ir": 77,
    "Pt": 78,
    "Au": 79,
    "Hg": 80,
    "Tl": 81,
    "Pb": 82,
    "Bi": 83,
    "Po": 84,
    "At": 85,
    "Rn": 86,
    "Fr": 87,
    "Ra": 88,
    "Ac": 89,
    "Th": 90,
    "Pa": 91,
    "U": 92,
    "Np": 93,
    "Pu": 94,
    "Am": 95,
    "Cm": 96
}

def _covalent_radius_angstrom(elem: Union[str,int]) -> Optional[float]:
    """return covalent radius (Å) for a small set of elements used in neighbour checks."""
    if isinstance(elem,(int,np.integer)):
        return _COVALENT_RADII_ANGSTROM.get(int(elem))
    e = str(elem).strip().capitalize()
    z = _ELEMENT_TO_Z.get(e)
    if z is None:
        return None
    return _COVALENT_RADII_ANGSTROM.get(z)

def _is_hydrogen_symbol(symbol: str) -> bool:
    """Match the current sanitizer behaviour: explicit H only."""
    return str(symbol).strip().upper() == "H"

__all__ = [
    "_covalent_radius_angstrom",
    "_COVALENT_RADII_ANGSTROM",
    "_ELEMENT_TO_Z",
    "_is_hydrogen_symbol"
]
