"""
Utility modules for In-Pocket Analysis framework.
This package contains various utility functions and classes for quantum mechanical
calculations, structure manipulation, and protein-ligand interaction analysis.
"""

# Re-export all utility modules to make them available when importing from in_pocket.util
from util.BaseAminoAcids import *
from util.Formaldehyde import *
from util.Geometry import *
from util.Ions import *
from util.Protoro import *
from util.QMwrapper import *
from util.RMSD import *
from util.Tautoro import *
from util.VSEPR import *

# Define package-wide constants
VERSION = "1.0.0"

# Define what gets imported with "from in_pocket.util import *"
__all__ = [
    # BaseAminoAcids
    "IsBaseAminoAcid",
    # Formaldehyde functions
    "Eformaldehyde",
    # Geometry functions
    "ReadGeometry",
    "Write2XYZ",
    "Write2PDB",
    "PDB2XYZ",
    "WriteGemmiPDB",
    "ReadPDB_base",
    "ReadCIFStructure",
    "WriteSDF",
    "MOL2toXYZ",
    "SelectOneAltLocOnly",
    "ReCenterAndWrite",
    # Ion utilities
    "IsIon",
    "IsMetal",
    # Protonation utilities
    "print_header_protoro",
    "RDKitProtonation",
    "ObabelProtonation",
    "ProteinProtonation",
    # QM wrapper functions
    "RunIP",
    "ParseOutput",
    # RMSD calculation
    "RMSDcalculation",
    # Tautomer utilities
    "process_tautomers",
    # VSEPR utilities
    "VSEPRProtonation",
    "CCDProtonation",
    "ExtractLigandConnectivity",
]
