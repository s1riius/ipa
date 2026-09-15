"""
In-Pocket Analysis: A protein-ligand interaction analysis tool.
v260402 - Added mmCIF altloc selection using Gemmi altloc/occupancy metadata
v260401 - Added pre-QM structure validation, SDF H-distance checks, and post-QM H-connectivity warnings.
v260331 - VSEPR for ligand charge calculation in the native mode.
v260330 - Added RDKit-indepentent util/sanitizer module.
v260327 - VSEPR for protein-like residues only, no tautomerisation support yet.
v260317 - Added --protonation original/native to preserve the input chemical state.
        - Fixed missing explicit H-bond connectivity in native SDF output.
        - Added --force-user-charge with unified active-moiety charge override handling.
v260225 - Ligand IP energies: aborts with a clear error if any energy_ip is missing/zero.
v260224 - New covalent-exclusion logic in find_neighbour_residues
v260223 - Modified find_neighbour_residues to skip +/- X residues (default 1) with any atom within covalent_threshold of the ligand.
        - Updated best_ligand.py.
v260220 - Modified --no-cap "standard residue" selection with fallback for total energy calculation.
v260218 - Switching from protonation-state tags to explicit integer charges
        - Restored debug messaging
v260217 - Auto-path for util/monomers_library
        - Removed Selenium Hack
        - single-atom fast-path
v260129 - Strip external/artificial caps present in input PDB.
v260127 - Improved automatic removal of input hydrogens.
        - V3000, non-kekulizing SDF handling.
        - Introduced fallback diagnostics.
v251227 - RDKit preprocessing: automatic removal of input hydrogens (optional --keep-input-h).
        - Selenium handling for RDKit (MSE/SEC → MET/CYS internally, reversible).
        - Improved error handling and logging with --debug.
v251219 - Clash analysis with EDDA (--edda).
        - PNG output (--png and --png-hq).
v251212 - Added residue-state export (CSV) for IPA reuse, with support for flips/protonation.
        - Stable ligand_states cache; --best-unique selects lowest-energy residue variants.
        - Per-ligand IP-optimised SDF output.
v251210 - IP geometry translation to match the COM of heavy atoms of starting geometry.
        - Cleaned QM parameter pipeline, removed redundant config.* fields.
v251209 - Added --best flags post-processing for generating summary files.
This module performs quantum mechanical calculations on protein-ligand systems
to analyze binding interactions, deformation energies, and structural changes.
"""

import argparse
import logging
import copy
import os
import re
import sys
import time
import math
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from enum import Enum, auto
from pathlib import Path
from typing import Dict,List,Optional,Tuple,Union
import numpy as np
import shutil

#import utility modules
import util.Residues as RES
import util.AtomicData as AD
import util.BaseAminoAcids as AA
import util.BaseNucleicAcids as NA
import util.Formaldehyde as FORM
import util.Geometry as GEOM
import util.Connectivity as CONNECT
import util.Ions as ION
import util.Formatting as FORMAT
import util.QMwrapper as QM
import util.RMSD as RMSD
import util.SolventList as SOLV
import util.VSEPR as VSEPR
import util.pKA as PKA
import util.Scoring as SCORE
import util.Sanitiser as SANITY
import util.UnitConversion as UC
import util.Protoro as Protoro
import util.Tautoro as Tautoro
import util.SmartsClean as SC
from util.cub2mrc import convert_cub_tree_to_mrc
from util.best_ligand import generate_best_outputs,cleanup_outputs

# ------------------------------------------------------------------------------
# resolve util/monomers_library paths independent of current working dir
# ------------------------------------------------------------------------------
import util.Monomers as MONOMER

_IPA_REPO_DIR = Path(__file__).resolve().parent
_orig_GetMonomerPath = getattr(MONOMER,"GetMonomerPath",None)
_orig_MonomerExists = getattr(MONOMER,"MonomerExistsInLibrary",None)

def _ipa_resolve_monomer_path(monomer_name: str) -> str:
    """Returns absolute path to CCD CIF in util/monomers_library for monomer_name."""
    #1) ask original implementation first (often returns util/monomers_library/..)
    if callable(_orig_GetMonomerPath):
        p = Path(_orig_GetMonomerPath(monomer_name))
        if p.is_absolute():
            return str(p)
        cand = (_IPA_REPO_DIR / p).resolve()
        if cand.exists():
            return str(cand)
    #2) deterministic fallback (case-insensitive folder letter)
    first = (monomer_name[:1] or "x")
    for letter in (first.lower(), first.upper()):
        cand = (_IPA_REPO_DIR / "util" / "monomers_library" / letter / f"{monomer_name}.cif").resolve()
        if cand.exists():
            return str(cand)
    #3) last resort: return the "repo-rooted" expected path (even if missing)
    return str((_IPA_REPO_DIR / "util" / "monomers_library" / first.lower() / f"{monomer_name}.cif").resolve())

def _ipa_monomer_exists(monomer_name: str) -> bool:
    try:
        return Path(_ipa_resolve_monomer_path(monomer_name)).is_file()
    except Exception:
        if callable(_orig_MonomerExists):
            return bool(_orig_MonomerExists(monomer_name))
        return False

#patch module functions (for other util code that imports util.Monomers.*)
MONOMER.GetMonomerPath = _ipa_resolve_monomer_path
MONOMER.MonomerExistsInLibrary = _ipa_monomer_exists

#patch local names imported via "from util.Monomers import *"
GetMonomerPath = _ipa_resolve_monomer_path
MonomerExistsInLibrary = _ipa_monomer_exists

# ------------------------------------------------------------------------------
# versioning info
# ------------------------------------------------------------------------------
def _extract_ipa_version_from_docstring(doc: str) -> str:
    """Extract the first version token 'vYYMMDD' (or 'v########') from the module docstring.
       Assumes the changelog lines start with something like: 'v251227 - ...'"""
    if not doc:
        return "vUNKNOWN"
    for line in doc.splitlines():
        line = line.strip()
        #match 'v' + digits at beginning of the line
        m = re.match(r"^(v\d+)\b", line)
        if m:
            return m.group(1)
    return "vUNKNOWN"

IPA_VERSION = _extract_ipa_version_from_docstring(__doc__)

# ------------------------------------------------------------------------------
# IPA supporting data classes
# ------------------------------------------------------------------------------
class ProtonationMode(Enum):
    """Enumeration for protonation algorithms."""
    RDKIT = "rdkit"
    OBABEL = "obabel"
    VSEPR = "vsepr"
    CCD = "ccd"
    ORIGINAL = "original"

class LigandTopologyType(Enum):
    """Enumeration for ligand topology types."""
    NONE = auto()
    SDF_FILE = auto()
    SMI_FILE = auto()
    SMILES_STRING = auto()

class CalculationType(Enum):
    """Types of QM calculations."""
    IN_POCKET = "1"
    FULL_OPT = "0"
    BASE_ENERGY = "-1"

@dataclass
class ClashAnalysis:
    """Data object allocating EDDA results for QM clash analysis.
       Notes
       -----
       - short_contact_indices is the per-atom EDDA local score list (ligand atoms first, then residue atoms), matching atom_labels.
       - split_index gives the index where residue atoms start in atom_labels / short_contact_indices."""
    repulsion_energy: List = field(default_factory = list)
    EDDA_binding_energy: List = field(default_factory = list)
    heavy_atom_list: List = field(default_factory = list)
    repulsion_contacts: List = field(default_factory = list)
    interaction_contacts: List = field(default_factory = list)
    short_contact_indices: List = field(default_factory = list)
    #EDDA local-score labeling helpers (for downstream outputs)
    split_index: int = 0
    atom_labels: List[str] = field(default_factory = list)
    residue_id: str = ""
    
@dataclass(frozen = True)
class QMParameters:
    """Immutable QM calculation parameters."""
    alpha: float = 1.0
    kappa: float = 0.025
    ignore_hydrogens: bool = True
    electronic_temperature: float = 300.0
    use_solvation: bool = True
    solvent_name: str = "water"
    calculate_density: bool = False
    energy_threshold: float = 5.0e-6
    gradient_threshold: float = 5.0e-3

@dataclass
class AnalysisConfig:
    """Configuration class for In-Pocket Analysis parameters."""
    #file paths
    pdb_file: Path
    ligand_file: Optional[Path] = None
    output_folder: Path = field(default_factory=lambda: Path.cwd() / "analysis_output")
    geometry_folder: str = "geometries"
    output_folder_qm: str = "output"
    keep_png: bool = False  #keep PNG previews
    keep_png_hq: bool = False
    #ligand parameters
    charge: int = 0
    user_defined_charge: bool = False
    force_user_charge: bool = False
    active_residue: Union[int, str] = "NAN"
    active_residue_chain: str = "A"
    active_residue_name: str = "0"
    active_residue_is_like_aa: bool = False
    #analysis parameters
    critical_distance: float = 5.0
    cap_residues: bool = True
    remove_cap_effect: bool = True
    ph: float = 7.0
    ignore_hydrogens: bool = True
    #input sanitisation
    #    default behaviour: remove H atoms already present in the input PDB/mmCIF (RDKit often misbehaves)
    keep_input_h: bool = False                #if True, do NOT strip input H
    keep_external_caps: bool = False          #if True, do NOT strip external Maestro-like caps from input PDB
    #QM parameters
    qm_params: QMParameters = field(default_factory = QMParameters)
    do_in_pocket: bool = True
    evaluate_only: bool = False
    skip_deformation_energy: bool = False
    use_ip_structure: bool = False
    use_ip_residues: bool = False
    #processing options
    skip_water: bool = True
    flip_residues: bool = True
    flip_ASN: bool = True
    flip_GLN: bool = True
    flip_HIS: bool = True
    flip_SER: bool = True
    flip_THR: bool = True
    flip_TYR: bool = True
    flip_NA: bool = True
    edda: bool = False
    debug: bool = False  #enable very verbose logging (per-residue/per-cap)
    exclude_rdkit_bad_structures: bool = True        #exclude bad tautomers and protonation states originates by RDKit
    #protonation and tautomerisation
    protonation_mode: ProtonationMode = ProtonationMode.RDKIT
    use_ligand_topology: LigandTopologyType = LigandTopologyType.NONE
    ligand_topology: str = ""
    tautomerise: bool = False
    max_tautomers: int = 10
    rmsd_tautomer_threshold: float = 0.05
    check_residue_structures: bool = False
    #file format options
    mmcif_file: bool = False
    ccd_file: str = ""
    covalent_threshold: float = 1.45            #increased by 0.05 to solve issues with PDB 7NRE
    altloc: str = "A"
    check_altloc: bool = False
    #internal lists
    neighbor_residues: List = field(default_factory = list)
    #electron density options
    keep_cub: bool = False
    #performance options
    parallel_processing: bool = False
    max_workers: Optional[int] = None
    #best-ligand post-processing
    #best_mode: enable generation of best_interaction_energies.csv + best_lig.sdf
    #best_debug: if True, keep all original outputs; if False, remove everything
    #            except the two "best" summary files.
    best_mode: bool = False
    best_debug: bool = False
    #best unique residue list (only for the best ligand):
    #write best_list_residues_unique.csv with one entry per residue (best Eint variant)
    best_unique: bool = False

    def __post_init__(self):
        """Post-initialization processing."""
        #convert paths to Path objects
        if isinstance(self.pdb_file,str):
            self.pdb_file = Path(self.pdb_file)
        if isinstance(self.output_folder,str):
            self.output_folder = Path(self.output_folder)
        #auto-detect file format
        if self.pdb_file.suffix.lower() in [".cif",".mmcif"]:
            self.mmcif_file = True

@dataclass
class LigandState:
    """Represents a single ligand protonation/tautomer state."""
    file_path: Path
    charge: int
    ligand_index: int
    smiles: str = ""
    energies: Dict[str, Optional[float]] = field(default_factory = dict)
    converged: Dict[str, bool] = field(default_factory = dict)
    rmsd: Dict[str, Optional[float]] = field(default_factory = dict)
    acid_base_correction: float = 0.0

    @property
    def energy_ip(self) -> Optional[float]:
        """Get IP energy for compatibility."""
        return self.energies.get("ip")

    @property
    def energy_opt(self) -> Optional[float]:
        """Get OPT energy for compatibility."""
        return self.energies.get("opt")

    @property
    def energy_base(self) -> Optional[float]:
        """Get base energy for compatibility."""
        return self.energies.get("base")

    @property
    def deformation_energy(self) -> Optional[float]:
        """Calculate deformation energy in kcal/mol."""
        ip = self.energies.get("ip")
        opt = self.energies.get("opt")
        if ip is not None and opt is not None:
            return (ip - opt)*UC.Hartree2Kcalmol
        return None

    @property
    def distortion_energy(self) -> Optional[float]:
        """Calculate distortion energy in kcal/mol."""
        base = self.energies.get("base")
        ip = self.energies.get("ip")
        if base is not None and ip is not None:
            return (base - ip)*UC.Hartree2Kcalmol
        return None

@dataclass
class BindingResult:
    """Results from binding energy calculation."""
    residue_name: str
    residue_state: str
    binding_energy: float
    cap_correction: float = 0.0
    converged: bool = True
    components: Dict[str,float] = field(default_factory = dict)
    clash: Optional[ClashAnalysis] = None

class InPocketAnalyzerError(Exception):
    """Custom exception for InPocketAnalyzer errors."""
    pass

class FileManager:
    """Manages file operations and directory structure."""

    def __init__(self, output_folder: Path, geometry_folder: str, output_folder_qm: str):
        """Parameters
            ----------
            output_folder : Path
                Root output directory (self.config.output_folder).
            geometry_folder : str
                Name of the geometry subfolder (self.config.geometry_folder,
                typically "geometries").
            output_folder_qm : str
                Name of the QM-output subfolder (self.config.output_folder_qm,
                typically "output")."""
        self.output_folder = output_folder
        self.geometry_folder = geometry_folder
        self.output_folder_qm = output_folder_qm
        root_geom = output_folder / geometry_folder
        root_qm = output_folder / output_folder_qm
        self.dirs = {
            #root geometry and QM-output folders
            "geometries": root_geom,
            "output": root_qm,
            #geometry subfolders (match serial layout)
            "start_xyz": root_geom / "start_xyz",
            "start_sdf": root_geom / "start_sdf",
            "start_pdb": root_geom / "start_pdb",
            "smiles": root_geom / "smiles",
            "ip": root_geom / "ip",
            "opt": root_geom / "opt"
        }
        self._setup_directories()

    def _setup_directories(self):
        """Create necessary directories (geometry + QM output)."""
        for directory in self.dirs.values():
            directory.mkdir(parents = True,exist_ok = True)

    def get_path(self, category: str, filename: str) -> Path:
        """Get full path for a file in a specific category."""
        if category not in self.dirs:
            raise ValueError(f"Unknown category: {category}")
        return self.dirs[category] / filename

    def move_file(
        self, source: Path, category: str, new_name: Optional[str] = None
    ) -> Optional[Path]:
        """Move file to appropriate directory."""
        if not source.exists():
            return None
        dest_name = new_name or source.name
        dest = self.get_path(category,dest_name)
        source.rename(dest)
        return dest

class QMCalculator:
    """Handles quantum mechanical calculations."""

    def __init__(self, qm_params: QMParameters, file_manager: FileManager):
        self.qm_params = qm_params
        self.file_manager = file_manager
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        #RDKit residue protonation summary counters (kept lightweight)
        self._rdkit_res_calls = 0
        self._rdkit_res_fail = 0

    def calculate_energy(self, input_file: str, charge: int, calc_type: CalculationType, output_base: str) -> Tuple[
        Optional[float], Optional[float], Optional[float], Optional[list], Optional[list], Optional[list], bool]:
        #run QM calculation and return full ParseOutput tuple:
        #    (energy, RMSD, totalpol, atmlist, charges, polarisabilities, converged).
        try:
            out,_ = QM.RunIP(input_file,
                             str(charge),
                             str(self.qm_params.alpha),
                             str(self.qm_params.kappa),
                             str(int(self.qm_params.ignore_hydrogens)),
                             calc_type.value,
                             str(self.qm_params.electronic_temperature),
                             str(int(self.qm_params.use_solvation)),
                             self.qm_params.solvent_name,
                             output_base,
                             str(self.qm_params.energy_threshold),
                             str(self.qm_params.gradient_threshold))
            (energy,rmsd,totalpol,atmlist,charges,polarisabilities,converged) = QM.ParseOutput(out)
            return (energy,rmsd,totalpol,atmlist,charges,polarisabilities,converged)
        except Exception as e:
            self.logger.error(f"QM calculation failed for {input_file}: {e}")
            return None,None,None,None,None,None,False

class ResidueProcessor:
    """Handles residue processing and protonation."""

    def __init__(self, config: AnalysisConfig, file_manager: FileManager):
        self.config = config
        self.file_manager = file_manager
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        #RDKit residue protonation summary counters (lightweight, used for INFO summary)
        self._rdkit_res_calls = 0
        self._rdkit_res_fail = 0

    def find_phosphorous(self, atoms):
        #find the phosphorous in the system
        p_type = ""
        p_idx = -1
        for idatm in range(len(atoms)):
            if atoms[idatm] == "P":
                p_type = atoms[idatm]
                p_idx = idatm
                break
        return p_type,p_idx

    def modify_phosphate(self,
                         residue: RES.ResidueData,
                         action: str,
                         idx_modify: int,
                         coordinates: list[float],
                         old_e: str,
                         new_e: str):
        #function to apply simple modifications to phosphates, either by replacing or removing atoms
        new_residue = copy.deepcopy(residue)
        modified = 0
        if action.lower() == "replace":
            if (old_e == "H") and (new_e == "O"): 
                #in this scenario, the phosphate was missing a terminal oxygen, so protonation algorithms added a proton
                #we fix this by replacing the proton with an oxygen and shift the distance
                new_residue.atoms.append("O")
                new_residue.atom_types.append("OP3")
                new_residue.geometry.append([coordinates[0],coordinates[1],coordinates[2]])
                new_residue.charge -= 1
                new_residue.protonation_state = str(new_residue.charge)
                new_residue.terminal = "C"
                new_atoms = new_residue.atoms
                new_geom = new_residue.geometry
                for idatm in range(len(new_residue.processed_atoms)):
                    if new_residue.processed_atoms[idatm] != "H": continue
                    elif idatm == idx_modify: 
                        continue
                    new_atoms.append(new_residue.processed_atoms[idatm])
                    new_geom.append(new_residue.processed_geometry[idatm])
                new_residue.processed_geometry = copy.deepcopy(new_geom)
                new_residue.processed_atoms = copy.deepcopy(new_atoms)
                #now replace the protonated SDF
                connectivity,atom_types,new_atoms,new_geom,local_charges = GEOM.ReadSDFWithTopology(new_residue.protonated_file)
                GEOM.WriteSDF(atom_types,new_residue.processed_atoms,new_residue.processed_geometry,connectivity,local_charges,new_residue.protonated_file,False)
                modified = 1
        elif action.lower() == "remove":
            if (old_e == "H") and (new_e == ""):
                #in this scenario we remove a proton from a phosphate
                #assumes that proton goes to OP2
                connectivity,atom_types,new_atoms,new_geom,local_charges = GEOM.ReadSDFWithTopology(new_residue.protonated_file)
                #find the phosphorous atom
                p_type,p_idx = self.find_phosphorous(new_atoms)
                p_type = atom_types[p_idx]
                if p_type + "blahblah" != "blahblah":
                    idx_H = -1
                    #now find its neighbour with 2 bonding partners, one of which, a proton
                    for icnn in connectivity:
                        if (icnn[0] == p_type) or (icnn[1] == p_type):
                            o_type = icnn[1] if icnn[0] == p_type else icnn[0]
                            the_other_type = ""
                            n_bond_partners = 0
                            for jcnn in connectivity:
                                if (jcnn[0] == o_type) or (jcnn[1] == o_type): 
                                    n_bond_partners += 1
                                    if (jcnn[0] == o_type) and (jcnn[1] != p_type): the_other_type = jcnn[1]
                                    elif (jcnn[1] == o_type) and (jcnn[0] != p_type): the_other_type = jcnn[0]
                            if n_bond_partners == 2 and the_other_type.startswith("H"): 
                                for idatm in range(len(atom_types)):
                                    if atom_types[idatm] == the_other_type:
                                        idx_H = idatm
                                        break
                                if idx_H >= 0: break
                    if idx_H >= 0:
                        #now we can safely remove the proton
                        modified = 1
                        #no need to change: atoms, geometry, atom_types, name, number, chain_id, residue_type, pdb_lines, cap_geometry, cap_atoms, terminal, residue_id, ip_file
                        new_residue.charge -= 1
                        new_residue.protonation_state = str(new_residue.charge)
                        #reset this to contain the new data
                        new_atoms = []
                        new_geom = []
                        new_types = []
                        for iadtm in range(len(new_residue.processed_atoms)):
                            if iadtm == idx_H: continue
                            new_atoms.append(new_residue.processed_atoms[iadtm])
                            new_geom.append(new_residue.processed_geometry[iadtm])
                            new_types.append(atom_types[iadtm])
                        #now the new connectivity:
                        new_connectivity = []
                        for icnn in connectivity:
                            if (icnn[0] == atom_types[idx_H]) or (icnn[1] == atom_types[idx_H]): continue
                            new_connectivity.append(icnn)
                        #update geometry containers
                        new_residue.processed_atoms = copy.deepcopy(new_atoms)
                        new_residue.processed_geometry = copy.deepcopy(new_geom)
                        #modify the sdf
                        GEOM.WriteSDF(new_types,new_atoms,new_geom,new_connectivity,local_charges,new_residue.protonated_file,False)
        return new_residue,modified

    def protonate_phosphate(self,
                         residue: RES.ResidueData,
                         action: str,
                         start_idx: int):
        #function to modify phosphates, either by replacing or adding new atoms
        new_residue = copy.deepcopy(residue)
        connectivity,atom_types,atoms,geom,local_charges = GEOM.ReadSDFWithTopology(new_residue.protonated_file)
        #find the phosphorous in the system
        p_type,p_idx = self.find_phosphorous(atoms)
        p_type = atom_types[p_idx]
        if p_type + "blahblah" != "blahblah":
            idx_H = []
            o_add_H = []
            o_remove_H = []
            o_other = []
            for icnn in connectivity:
                #search for phosphorous neighbours
                if (icnn[0] == p_type) or (icnn[1] == p_type):
                    o_temp = None
                    if icnn[0].startswith("O"): o_temp = icnn[0]
                    elif icnn[1].startswith("O"): o_temp = icnn[1]
                    if o_temp is None: continue
                    #here we have an oxygen neighbour of the phosphate
                    idx_temp = -1
                    #find its index
                    for idatm in range(len(atom_types)):
                        if atom_types[idatm] == o_temp:
                            idx_temp = idatm
                            break
                    if idx_temp == -1: continue
                    nneighbours = 0
                    foundH = False
                    for jcnn in connectivity:
                        if (jcnn[0] == o_temp) or (jcnn[1] == o_temp):
                            nneighbours += 1
                            the_other_atom = jcnn[0]
                            if jcnn[0] == o_temp: the_other_atom = jcnn[1]
                            if the_other_atom.startswith("H"):
                                foundH = True
                                #if we found a proton, get its index
                                for idatm in range(len(atom_types)):
                                    if atom_types[idatm] == the_other_atom:
                                        idx_H.append(idatm)
                                        break
                    if nneighbours == 1: o_add_H.append(idx_temp)
                    elif (nneighbours == 2) and foundH: o_remove_H.append(idx_temp)
                    elif (nneighbours == 2) and (not foundH): o_other.append(idx_temp)
            #print("proton indices",idx_H)
            #print("oxygens to protonate",o_add_H)
            #print("oxygens to deprotonate",o_remove_H)
            #print("oxygens other",o_other)
            if action.lower() == "add":
                #adding protons to all phosphate oxygens that can receive a proton
                counter = start_idx
                outs = []
                for o_idx in o_add_H:
                    #we always take as a reference the oxygen atom connected to the sugar
                    H_coord = NA.place_H_from_dihedral(geom[o_other[0]],geom[p_idx],geom[o_idx],0.96,108.5,180.0)
                    new_residue2 = copy.deepcopy(new_residue)
                    connectivity_new = copy.deepcopy(connectivity)
                    geom_new = copy.deepcopy(geom)
                    atoms_new = copy.deepcopy(atoms)
                    atom_types_new = copy.deepcopy(atom_types)
                    #now we add the proton to the residue
                    new_residue2.geometry.append(H_coord)
                    new_residue2.processed_geometry.append(H_coord)
                    geom_new.append(H_coord)
                    new_residue2.atoms.append("H")
                    new_residue2.processed_atoms.append("H")
                    atoms_new.append("H")
                    H_name = ""
                    for atom_tp in atom_types_new:
                        if atom_tp.startswith("H"): H_name = atom_tp
                    Hindex = int(re.search(r'\d+',H_name).group()) + 1
                    atom_types_new.append("H" + str(Hindex))
                    new_residue2.atom_types.append("H" + str(Hindex))
                    connectivity_new.append([atom_types_new[o_idx],"H" + str(Hindex),"sing"])
                    new_residue2.charge += 1
                    new_residue2.protonation_state = str(new_residue2.charge)
                    new_residue2.residue_id = self.set_residue_id(new_residue2.residue_id,counter)
                    new_residue2.protonated_file = new_residue2.protonated_file.replace(".sdf","_prot_" + str(counter) + ".sdf")
                    for icnn in range(len(connectivity_new)):
                        iconnect = connectivity_new[icnn]
                        #check whether protonated oxygen was involved in double bond with phosphorous and correct
                        if ((iconnect[0] == p_type) and (iconnect[1] == atom_types_new[o_idx])) or ((iconnect[1] == p_type) and (iconnect[0] == atom_types_new[o_idx])):
                            if (iconnect[2] == "doub"):
                                connectivity_new[icnn] = (iconnect[0],iconnect[1],"sing")
                                for jcnn in range(len(connectivity_new)):
                                    #find phosphorous bond to deprotonated oxygen
                                    o_add_double_bond = o_add_H[0]
                                    if (o_add_double_bond == o_idx) or (o_add_double_bond == iconnect[0]) or (o_add_double_bond == iconnect[1]): o_add_double_bond = o_add_H[1]
                                    if (o_add_double_bond == o_idx) or (o_add_double_bond == iconnect[0]) or (o_add_double_bond == iconnect[1]) and (len(o_add_H) > 2): o_add_double_bond = o_add_H[2]
                                    jconnect = connectivity_new[jcnn]
                                    if ((jconnect[0] == p_type) and (jconnect[1] == atom_types_new[o_add_double_bond])) or ((jconnect[1] == p_type) and (jconnect[0] == atom_types_new[o_add_double_bond])):
                                        connectivity_new[jcnn] = (jconnect[0],jconnect[1],"doub")
                                        break
                    GEOM.WriteSDF(atom_types_new,atoms_new,geom_new,connectivity_new,local_charges,new_residue2.protonated_file,False)
                    counter += 1
                    outs.append(new_residue2)
                    if not self.config.flip_NA: continue
                    #now we flip
                    H_coord_2 = NA.place_H_from_dihedral(geom[o_other[0]],geom[p_idx],geom[o_idx],0.96,108.5,0.0)
                    if not np.array_equal(H_coord_2,H_coord):
                        new_residue3 = copy.deepcopy(new_residue2)
                        if np.array_equal(new_residue3.geometry[-1],H_coord):
                            new_residue3.geometry[len(new_residue3.geometry) - 1] = H_coord_2
                        else: 
                            for idatm in range(len(new_residue3.geometry)):
                                if np.array_equal(new_residue3.geometry[idatm],H_coord):
                                    new_residue3.geometry[idatm] = H_coord_2
                        if np.array_equal(new_residue3.processed_geometry[-1],H_coord): 
                            new_residue3.processed_geometry[len(new_residue3.processed_geometry) - 1] = H_coord_2
                        else: 
                            for idatm in range(len(new_residue3.processed_geometry)):
                                if np.array_equal(new_residue3.processed_geometry[idatm],H_coord):
                                    new_residue3.processed_geometry[idatm] = H_coord_2
                        new_residue3.residue_id = self.set_residue_id(new_residue3.residue_id,counter)
                        new_residue3.protonated_file = new_residue3.protonated_file.replace("_prot_" + str(counter - 1) + ".sdf","_prot_" + str(counter) + ".sdf")
                        GEOM.WriteSDF(atom_types_new,atoms_new,new_residue3.processed_geometry,connectivity_new,local_charges,new_residue3.protonated_file,False)
                        counter += 1
                        outs.append(new_residue3)
                return outs

    def compare_connectivity(self, state1, state2, tolerance = 1.0e-4):
        #function comparing the protonation states of two "different" (?) protonated molecules
        matching = True
        _,_,atoms1,geometry1,_ = GEOM.ReadSDFWithTopology(state1.protonated_file)
        _,_,atoms2,geometry2,_ = GEOM.ReadSDFWithTopology(state2.protonated_file)
        if len(atoms1) == len(atoms2):
            used2 = set()
            unique1_atoms = []
            unique1_geom = []
            #compare geometry1 against geometry2
            for i,(a1,c1) in enumerate(zip(atoms1,geometry1)):
                found = False
                for j,(a2,c2) in enumerate(zip(atoms2,geometry2)):
                    if j in used2:
                        continue
                    if a1.startswith("H") and a2.startswith("H"):
                        match = np.allclose(c1,c2,atol = tolerance)
                    elif (not a1.startswith("H")) and (not a2.startswith("H")):
                        match = (a1 == a2) and np.allclose(c1,c2,atol = tolerance)
                    else:
                        match = False
                    if match:
                        used2.add(j)
                        found = True
                        break
                if not found:
                    unique1_atoms.append(a1)
                    unique1_geom.append(c1)
            # remaining unmatched atoms in geometry 2 are unique to geometry 2
            unique2_atoms = []
            unique2_geom = []
            for j,(a2,c2) in enumerate(zip(atoms2,geometry2)):
                if j not in used2:
                    unique2_atoms.append(a2)
                    unique2_geom.append(c2)
            if (len(unique1_geom) != 0) or (len(unique2_geom) != 0): matching = False
        else: matching = False
        return matching
    
    def set_residue_id(self, string: str, value: int):
        #function that increases the residue identifier
        base,num = string.rsplit("_",1)
        return f"{base}_{value}"

    def process_residues(self, residues: List[RES.ResidueData]) -> List[RES.ResidueData]:
        """Process residues with protonation and multiple conformations."""
        #first handle multiple conformations
        separated_residues = self._separate_multiple_conformations(residues)
        #then protonate
        processed_residues = []
        for residue in separated_residues:
            protonated = self._protonate_residue(residue)
            #extra work needed for nucleic acids
            if NA.IsBaseNucleicAcid(residue.name) or NA.IsExtendedNucleicAcid(residue.name):
                new_states = []
                for prot in protonated:
                    phosphate_charge,number_oxygen_neighbours = NA.get_phosphate_charge(prot.protonated_file)
                    add_state = False
                    prot_residue = None
                    prot_state_extra = None
                    add_deprotonated = False
                    #checking states
                    if (number_oxygen_neighbours != 4) and (phosphate_charge is not None):
                        self.logger.warning(f"Phosphate with incomplete octet misses an oxygen.")
                        if phosphate_charge == prot.charge: 
                            #we are safe to replace atoms
                            self.logger.warning(f"Adding oxygen atom.")
                            add_state = True
                            #identify which atoms to swap
                            idx_replace,coordinates_O,old_element,new_element = NA.fix_phosphate_incomplete_octet(prot.protonated_file)
                            #perform the replacement
                            prot_residue_tmp,modified1 = self.modify_phosphate(prot,"replace",idx_replace,coordinates_O,old_element,new_element)
                            #we add the residue without phosphate protons, because it is easier to handle this.
                            prot_residue,_ = self.modify_phosphate(prot_residue_tmp,"remove",-1,[0.0,0.0,0.0],"H","")
                            if modified1: add_deprotonated = True
                    elif (phosphate_charge is None) and (prot.charge == 0): 
                        add_state = True
                        prot_residue = prot
                    elif (phosphate_charge is not None) and (phosphate_charge == prot.charge): 
                        add_state = True
                        prot_residue = prot
                    #only take those states with acid-base in the phosphate (cleaning)
                    if add_state:
                        #add the base state
                        new_states.append(prot_residue)
                        #add a singly neutralised protonation state
                        idx_start = 2
                        prot_residue_2 = self.protonate_phosphate(prot_residue,"add",idx_start)
                        new_states.extend(prot_residue_2)
                        if add_deprotonated: 
                            #since we added previously the -2 protonation state of the phosphate to the CSs, now we add the first proton to generate a common base with any other nucleotide
                            idx_start += len(prot_residue_2)
                            for pres in prot_residue_2:
                                prot_residue_3 = self.protonate_phosphate(pres,"add",idx_start)
                                for pres_3 in prot_residue_3:
                                    repeated = False
                                    for state in new_states:
                                        match = self.compare_connectivity(pres_3,state)
                                        if (pres_3.protonated_file.rsplit("_",1)[0] == state.protonated_file.rsplit("_",1)[0]): match = 0
                                        if match == 1:
                                            self.logger.warning(f"Found repeated pair: {pres_3.residue_id} - {state.residue_id}. Excluding {pres_3.residue_id}")
                                            repeated = True
                                            break
                                    if not repeated:
                                        new_states.append(pres_3)
                                idx_start += len(prot_residue_3)
                protonated = copy.deepcopy(new_states)
            processed_residues.extend(protonated)
            if ((protonated[0].name == "SER") and (self.config.flip_SER)) or ((protonated[0].name == "THR") and (self.config.flip_THR)) or ((protonated[0].name == "TYR") and (self.config.flip_TYR)):
                for cs in protonated:
                    if (cs.charge == -1): continue
                    new_residue = copy.deepcopy(cs)
                    new_residue = RES.flip_hydroxyl_groups(new_residue)
                    if new_residue is not None: processed_residues.append(new_residue)
        #compact summary at INFO (details are in DEBUG)
        if self.config.check_residue_structures: input("Please check the structures of protonated residues and press enter to continue")
        calls = getattr(self,"_rdkit_res_calls", 0)
        if calls:
            self.logger.info("RDKit residue protonation summary: calls = %d failed = %d",
                             calls,
                             getattr(self,"_rdkit_res_fail",0))
        return processed_residues

    def _separate_multiple_conformations(self, residues: List[RES.ResidueData]) -> List[RES.ResidueData]:
        """Handle residues with multiple conformations."""
        separated = []
        for residue in residues:
            #check for repeated atoms
            atom_counts = {}
            for atp in residue.atom_types:
                atom_counts[atp] = atom_counts.get(atp, 0) + 1
            max_repeat = max(atom_counts.values())
            if max_repeat > 1:
                self.logger.info(f"Found {max_repeat} conformations in {residue.residue_id}")
                #separate conformations
                conformations = self._split_conformations(residue,max_repeat)
                separated.extend(conformations)
            else:
                separated.append(residue)
        return separated

    def _split_conformations(self, residue: RES.ResidueData, n_conf: int) -> List[RES.ResidueData]:
        """Split a residue with multiple conformations."""
        conformations = []
        #group atoms by conformation
        conf_data = [{"geom": [],"atp": [],"atoms": [],"pdb": []} for _ in range(n_conf)]
        atom_indices = {}
        for i,atp in enumerate(residue.atom_types):
            if atp not in atom_indices:
                atom_indices[atp] = []
            atom_indices[atp].append(i)
        #distribute atoms to conformations
        for atp,indices in atom_indices.items():
            for conf_idx,atom_idx in enumerate(indices):
                if conf_idx < n_conf:
                    conf_data[conf_idx]["geom"].append(residue.geometry[atom_idx])
                    conf_data[conf_idx]["atp"].append(residue.atom_types[atom_idx])
                    conf_data[conf_idx]["atoms"].append(residue.atoms[atom_idx])
                    conf_data[conf_idx]["pdb"].append(residue.pdb_lines[atom_idx])
        #create separate residue objects
        for conf_idx,data in enumerate(conf_data):
            if data["geom"]:  #only create if has atoms
                new_residue = RES.ResidueData(geometry = data["geom"],
                                              atoms = data["atoms"],
                                              atom_types=data["atp"],
                                              name = residue.name,
                                              number = residue.number,
                                              chain_id = residue.chain_id,
                                              residue_type = residue.residue_type,
                                              pdb_lines = data["pdb"],
                                              terminal = residue.terminal,
                                              residue_id = f"{residue.residue_id}_conf{conf_idx + 1}",
                                              cap_geometry = residue.cap_geometry.copy(),
                                              cap_atoms = residue.cap_atoms.copy())
                conformations.append(new_residue)
        return conformations

    def _protonate_residue(self, residue: RES.ResidueData) -> List[RES.ResidueData]:
        """Protonate a single residue and return all protonation states."""
        #write residue files
        res_file = self.file_manager.get_path("geometries",residue.residue_id)
        GEOM.Write2XYZ(residue.atoms,residue.geometry,str(res_file))
        if not self.config.mmcif_file:
            pdb_lines = residue.pdb_lines
            #sanitise structure?
            if not self.config.keep_input_h:
                pdb_lines,n_removed = SANITY.InputSanitizer.strip_h_from_pdb_lines(pdb_lines)
                if n_removed:
                    #keep this quiet by default; detailed work remains for --debug
                    self.logger.debug("Stripped %d input H atoms from residue container: %s (use --keep-input-h to disable).",n_removed,residue.residue_id)
            GEOM.Write2PDB(pdb_lines,str(res_file))
        else:
            self._write_mmcif_residue(residue,str(res_file))
        #single-atom handling:
        # - Ions/metals: keep the fast-path (no protonation backend; charge inferred from mapping)
        # - Everything else (e.g. waters stripped to a single O): DO go through the protonation backend
        #   so RDKit/OBabel/VSEPR can rebuild missing hydrogens. The "single-atom fast-path" is then
        #   applied later *only if* the protonated structure is still mono-atomic.        
        n_atoms = len(residue.atoms)
        if n_atoms == 0: 
            self.logger.warning("Residue %s detected with no atoms",residue.residue_id)
        elif n_atoms == 1:
            #IMPORTANT: charge must be inferred here for ions/metals, otherwise QM runs become neutral
            #only use the single-atom fast-path for true ions and metals, the latter included in the ion pack
            is_ion,charge = ION.IsIon(residue.name)
            #non-ions (e.g. HOH after input-H stripping) must still go through protonation.
            if is_ion:
                #ensure the XYZ uses a valid element symbol as atom label for QM
                # since some PDB parsers can leave atom names like "MN1", then residue.name is the safest label
                GEOM.Write2XYZ([str(residue.name).strip()],[residue.geometry[0]],str(res_file))
                if int(charge) == 0:
                    self.logger.warning("Ion-like single-atom residue %s detected but charge resolved to 0. "
                                        "Check util/Ions.py (and/or BaseAminoAcids.IsMetal) mapping.",residue.residue_id)
                prot_file = str(res_file) + ".xyz"
                return [self._create_protonated_residue(residue,prot_file,int(charge),0)]
        #check if ion
        is_ion,charge = ION.IsIon(residue.name)
        print("CHECK HERE, THIS IS SUPPOSED TO DO LITERALLY NOTHING")
        if is_ion:
            return self._process_ion(residue,charge,str(res_file))
        #generate SDF
        connectivity = []
        terminal_type = "ipa_termini"
        if residue.terminal == "N": terminal_type = "0"
        elif residue.terminal == "C": terminal_type = "ipa_termini_C-term"
        residue_name = residue.name
        if "_flip" in residue.name: residue_name = residue.name.replace("_flip","")
        #amino acids and PTMs
        if AA.IsBaseAminoAcid(residue_name): connectivity = AA.BaseAminoAcidConnectivity(residue_name,terminal_type)
        elif AA.IsExtendedAminoAcid(residue_name): connectivity = AA.ExtendedAminoAcidConnectivity(residue_name,terminal_type)
        elif AA.IsPostTranslationalModification(residue_name): connectivity = AA.PTMConnectivity(residue_name,terminal_type)
        #nucleic acids
        elif NA.IsBaseNucleicAcid(residue_name): connectivity = NA.BaseNucleicAcidConnectivity(residue_name,terminal_type)
        elif NA.IsExtendedNucleicAcid(residue_name): connectivity = NA.ExtendedNucleicAcidConnectivity(residue_name,terminal_type)
        elif MonomerExistsInLibrary(residue_name): 
            path2monomerfile = GetMonomerPath(residue_name)
            connectivity = MONOMER.ExtractCCDconnectivityGEMMI(path2monomerfile)
        local_charges = AA.GetLocalCharges(residue_name,residue.atom_types)
        GEOM.WriteSDF(residue.atom_types,residue.atoms,residue.geometry,connectivity,local_charges,str(res_file))
        #protonate based on mode
        if self.config.protonation_mode == ProtonationMode.RDKIT:
            return self._rdkit_protonate_residue(residue,str(res_file))
        elif self.config.protonation_mode == ProtonationMode.OBABEL:
            return self._obabel_protonate_residue(residue,str(res_file))
        elif self.config.protonation_mode == ProtonationMode.ORIGINAL:
            #in original/native mode we keep the residue as provided, preserve any existing input hydrogens, and only use the protein/VSEPR path
            #to complete missing cap/terminal hydrogens introduced by capping
            return self._vsepr_protonate_residue(residue,str(res_file))
        else:
            return self._vsepr_protonate_residue(residue,str(res_file))

    def _rdkit_protonate_residue(self, residue: RES.ResidueData, res_file: str) -> List[RES.ResidueData]:
        """Protonate residue using RDKit."""
        self._ensure_worker_logger()
        #use residue-specific pH
        local_ph = self._get_residue_ph(residue.name)
        file_in = res_file
        file_for_rdkit = file_in
        #file_for_rdkit = self._prepare_rdkit_residue(file_in,residue.name)
        #debug: confirm what we pass to RDKit/Protoro and how many H atoms remain
        h_count = 0
        atoms_sdf,_ = GEOM.ReadSDF(file_for_rdkit)
        for atom in atoms_sdf:
            if atom == "H": h_count += 1
        #counters + verbose line only in DEBUG
        self._rdkit_res_calls += 1
        self.logger.debug("RDKit residue protonation input: %s -> %s | residue=%s | pH=%.2f | H_atoms=%d | keep_input_h=%s",
                          file_in,
                          file_for_rdkit,
                          residue.residue_id,
                          local_ph,
                          h_count,
                          getattr(self.config,"keep_input_h",False))
        #required for --no cap residues, otherwise some protonation states can be missed due to the lack of capping groups (GLU/ASP)
        local_ph_spread = 0.5 
        if ((residue.name == "GLU") or (residue.name == "ASP")) and (not self.config.cap_residues): local_ph_spread = 3.0
        try:
            res_files,res_charges,_ = Protoro.RDKitProtonation(file_for_rdkit + ".sdf",local_ph,"none",pH_spread = local_ph_spread)
        except Exception as e:
            self.logger.error(f"RDKit protonation failed for {residue.residue_id}: {e}")
            self.logger.error(f"Attempting VSEPR protonation.")
            self._rdkit_res_fail += 1
            local_charge = VSEPR.VSEPRSDFMOL2Protonation(file_for_rdkit + ".sdf",file_for_rdkit + "_H")
            res_charges = [local_charge]
            res_files = [file_for_rdkit + "_H.sdf"]
        protonated_residues = []
        for prot_idx,(prot_file,charge) in enumerate(zip(res_files,res_charges)):
            new_residue = self._create_protonated_residue(residue,prot_file,charge,prot_idx)
            protonated_residues.append(new_residue)
        return protonated_residues

    def _ensure_worker_logger(self) -> None:
        """Ensure logger has at least one handler.
           This is important when methods are executed in multiprocessing workers,
           where the parent logging config may not be inherited."""
        import logging
        try:
            logger = getattr(self,"logger",None)
            if logger is None: return
            if getattr(logger,"handlers",None) and len(logger.handlers) > 0: return
            #minimal safe handler so INFO/WARNING/ERROR show up in worker stdout
            handler = logging.StreamHandler()
            formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            logger.setLevel(logging.INFO)
            logger.propagate = False
        except Exception:
            #never break the workflow due to logging
            return

    def _get_residue_ph(self, residue_name: str) -> float:
        """Get appropriate pH for specific residue types."""
        if residue_name == "HIS": return 6.0
        elif residue_name in ["GLU","ASP"]: return 4.0
        else: return 7.0

    def _create_protonated_residue(self, original: RES.ResidueData, prot_file: str, charge: int, prot_idx: int) -> RES.ResidueData:
        """Create a new residue object for a protonation state."""
        pstate = FORMAT.Format_Charge_State(charge)
        new_residue = RES.ResidueData(geometry = original.geometry.copy(),
                                      atoms = original.atoms.copy(),
                                      atom_types = original.atom_types.copy(),
                                      name = original.name,
                                      number = original.number,
                                      chain_id = original.chain_id,
                                      residue_type = original.residue_type,
                                      pdb_lines = original.pdb_lines.copy(),
                                      cap_geometry = original.cap_geometry.copy(),
                                      cap_atoms = original.cap_atoms.copy(),
                                      charge = charge,
                                      protonation_state = pstate,
                                      terminal = original.terminal,
                                      residue_id = f"{original.name}_{original.chain_id}{original.number}_{prot_idx + 1}",
                                      protonated_file = prot_file)
        #load processed geometry
        try:
            new_residue.processed_atoms,new_residue.processed_geometry = GEOM.ReadGeometry(prot_file)
        except Exception as e:
            self.logger.warning(f"Could not read processed geometry: {e}")
            new_residue.processed_atoms = original.atoms.copy()
            new_residue.processed_geometry = original.geometry.copy()
        return new_residue

    def _create_fallback_residue(self, original: RES.ResidueData, res_file: str) -> RES.ResidueData:
        """Create fallback residue when protonation fails."""
        fallback = self._create_protonated_residue(original,res_file + ".pdb",0,0)
        fallback.residue_id = original.residue_id
        return fallback

    def _write_mmcif_residue(self, residue: RES.ResidueData, res_file: str):
        """Write residue in mmCIF format."""
        residue_nam = [residue.name]*len(residue.atoms)
        residue_idx = [residue.number]*len(residue.atoms)
        residue_typ = [residue.residue_type]*len(residue.atoms)
        chain_ident = [residue.chain_id]*len(residue.atoms)
        GEOM.WriteGemmiPDB(residue.geometry,residue.atoms,residue.atom_types,residue_nam,residue_idx,chain_ident,residue_typ,res_file)

    def _process_ion(self, residue: RES.ResidueData, charge: int, res_file: str) -> List[RES.ResidueData]:
        """Process ion residue."""
        try:
            res_file_xyz = GEOM.PDB2XYZ(res_file + ".pdb")
            return [self._create_protonated_residue(residue,res_file_xyz,charge,0)]
        except Exception as e:
            self.logger.error(f"Ion processing failed: {e}")
            return [self._create_fallback_residue(residue,res_file)]

    def _obabel_protonate_residue(self, residue: RES.ResidueData, res_file: str) -> List[RES.ResidueData]:
        """Protonate residue using OpenBabel."""
        try:
            charge = GEOM.ObabelProtonation(res_file + ".sdf",self.config.ph)
            GEOM.MOL2toXYZ(res_file + "H.mol2",residue.name)
            #adjust charges for specific residues
            if residue.name == "ARG" and charge == 0: charge = 1
            elif residue.name in ["ASP","GLU"] and charge == 0: charge = -1
            return [self._create_protonated_residue(residue, res_file + "H.xyz",charge,0)]
        except Exception as e:
            self.logger.error(f"OpenBabel protonation failed: {e}")
            return [self._create_fallback_residue(residue,res_file)]

    def _vsepr_protonate_residue(self, residue: RES.ResidueData, res_file: str) -> List[RES.ResidueData]:
        """Protonate residue using VSEPR."""
        try:
            is_ion,charge = ION.IsIon(residue.name)
            print(residue.name,is_ion)
            if not is_ion:
                if AA.IsBaseAminoAcid(residue.name) or AA.IsExtendedAminoAcid(residue.name)or AA.IsPostTranslationalModification(residue.name): 
                    print("base")
                    charge = VSEPR.ProteinProtonation(res_file + ".pdb",res_file + "H")
                    prot_file = res_file + "H.xyz"
                elif residue.name.upper() in ["HOH","SOL","WAT","H2O","TP3","TIP"]:
                    print("water")
                    charge = VSEPR.VSEPRProtonation(res_file + ".xyz",res_file)
                    prot_file = res_file + ".xyz"
                elif os.path.exists(res_file + ".sdf"):
                    print("sdf")
                    charge = VSEPR.VSEPRSDFMOL2Protonation(res_file + ".sdf", res_file + "H")
                    prot_file = res_file + "H.xyz"
                else:
                    print("xyz")
                    self.logger.warning("Using the generic VSEPR protonation.",residue.name,e)
                    charge = VSEPR.VSEPRProtonation(res_file + ".pdb",res_file)
                    prot_file = res_file + ".xyz"
            else:
                prot_file = GEOM.PDB2XYZ(res_file + ".pdb")
            return [self._create_protonated_residue(residue, prot_file,charge,0)]
        except Exception as e:
            self.logger.error(f"VSEPR protonation failed: {e}")
            return [self._create_fallback_residue(residue,res_file)]

class BindingCalculator:
    """Handles binding energy calculations."""

    def __init__(self, config: AnalysisConfig, qm_calculator: QMCalculator, file_manager: FileManager):
        self.config = config
        self.qm_calc = qm_calculator
        self.file_manager = file_manager
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        self.Eh2kcalmol = UC.Hartree2Kcalmol
        self.clash_analysis_list: List[List[ClashAnalysis]] = []

    def _density_prefix(self, base_tag: str) -> str:
        """Return the output prefix (no extension) for density files.
           We ALWAYS write .cub under <output>/cub_files to avoid polluting geometries/.
           base_tag is typically:
             - '<...>/geometries/TYR_B253_ip'
             - '<...>/geometries/<pair>_ip'"""
        cub_dir = self.config.output_folder / "cub_files"
        cub_dir.mkdir(parents = True,exist_ok = True)
        name = Path(base_tag).name
        return str(cub_dir / f"{name}_dens")

    def calculate_binding_energies(self, ligand_states: List[LigandState], residues: List[RES.ResidueData]) -> List[List[BindingResult]]:
        """Calculate all binding energies."""
        #always reset to exactly one list per ligand state (avoid accumulation across runs)
        self.clash_analysis_list = [[] for _ in ligand_states]
        #calculate individual energies
        residue_energies = self._calculate_residue_energies(residues)
        #calculate pair energies
        pair_energies = self._calculate_pair_energies(ligand_states,residues)
        #any missing pair energy is a failure.
        #if --edda was requested, any EDDA failure is also a failure.
        self.pair_failures = 0
        self.edda_failures = 0
        for lig_idx in range(len(ligand_states)):
            for res_idx in range(len(residues)):
                e = pair_energies[lig_idx][res_idx]
                if e is None:
                    self.pair_failures += 1
                    continue
                if self.config.edda:
                    clash = None
                    try: clash = self.clash_analysis_list[lig_idx][res_idx]
                    except Exception: clash = None
                    #EDDA requested but missing/failed clash -> count as EDDA failure
                    if clash is None or getattr(clash,"edda_failed",False): self.edda_failures += 1
        self.analysis_failed = (self.pair_failures > 0) or (self.config.edda and self.edda_failures > 0)
        #calculate cap corrections if needed
        cap_corrections = {}
        if self.config.cap_residues and self.config.remove_cap_effect:
            cap_corrections = self._calculate_cap_corrections(ligand_states,residues)
        #combine results
        binding_results = []
        for lig_idx,ligand in enumerate(ligand_states):
            lig_results = []
            for res_idx,residue in enumerate(residues):
                cap_corr = cap_corrections.get(self._get_base_residue_id(residue),0.0)
                cap_corr = (cap_corr.get(lig_idx,0.0)
                            if isinstance(cap_corr,dict)
                            else cap_corr)
                result = self._calculate_single_binding(ligand,residue,residue_energies[res_idx],pair_energies[lig_idx][res_idx],cap_corr)
                if result:
                    #attach EDDA clash object aligned by residue index (NOT by filtered result index)
                    if self.config.edda:
                        try:
                            result.clash = self.clash_analysis_list[lig_idx][res_idx]
                        except Exception:
                            result.clash = None
                    lig_results.append(result)
            binding_results.append(lig_results)
        return binding_results

    def _calculate_residue_energies(self, residues: List[RES.ResidueData]) -> List[Optional[float]]:
        """Calculate energies for all residues."""
        self.logger.info(f"Calculating energies for {len(residues)} residues")
        if self.config.parallel_processing:
            return self._calculate_residue_energies_parallel(residues)
        else:
            return self._calculate_residue_energies_sequential(residues)

    def _calculate_residue_energies_sequential(self, residues: List[RES.ResidueData]) -> List[Optional[float]]:
        """Calculate residue energies sequentially."""
        energies = []
        for residue in residues:
            energy = self._calculate_single_residue_energy(residue)
            energies.append(energy)
        return energies

    def _calculate_residue_energies_parallel(self, residues: List[RES.ResidueData]) -> List[Optional[float]]:
        """Calculate residue energies in parallel."""
        energies = [None]*len(residues)
        with ProcessPoolExecutor(max_workers = self.config.max_workers) as executor:
            #create futures with their indices
            futures = [
                executor.submit(self._calculate_single_residue_energy,res)
                for res in residues
            ]
            #wait for all futures to complete and collect results in order
            for idx,future in enumerate(futures):
                try:
                    energies[idx] = future.result()
                except Exception as e:
                    self.logger.error(f"Failed to calculate residue energy: {e}")
                    energies[idx] = None
        #ensure residue.ip_file is populated in the parent process
        for res in residues:
            if res.protonated_file and (not res.ip_file):
                res_ip = res.protonated_file.replace(".xyz","_ip").replace(".sdf","_ip")
                res.ip_file = f"{res_ip}.xyz"
        return energies

    def _calculate_single_residue_energy(self, residue: RES.ResidueData) -> Optional[float]:
        """Calculate energy for a single residue and write <RESID>_ip.out."""
        if not residue.protonated_file:
            return None
        #use IN_POCKET calculation for residues
        res_ip = residue.protonated_file.replace(".xyz","_ip").replace(".sdf","_ip")
        #check whether to run proper refinement or whether we just evaluate
        #single-atom residues have no internal degrees of freedom -> skip refinement
        atoms = getattr(residue,"processed_atoms",None) or residue.atoms
        n_atoms = len(atoms)
        if n_atoms == 0:
            self.logger.warning("Residue %s has 0 atoms (unexpected)",residue.residue_id)
        single_atom = (n_atoms == 1)
        calc_mode = CalculationType.IN_POCKET
        if self.config.evaluate_only or single_atom:
            calc_mode = CalculationType.BASE_ENERGY
        (energy,rmsd,totalpol,atmlist,charges,polarisabilities,converged) = self.qm_calc.calculate_energy(residue.protonated_file,
                                                                                                          residue.charge,
                                                                                                          calc_mode,
                                                                                                          res_ip)
        if energy is None:
            return None
        #ip geometry: for single-atom residues, just copy the starting XYZ to *_ip.xyz
        if single_atom:
            try: shutil.copyfile(residue.protonated_file,res_ip + ".xyz")
            except Exception as e: self.logger.debug("Failed to copy single-atom residue ip geometry (%s): %s",residue.residue_id,e)
        else:
            #recentering ip geometry
            GEOM.ReCenterAndWrite(residue.protonated_file,res_ip + ".xyz",res_ip)
        if self.config.qm_params.calculate_density:
            density_file = self._density_prefix(res_ip)
            QM.GetDensity(res_ip + ".xyz",
                       str(residue.charge),
                       str(self.config.qm_params.electronic_temperature),
                       str(int(self.config.qm_params.use_solvation)),
                       self.config.qm_params.solvent_name,
                       density_file)
        #write QM output exactly as in the serial code
        out_path = (self.config.output_folder / self.config.output_folder_qm / f"{residue.residue_id}_ip.out")
        QM.Write2File(energy,rmsd,totalpol,atmlist,charges,polarisabilities,converged,out_path)
        if not converged:
            self.logger.warning(f"Convergence issue for {residue.residue_id}")
        #store IP file location
        residue.ip_file = f"{res_ip}.xyz"
        return energy

    def _calculate_pair_energies(self, ligand_states: List[LigandState], residues: List[RES.ResidueData]) -> List[List[Optional[float]]]:
        """Calculate energies for ligand-residue pairs."""
        self.logger.info("Calculating pair energies")
        if self.config.edda:
            self.logger.info("EDDA analysis started")
        pair_energies = []
        for lig_idx,ligand in enumerate(ligand_states):
            ligand_file = self._get_ligand_file(ligand)
            if self.config.parallel_processing:
                lig_energies = self._calculate_ligand_pairs_parallel(ligand,ligand_file,residues)
            else:
                lig_energies = self._calculate_ligand_pairs_sequential(ligand,ligand_file,residues)
            pair_energies.append(lig_energies)
        return pair_energies

    def _get_ligand_file(self, ligand: LigandState) -> str:
        """Get appropriate ligand file path."""
        ligand_file = str(ligand.file_path)
        if self.config.use_ip_structure:
            #look for IP structure first
            if ligand_file.endswith(".xyz"):
                ip_name = Path(ligand_file).name.replace(".xyz","_ip.xyz")
            elif ligand_file.endswith(".sdf"):
                ip_name = Path(ligand_file).name.replace(".sdf","_ip.xyz")
            else:
                ip_name = Path(ligand_file).name + "_ip.xyz"
            ip_folder_path = self.file_manager.get_path("ip",ip_name)
            if ip_folder_path.exists():
                return str(ip_folder_path)
        #FIXED: Check if the original file has been moved to start folders
        original_name = Path(ligand_file).name
        #check start_sdf folder first (since files are moved there)
        if original_name.endswith(".sdf"):
            start_sdf_path = self.file_manager.get_path("start_sdf",original_name)
            if start_sdf_path.exists():
                return str(start_sdf_path)
        #check start_xyz folder
        if original_name.endswith(".xyz"):
            start_xyz_path = self.file_manager.get_path("start_xyz",original_name)
            if start_xyz_path.exists():
                return str(start_xyz_path)
        #check start_pdb folder
        pdb_name = original_name.replace(".sdf",".pdb").replace(".xyz",".pdb")
        start_pdb_path = self.file_manager.get_path("start_pdb",pdb_name)
        if start_pdb_path.exists():
            return str(start_pdb_path)
        #fallback to original path if it still exists
        if Path(ligand_file).exists():
            return ligand_file
        #last resort: look in geometries folder
        geom_path = self.file_manager.get_path("geometries",original_name)
        if geom_path.exists():
            return str(geom_path)
        raise FileNotFoundError(f"Cannot find ligand file: {ligand_file}")

    def _calculate_ligand_pairs_sequential(self, ligand: LigandState, ligand_file: str, residues: List[RES.ResidueData]) -> List[Optional[float]]:
        """Calculate pair energies for one ligand sequentially."""
        energies = []
        for residue in residues:
            energy,clash = self._calculate_single_pair_energy(ligand,ligand_file,residue)
            energies.append(energy)
            self.clash_analysis_list[ligand.ligand_index].append(clash)
        return energies

    def _calculate_ligand_pairs_parallel(self, ligand: LigandState, ligand_file: str, residues: List[RES.ResidueData]) -> List[Optional[float]]:
        """Calculate pair energies for one ligand in parallel."""
        energies = [None]*len(residues)
        with ProcessPoolExecutor(max_workers = self.config.max_workers) as executor:
            #create futures in order
            futures = [executor.submit(self._calculate_single_pair_energy,ligand,ligand_file,res)
                       for res in residues]
            #collect results in order
            for idx,future in enumerate(futures):
                try:
                    energies[idx],clash = future.result()
                    self.clash_analysis_list[ligand.ligand_index].append(clash)
                except Exception as e:
                    self.logger.error(f"Failed to calculate pair energy: {e}")
                    energies[idx] = None
                    self.clash_analysis_list[ligand.ligand_index].append(None)
        return energies

    def _calculate_single_pair_energy(self, ligand: LigandState, ligand_file: str, residue: RES.ResidueData) -> Tuple[float, ClashAnalysis]:
        """Calculate energy for a single ligand-residue pair."""
        try:
            #read geometries with a simple per-file cache to avoid repeated I/O
            #and parsing when the same ligand_file is reused across residues.
            if not hasattr(self,"_ligand_geometry_cache"): self._ligand_geometry_cache = {}
            cache = self._ligand_geometry_cache
            if ligand_file in cache:
                cached_atoms,cached_geom = cache[ligand_file]
                #recreate fresh mutable lists so no call can accidentally
                #modify the cached data.
                ligand_atoms = list(cached_atoms)
                ligand_geom = [list(coords) for coords in cached_geom]
            else:
                ligand_atoms,ligand_geom = GEOM.ReadGeometry(ligand_file)
                #store an immutable snapshot in the cache
                cache[ligand_file] = (tuple(ligand_atoms),
                                      tuple(tuple(coords) for coords in ligand_geom))
            #read residue ip structure
            ip_res_atoms = residue.processed_atoms
            ip_res_geom = residue.processed_geometry
            if self.config.use_ip_residues:
                ip_res_atoms,ip_res_geom = GEOM.ReadGeometry(residue.ip_file)
                residue.processed_atoms = ip_res_atoms.copy()
                residue.processed_geometry = ip_res_geom.copy()
            if not ip_res_atoms:
                return None,None
            #create pair file
            ligand_name = Path(ligand_file).stem
            pair_file = self.file_manager.get_path("geometries",f"{ligand_name}_{residue.residue_id}")
            #combine and write structures
            pair_atoms = ligand_atoms + ip_res_atoms
            pair_geom = ligand_geom + ip_res_geom
            GEOM.Write2XYZ(pair_atoms,pair_geom,str(pair_file))
            #use IN_POCKET for pair calculations
            total_charge = ligand.charge + residue.charge
            calc_type = (
                CalculationType.IN_POCKET
                if self.config.do_in_pocket
                else CalculationType.FULL_OPT
            )
            if self.config.evaluate_only: 
                calc_type = CalculationType.BASE_ENERGY
            (energy,_,_,_,_,_,converged) = self.qm_calc.calculate_energy(str(pair_file) + ".xyz",
                                                                         total_charge,
                                                                         calc_type,  #IN_POCKET when do_in_pocket is True
                                                                         str(pair_file) + "_ip")
            if not converged:
                #use file name or SMILES for identification
                lig_label = (ligand.smiles
                             if ligand.smiles
                             else Path(ligand.file_path).name)
                self.logger.warning("Convergence issue for pair %s-%s",lig_label,residue.residue_id)
            #recentering ip geometry
            GEOM.ReCenterAndWrite(str(pair_file) + ".xyz",str(pair_file) + "_ip" + ".xyz",str(pair_file) + "_ip")
            if self.config.qm_params.calculate_density:
                density_file = self._density_prefix(str(pair_file) + "_ip")
                QM.GetDensity(str(pair_file) + "_ip.xyz",
                              str(total_charge),
                              str(self.config.qm_params.electronic_temperature),
                              str(int(self.config.qm_params.use_solvation)),
                              self.config.qm_params.solvent_name,
                              density_file)
            clash = ClashAnalysis()
            clash.residue_id = residue.residue_id
            edda_ok = False
            if self.config.edda:
                edda_ok = True
                try:
                    if getattr(self.config, "debug", False): self.logger.info(f"EDDA processing for {str(pair_file)}_ip")
                    edda_out,edda_err = QM.RunEDDA(str(pair_file) + "_ip" + ".xyz",
                                                   str(total_charge),
                                                   ligand_file,
                                                   str(ligand.charge),
                                                   residue.ip_file,
                                                   str(residue.charge))
                    RepulsionE,InteractionE,LocalRepulsionA,LocalRepulsionB,LocalInteractionA,LocalInteractionB = QM.ParseEDDAOutput(edda_out,ligand_atoms,ip_res_atoms)
                    clash.repulsion_energy = RepulsionE
                    clash.EDDA_binding_energy = InteractionE
                    heavyatm_list_A, heavyatm_rep_A, heavyatm_int_A, shortcontact_A = SCORE.ShortContactScore(ligand_atoms,LocalRepulsionA,LocalInteractionA,ligand_geom)
                    heavyatm_list_B, heavyatm_rep_B, heavyatm_int_B, shortcontact_B = SCORE.ShortContactScore(ip_res_atoms,LocalRepulsionB,LocalInteractionB,ip_res_geom)
                    clash.heavy_atom_list = [*heavyatm_list_A,*heavyatm_list_B]
                    clash.repulsion_contacts = [*heavyatm_rep_A,*heavyatm_rep_B]
                    clash.interaction_contacts = [*heavyatm_int_A,*heavyatm_int_B]
                    clash.short_contact_indices = [*shortcontact_A,*shortcontact_B]
                    clash.split_index = len(shortcontact_A)
                    clash.residue_id = residue.residue_id
                except Exception as _edda_e:
                    import traceback
                    #mark EDDA failure on the clash object (so parent can detect it)
                    clash.edda_failed = True
                    clash.edda_error = str(_edda_e)
                    edda_ok = False
                    self.logger.error(f"EDDA failed for lig = {ligand.ligand_index} residue = {residue.residue_id}: {_edda_e}")
                    #dump debug bundle ONLY on EDDA error
                    dbg_dir = Path(self.file_manager.output_folder) / "debug"
                    dbg_dir.mkdir(parents = True,exist_ok = True)
                    tag = f"edda_fail__lig{ligand.ligand_index}__{residue.residue_id}"
                    #1) traceback + key sizes
                    with open(dbg_dir / f"{tag}.log", "w") as fh:
                        fh.write(f"PAIR_FILE_BASE: {str(pair_file)}\n")
                        fh.write(f"LIGAND_FILE: {str(ligand_file)}\n")
                        fh.write(f"RES_IP_FILE: {str(residue.ip_file)}\n")
                        fh.write(f"TOTAL_CHARGE: {total_charge} (lig = {ligand.charge}, res = {residue.charge})\n")
                        fh.write(f"EXCEPTION: {_edda_e}\n\n")
                        fh.write(traceback.format_exc())
                        fh.write("\n\n--- SHAPES ---\n")
                        fh.write(f"len(ligand_atoms) = {len(ligand_atoms)}\n")
                        fh.write(f"len(ligand_geom) = {len(ligand_geom)}\n")
                        fh.write(f"len(residue.processed_atoms) = {len(ip_res_atoms) if ip_res_atoms else 0}\n")
                        fh.write(f"len(residue.processed_geometry) = {len(ip_res_geom) if ip_res_geom else 0}\n")
                        fh.write(f"heavy(ligand_atoms) = {sum(1 for el in ligand_atoms if el != 'H')}\n")
                        fh.write(f"heavy(residue.processed_atoms) = {sum(1 for el in (ip_res_atoms or []) if el != 'H')}\n")
                    #2) raw EDDA stdout/stderr (if available)
                    try:
                        with open(dbg_dir / f"{tag}.edda.out", "w") as fh:
                            fh.write(edda_out if isinstance(edda_out,str) else str(edda_out))
                    except Exception: pass
                    try:
                        with open(dbg_dir / f"{tag}.edda.err","w") as fh:
                            fh.write(edda_err if isinstance(edda_err,str) else str(edda_err))
                    except Exception: pass
                #atom labels for downstream reporting
                #only attempt label alignment if EDDA succeeded and the score arrays exist.
                if not edda_ok:
                    return energy,clash
                lig_labels: List[str] = []
                try:
                    if str(ligand_file).lower().endswith(".pdb") and os.path.exists(str(ligand_file)):
                        with open(str(ligand_file),"r") as _fh:
                            for _ln in _fh:
                                if _ln.startswith(("ATOM","HETATM")):
                                    _name = _ln[12:16].strip()
                                    if _name:
                                        lig_labels.append(_name)
                    #if parsing failed or file type doesn't carry names, fall back to element+index
                    if len(lig_labels) != len(ligand_atoms):
                        lig_labels = [f"{el}{i+1}" for i,el in enumerate(ligand_atoms)]
                except Exception:
                    lig_labels = [f"{el}{i+1}" for i,el in enumerate(ligand_atoms)]
                res_labels: List[str] = []
                try:
                    #residue.atom_types carries PDB atom names for the original residue atoms
                    base = list(getattr(residue,"atom_types",[]) or [])
                    n_proc = len(ip_res_atoms) if ip_res_atoms else 0
                    if base and n_proc:
                        if len(base) < n_proc:
                            base = base + [f"CAP{i + 1}" for i in range(n_proc - len(base))]
                        res_labels = base[:n_proc]
                    else:
                        res_labels = [f"{el}{i + 1}" for i,el in enumerate(ip_res_atoms)]
                except Exception:
                    res_labels = [f"{el}{i + 1}" for i,el in enumerate(ip_res_atoms)]
                #prefix residue labels with residue id so output is unambiguous
                res_labels = [f"{residue.residue_id}:{nm}" for nm in res_labels]
                #build labels aligned to EDDA local scores (assume scores are heavy-atom-based)
                lig_heavy = [(i,el) for i,el in enumerate(ligand_atoms) if el != "H"]
                res_heavy = [(i,el) for i,el in enumerate(ip_res_atoms) if el != "H"]
                #if EDDA scores are heavy-atom based, these lengths should match
                if len(lig_heavy) != len(shortcontact_A):
                    #fallback: trust EDDA length, label as element + rank
                    lig_labels = [f"LIG:{el}{k + 1}" for k, el in enumerate(heavyatm_list_A)]
                else:
                    lig_labels = [f"LIG:{el}{i + 1}" for (i, el) in lig_heavy]
                #residue atomtype names from residue.atom_types if possible
                res_atomtypes = list(getattr(residue,"atom_types",[]) or [])
                res_proc_n = len(ip_res_atoms) if ip_res_atoms else 0
                if len(res_atomtypes) < res_proc_n:
                    res_atomtypes = res_atomtypes + [f"CAP{i + 1}" for i in range(res_proc_n - len(res_atomtypes))]
                if len(res_heavy) != len(shortcontact_B):
                    res_labels = [f"{residue.residue_id}:X{k + 1}" for k in range(len(shortcontact_B))]
                else:
                    res_labels = []
                    for (i,el) in res_heavy:
                        atn = res_atomtypes[i] if i < len(res_atomtypes) else f"X{i + 1}"
                        res_labels.append(f"{residue.residue_id}:{atn}")
                clash.atom_labels = [*lig_labels,*res_labels]
            return energy,clash
        except Exception as e:
            self.logger.error(f"Pair calculation failed: {e}")
            return None,None

    def _calculate_cap_corrections(self, ligand_states: List[LigandState], residues: List[RES.ResidueData]) -> Dict[str, float]:
        """Calculate cap corrections for unique base residues."""
        self.logger.info("Calculating cap corrections")
        _capcorr_skipped = 0
        _capcorr_done = 0
        #get formaldehyde reference
        formaldehyde_energy, formaldehyde_file = FORM.Eformaldehyde(self.config.qm_params.solvent_name
                                                                    if self.config.qm_params.use_solvation
                                                                    else "gas")
        #robust path resolution (avoid CWD-dependent relative paths) ---
        #Eformaldehyde() may return a relative path like "util/formaldehyde/CH2O_water"
        #but during the run the CWD is typically the system output folder.
        fpath = Path(formaldehyde_file)
        #if relative, resolve relative to the *ipa-l repo* (directory containing in_pocket.py)
        if not fpath.is_absolute():
            base_dir = Path(__file__).resolve().parent
            fpath = (base_dir / fpath).resolve()
        #util/Geometry.ReadXYZ() appends ".xyz" internally, so we check "<path>.xyz"
        if not Path(str(fpath) + ".xyz").is_file():
            #fallback to canonical location in case Eformaldehyde returned something odd
            base_dir = Path(__file__).resolve().parent
            fallback = (base_dir / "util" / "formaldehyde" / "CH2O_water").resolve()
            if Path(str(fallback) + ".xyz").is_file():
                fpath = fallback
            else:
                raise FileNotFoundError("Formaldehyde reference geometry not found.\n"
                                        f"Expected: {str(fpath)}.xyz\n"
                                        f"Fallback: {str(fallback)}.xyz\n"
                                        "Fix: ensure util/formaldehyde/CH2O_water.xyz exists in the IPA repo "
                                        "(same repo where in_pocket.py lives).")
        #keep the original convention: path without ".xyz"
        formaldehyde_file = str(fpath)
        #group residues by base identifier
        base_residues = {}
        for residue in residues:
            base_id = self._get_base_residue_id(residue)
            if base_id not in base_residues:
                base_residues[base_id] = residue
        #calculate corrections
        corrections = {}
        for base_id,residue in base_residues.items():
            cap_atoms_n = len(residue.cap_atoms) if residue.cap_atoms else 0
            cap_geom_n = len(residue.cap_geometry) if residue.cap_geometry else 0
            std_aa = residue.is_standard_amino_acid()
            std_na = NA.IsBaseNucleicAcid(residue.name) or NA.IsExtendedNucleicAcid(residue.name)
            self.logger.debug("CapCorr candidate: base_id = %s residue = %s standardAA = %s cap_atoms = %d cap_geom = %d protonated_file = %s",
                              base_id,
                              getattr(residue,"residue_id","<?>"),
                              std_aa,
                              cap_atoms_n,
                              cap_geom_n,
                              str(getattr(residue,"protonated_file","")))
            if ((not std_aa) and (not std_na)) or (cap_atoms_n == 0):
                _capcorr_skipped += 1
                self.logger.debug("CapCorr skip: base_id = %s residue = %s (standardAA = %s cap_atoms = %d)",
                                  base_id,
                                  getattr(residue,"residue_id","<?>"),
                                  std_aa,
                                  cap_atoms_n)
                corrections[base_id] = {i: 0.0 for i in range(len(ligand_states))}
                continue
            try:
                if std_aa: lig_correction = self._calculate_amino_acid_single_cap_correction(residue,ligand_states,formaldehyde_energy,formaldehyde_file)
                elif std_na: lig_correction = self._calculate_nucleic_acid_single_cap_correction(residue,ligand_states)
            except Exception as e:
                #HARD FAIL (for this in_pocket run) with maximal context for debugging
                self.logger.error("CapCorr FAILED: base_id = %s residue = %s cap_atoms = %d cap_geom = %d protonated_file = %s error = %s",
                                  base_id,
                                  getattr(residue,"residue_id","<?>"),
                                  cap_atoms_n,
                                  cap_geom_n,
                                  str(getattr(residue,"protonated_file","")),
                                  repr(e),
                                  exc_info = True)
                raise
            corrections[base_id] = lig_correction
            _capcorr_done += 1
        self.logger.info("CapCorr summary: done = %d skipped = %d total_unique_base = %d",
                         _capcorr_done,
                         _capcorr_skipped,
                         (_capcorr_done + _capcorr_skipped))
        return corrections

    def _get_base_residue_id(self, residue: RES.ResidueData) -> str:
        """Get base residue identifier without protonation state."""
        #extract base ID properly - should be "TYR_B253" not "TYR_B253_1"
        parts = residue.residue_id.split("_")
        if len(parts) >= 2:
            return f"{parts[0]}_{parts[1]}"  #e.g., "TYR_B253"
        return residue.residue_id
    
    def _generate_methanol_cap(self, residue: RES.ResidueData, rOH = 0.967, rCH = 1.3):
        """Function generating a methanol cap for NA structures"""
        atoms_ip,geometry_ip = GEOM.ReadXYZ(residue.ip_file)
        cap_protons = []
        phosphorous = []
        cap_C_coord = []
        cap_O_coord = []
        rCHrCH = rCH*rCH
        #find the methanol carbon and oxygen
        if residue.cap_atoms[0] == "C": 
            cap_C_coord = residue.cap_geometry[0]
            cap_O_coord = residue.cap_geometry[1]
        elif residue.cap_atoms[1] == "C": 
            cap_C_coord = residue.cap_geometry[1]
            cap_O_coord = residue.cap_geometry[0]
        #find the methanol cap protons and the phosphorous
        for idatm in range(len(atoms_ip)):
            if atoms_ip[idatm] == "P": phosphorous = geometry_ip[idatm]
            elif atoms_ip[idatm] != "H": continue
            d2 = 0.0
            for idcoord in range(3):
                aux = cap_C_coord[idcoord] - geometry_ip[idatm][idcoord]
                d2 += aux*aux
            if d2 < rCHrCH:
                cap_protons.append(geometry_ip[idatm])
        #if we found the correct number of atoms, then proceed
        if len(cap_protons) == 3:
            cap_geom = copy.deepcopy(residue.cap_geometry)
            cap_atoms = copy.deepcopy(residue.cap_atoms)
            #add protons
            for coord in cap_protons:
                cap_geom.append(coord)
                cap_atoms.append("H")
            #now the "phosphours"
            distance = 0.0
            rAB = [0.0,0.0,0.0]
            for idcoord in range(3):
                rAB[idcoord] = cap_O_coord[idcoord] - phosphorous[idcoord]
                distance += rAB[idcoord]*rAB[idcoord]
            renorm = rOH/math.sqrt(distance)
            for idcoord in range(3):
                rAB[idcoord] *= -renorm
                rAB[idcoord] += cap_O_coord[idcoord]
            cap_geom.append(rAB)
            cap_atoms.append("H")
            return cap_atoms,cap_geom
        else: return None,None

    def _calculate_cap_corrections_general(self, ligand_states: List[LigandState], cap_atoms, cap_geom, cap_base_id, calc_type, cap_ref_energy):
        """General function to calculate cap corrections for amino and nucleic acids"""
        corrections: Dict[int,float] = {}
        for i,ligand in enumerate(ligand_states):
            ligand_file = self._get_ligand_file(ligand)
            if not ligand_file or not Path(ligand_file).exists():
                continue
            #cached geometry block to replace ligand_atoms,ligand_geom = ReadGeometry(ligand_file)) ---
            #reuse ligand geometry cache to avoid repeatedly reading the same
            #ligand file from disk across residues / cap corrections.
            if not hasattr(self,"_ligand_geometry_cache"):
                self._ligand_geometry_cache = {}
            cache = self._ligand_geometry_cache
            if ligand_file in cache:
                cached_atoms,cached_geom = cache[ligand_file]
                #work on fresh mutable copies so callers cannot modify the cache
                ligand_atoms = list(cached_atoms)
                ligand_geom = [list(coords) for coords in cached_geom]
            else:
                ligand_atoms,ligand_geom = GEOM.ReadGeometry(ligand_file)
                #store an immutable snapshot in the cache
                cache[ligand_file] = (tuple(ligand_atoms),tuple(tuple(coords) for coords in ligand_geom))
            #create cap–ligand pair:
            #<ligand_name>_<base_id>_capH.xyz
            #e.g. "B9000_ip_THR_B315_capH.xyz" (matching the serial naming)
            ligand_name = Path(ligand_file).stem
            pair_file = self.file_manager.get_path("geometries",f"{ligand_name}_{cap_base_id}_capH")
            pair_atoms = ligand_atoms + cap_atoms
            pair_geom = ligand_geom + cap_geom
            GEOM.Write2XYZ(pair_atoms,pair_geom,str(pair_file))
            #calculate energy
            (energy_cappair,rmsd_cappair,totalpol_cappair,atmlist_cappair,charges_cappair,polarisabilities_cappair,converged_cappair) = self.qm_calc.calculate_energy(str(pair_file) + ".xyz",ligand.charge,calc_type,str(pair_file) + "_ip")
            if calc_type is CalculationType.IN_POCKET:
                #recentering ip geometry
                GEOM.ReCenterAndWrite(str(pair_file) + ".xyz",str(pair_file) + "_ip.xyz",str(pair_file) + "_ip")
                if self.config.qm_params.calculate_density:
                    density_file = self._density_prefix(str(pair_file) + "_ip")
                    QM.GetDensity(str(pair_file) + "_ip.xyz",
                                  str(ligand.charge),
                                  str(self.config.qm_params.electronic_temperature),
                                  str(int(self.config.qm_params.use_solvation)),
                                  self.config.qm_params.solvent_name,
                                  density_file)
            if energy_cappair is not None and ligand.energies.get("ip") is not None:
                #write QM output for this cap–ligand pair
                out_path = (self.config.output_folder / self.config.output_folder_qm / f"{cap_base_id}_{ligand_name}_capH_ip.out")
                QM.Write2File(energy_cappair,rmsd_cappair,totalpol_cappair,atmlist_cappair,charges_cappair,polarisabilities_cappair,converged_cappair,out_path)
                #cap correction: Ecapligand - Eformaldehyde - Eligand
                correction = (energy_cappair - cap_ref_energy - ligand.energies["ip"])
                corrections[i] = correction
        #return corrections mapping (or 0.0 if nothing computed)
        return corrections if corrections else 0.0
    
    def _calculate_nucleic_acid_single_cap_correction(self, residue: RES.ResidueData, ligand_states: List[LigandState], rOH = 0.967, rCH = 1.3):
        """Calculate cap correction for a single nucleic acid residue and write *_capH_ip.out."""
        #get the protonated cap model
        cap_atoms,cap_geom = self._generate_methanol_cap(residue,rOH,rCH)
        if cap_atoms is not None:
            #use base residue ID (e.g. TYR_B253) for BOTH geometry and output naming,
            base_id = self._get_base_residue_id(residue)
            residue.cap_atoms = cap_atoms
            residue.cap_geometry = cap_geom
            #create cap structure: <base_id>_capH.xyz
            cap_file = self.file_manager.get_path("geometries",f"{base_id}_capH")
            self.logger.debug("CapCorr output: residue = %s base_id = %s cap_atoms = %d cap_geom = %d last_atom = %s",
                              getattr(residue,"residue_id","<?>"),
                              base_id,
                              len(residue.cap_atoms) if residue.cap_atoms else 0,
                              len(residue.cap_geometry) if residue.cap_geometry else 0,
                              residue.cap_atoms[-1] if residue.cap_atoms else "<?>")
            #write the geometry
            GEOM.Write2XYZ(residue.cap_atoms,residue.cap_geometry,str(cap_file))
            #contrary to the amino-acid case, we have here already the oriented cap, se we can proceed with evaluation
            #lets get the base energy of this cap
            calc_type = (CalculationType.BASE_ENERGY
                         if self.config.use_ip_structure
                         else CalculationType.IN_POCKET)
            (methanol_energy,_,_,_,_,_,converged_cap) = self.qm_calc.calculate_energy(str(cap_file) + ".xyz","0",calc_type,str(cap_file) + "_ip")
            if converged_cap: 
                #calculate corrections over all ligand states
                return self._calculate_cap_corrections_general(ligand_states,cap_atoms,cap_geom,base_id,calc_type,methanol_energy)
            else:
                self.logger.error("No calculation possible on the methanol cap")
                raise ValueError("No calculation possible on the methanol cap")
        else: 
            self.logger.error("No cap atoms detected")
            raise ValueError("No cap atoms detected")

    def _calculate_amino_acid_single_cap_correction(self, residue: RES.ResidueData, ligand_states: List[LigandState], formaldehyde_energy: float, formaldehyde_file: str) -> float:
        """Calculate cap correction for a single amino acid residue and write *_capH_ip.out."""
        #use base residue ID (e.g. TYR_B253) for BOTH geometry and output naming,
        base_id = self._get_base_residue_id(residue)
        #create cap structure: <base_id>_capH.xyz
        cap_file = self.file_manager.get_path("geometries",f"{base_id}_capH")
        #add proton to cap
        self.logger.debug("CapCorr AddCapProton input: residue = %s base_id = %s cap_atoms = %d cap_geom = %d protonated_file = %s",
                          getattr(residue,"residue_id","<?>"),
                          base_id,
                          len(residue.cap_atoms) if residue.cap_atoms else 0,
                          len(residue.cap_geometry) if residue.cap_geometry else 0,
                          str(getattr(residue,"protonated_file","")))
        GEOM.AddCapProton(residue.cap_atoms,residue.cap_geometry,str(residue.protonated_file))
        self.logger.debug("CapCorr AddCapProton output: residue = %s base_id = %s cap_atoms = %d cap_geom = %d last_atom = %s",
                          getattr(residue,"residue_id","<?>"),
                          base_id,
                          len(residue.cap_atoms) if residue.cap_atoms else 0,
                          len(residue.cap_geometry) if residue.cap_geometry else 0,
                          residue.cap_atoms[-1] if residue.cap_atoms else "<?>")
        GEOM.Write2XYZ(residue.cap_atoms,residue.cap_geometry,str(cap_file))
        try:
            #cap_file is a prefix; the actual file is "<cap_file>.xyz"
            cap_xyz_path = str(cap_file) + ".xyz"
            cap_atoms_chk,cap_geom_chk = GEOM.ReadGeometry(cap_xyz_path)
            self.logger.debug("CapCorr cap.xyz check: residue = %s base_id = %s cap_xyz = %s atoms = %d geom = %d first_atoms = %s",
                              getattr(residue,"residue_id","<?>"),
                              base_id,
                              cap_xyz_path,
                              len(cap_atoms_chk) if cap_atoms_chk else 0,
                              len(cap_geom_chk) if cap_geom_chk else 0,
                              ",".join(cap_atoms_chk[:5]) if cap_atoms_chk else "<?>")
            #Formaldehyde.ReorientCap uses geometry[0], [1], [2] -> must have >= 3 atoms
            if (not cap_atoms_chk) or (len(cap_atoms_chk) < 3) or (not cap_geom_chk) or (len(cap_geom_chk) < 3):
                self.logger.error("CapCorr INVALID cap.xyz (needs >=3 atoms): residue = %s base_id = %s cap_xyz = %s atoms = %d geom = %d",
                                  getattr(residue,"residue_id","<?>"),
                                  base_id,
                                  cap_xyz_path,
                                  0 if not cap_atoms_chk else len(cap_atoms_chk),
                                  0 if not cap_geom_chk else len(cap_geom_chk))
                raise ValueError(f"Invalid cap xyz for ReorientCap (needs >=3 atoms): {cap_xyz_path}")
        except Exception as e:
            self.logger.error("CapCorr pre-ReorientCap validation failed: residue = %s base_id = %s cap_file = %s error = %s",
                              getattr(residue,"residue_id","<?>"),
                              base_id,
                              str(cap_file),
                              repr(e),
                              exc_info = True)
            raise
        FORM.ReorientCap(formaldehyde_file,str(cap_file))
        #calculate corrections over all ligand states
        cap_atoms,cap_geom = GEOM.ReadGeometry(str(cap_file) + ".xyz")
        calc_type = (CalculationType.BASE_ENERGY
                     if self.config.use_ip_structure
                     else CalculationType.IN_POCKET)
        return self._calculate_cap_corrections_general(ligand_states,cap_atoms,cap_geom,base_id,calc_type,formaldehyde_energy)

    def _calculate_single_binding(self,ligand: LigandState,residue: RES.ResidueData,residue_energy: Optional[float],pair_energy: Optional[float],cap_correction: float) -> Optional[BindingResult]:
        """Calculate binding energy for a single ligand-residue pair."""
        ligand_energy = ligand.energies.get("ip")
        if None in [ligand_energy,residue_energy,pair_energy]:
            return None
        binding_energy_hartree = (pair_energy - residue_energy - ligand_energy - cap_correction)
        binding_energy_kcal = binding_energy_hartree*self.Eh2kcalmol
        cap_correction_kcal = cap_correction*self.Eh2kcalmol
        return BindingResult(residue_name = residue.residue_id,
                             residue_state = residue.protonation_state,
                             binding_energy = binding_energy_kcal,
                             cap_correction = cap_correction_kcal,
                             converged = True,
                             components = {"pair": pair_energy,
                                           "residue": residue_energy,
                                           "ligand": ligand_energy,
                                           "cap": cap_correction}
                             #clash is now @ calculate_binding_energies() when edda == True
        )

class InPocketAnalyzer:
    """Main class for performing in-pocket analysis."""

    def __init__(self, config: AnalysisConfig):
        self.config = config
        self.file_manager = FileManager(config.output_folder,config.geometry_folder,config.output_folder_qm)
        self.logger = self._setup_logging()
        self.qm_calculator = QMCalculator(config.qm_params,self.file_manager)
        self.residue_processor = ResidueProcessor(config,self.file_manager)
        self.binding_calculator = BindingCalculator(config,self.qm_calculator,self.file_manager)
        self.qm_calc = self.qm_calculator
        self.Eh2kcalmol = UC.Hartree2Kcalmol

    def _density_prefix(self, base_tag: str) -> str:
        """Return the output prefix (no extension) for density files.
           We ALWAYS write .cub under <output>/cub_files to avoid polluting geometries/.
           base_tag is typically something like:
             - '<...>/geometries/B9000_ip' or
             - '<...>/geometries/TYR_B253_ip' or
             - '<...>/geometries/<pair>_ip'"""
        cub_dir = self.config.output_folder / "cub_files"
        cub_dir.mkdir(parents = True,exist_ok = True)
        #keep only the filename part from base_tag (strip directories)
        name = Path(base_tag).name
        return str(cub_dir / f"{name}_dens")

    def _log_timing(self, step_name: str, start_time: float) -> None:
        """Log the timing for a step."""
        elapsed = time.time() - start_time
        self.logger.info(f"{step_name} completed in {elapsed:.2f} seconds")

    def _log_interaction_summary(self,ligand_states: List[LigandState],binding_results: List[List[BindingResult]]) -> None:
        """Log interaction energy summary for each ligand state."""
        self.logger.info("\n" + "="*60)
        self.logger.info("INTERACTION ENERGY SUMMARY")
        self.logger.info("="*60)
        for lig_idx,(ligand,results) in enumerate(zip(ligand_states,binding_results)):
            #calculate total interaction energy (excluding flipped residues)
            total_interaction = 0.0
            non_flipped_results = []
            for res_j,result in enumerate(results):
                if "_flip" not in result.residue_name:
                    total_interaction += result.binding_energy
                    non_flipped_results.append(result)
            #sort by binding energy (most favorable first)
            sorted_results = sorted(non_flipped_results,key = lambda r: r.binding_energy)
            #log ligand state info
            self.logger.info(f"Ligand State {lig_idx + 1}:")
            if ligand.smiles:
                self.logger.info(f"  SMILES: {ligand.smiles}")
            self.logger.info(f"  Charge: {ligand.charge}")
            self.logger.info(f"  Total Interaction Energy: {total_interaction:.3f} kcal/mol")
            nshow = min(3,len(sorted_results))
            if nshow > 0:
                self.logger.info(f"  Top {nshow} Most Favorable Interactions:")
                for i,result in enumerate(sorted_results[:nshow],1):
                    self.logger.info(f"    {i}. {result.residue_name}: {result.binding_energy:.3f} kcal/mol")
            else:
                self.logger.info("  No interactions found")            

    def _setup_logging(self) -> logging.Logger:
        """Setup logging configuration."""
        logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        if not logger.handlers:  #avoid duplicate handlers
            #create formatter
            formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
            #console handler
            console_handler = logging.StreamHandler()
            console_handler.setFormatter(formatter)
            logger.addHandler(console_handler)
            #file handler
            log_file = self.config.output_folder / "in_pocket_analysis.log"
            file_handler = logging.FileHandler(log_file)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
            level = logging.DEBUG if getattr(self.config,"debug",False) else logging.INFO
            logger.setLevel(level)
            logger.propagate = False
            #ensure handlers follow the same level
            console_handler.setLevel(level)
            file_handler.setLevel(level)
            #propagate handlers to root so other class loggers (BindingCalculator, etc.)
            root = logging.getLogger()
            root.setLevel(level)
            for h in logger.handlers:
                if h not in root.handlers:
                    root.addHandler(h)
        return logger

    def validate_inputs(self) -> None:
        """Validate input parameters and files."""
        self.logger.info("Validating inputs...")
        #check PDB file
        if not self.config.pdb_file or not self.config.pdb_file.exists():
            #try with extensions
            if ((self.config.pdb_file.suffix != ".pdb") and (self.config.pdb_file.suffix != ".cif") and (self.config.pdb_file.suffix != ".mmcif")):
                pdb_with_ext = self.config.pdb_file.with_suffix(".pdb")
                pdb_with_ext_cif = self.config.pdb_file.with_suffix(".cif")
                pdb_with_ext_mmcif = self.config.pdb_file.with_suffix(".mmcif")
                if pdb_with_ext.exists():
                    self.config.pdb_file = pdb_with_ext
                elif pdb_with_ext_cif.exists():
                    self.config.pdb_file = pdb_with_ext_cif
                elif pdb_with_ext_mmcif.exists():
                    self.config.pdb_file = pdb_with_ext_mmcif
                else:
                    raise InPocketAnalyzerError(f"Structure file not found: {self.config.pdb_file}")
        #determine file format
        if self.config.pdb_file.suffix.lower() in [".cif",".mmcif"]:
            self.config.mmcif_file = True
            self.logger.info("Detected mmCIF format")
        #validate parameters
        if self.config.critical_distance <= 0:
            raise InPocketAnalyzerError("Critical distance must be positive")
        if self.config.qm_params.alpha <= 0:
            raise InPocketAnalyzerError("Alpha parameter must be positive")
        if not 0 <= self.config.ph <= 14:
            self.logger.warning(f"pH value {self.config.ph} is outside typical range [0-14]")
        if self.config.active_residue in ("NAN",None):
            if getattr(self.config,"active_residue_name",None) in ("0",None,""):
                raise InPocketAnalyzerError("Either --active-residue or --active-residue-name must be specified")
        self.logger.info("Input validation completed successfully")

    def read_structure(self) -> Tuple[RES.ResidueData, List[RES.ResidueData]]:
        """Read and parse the input structure file."""
        self.logger.info(f"Reading structure from {self.config.pdb_file}")
        try:
            if self.config.mmcif_file:
                geometry,atoms,atmtp,residues,resnumbr,chainid,restp,seqidcodes,pdblines = (GEOM.ReadCIFStructure(str(self.config.pdb_file)))
            else:
                geometry,atoms,atmtp,residues,resnumbr,chainid,restp,seqidcodes,pdblines = (GEOM.ReadPDB_base(str(self.config.pdb_file), self.config.altloc))
        except Exception as e:
            raise InPocketAnalyzerError(f"Failed to read structure file: {e}")
        #process altloc if needed
        if self.config.check_altloc and self.config.mmcif_file:
            try:
                geometry,atoms,atmtp,residues,resnumbr,chainid,restp,seqidcodes,pdblines = (
                    GEOM.SelectOneAltLocOnly(geometry,atoms,atmtp,residues,resnumbr,chainid,restp,seqidcodes,pdblines,self.config.altloc))
                altloc_warnings = GEOM.GetAltLocSelectionWarnings(pdblines)
                for msg in altloc_warnings:
                    self.logger.warning(msg)
            except Exception as e:
                self.logger.warning(f"Failed to process altloc: {e}")
        #strip hydrogens already present in the input structure if not requested otherwise
        #(RDKit protonation tool chokes on "pre-protonated" PDBs)
        if not getattr(self.config,"keep_input_h",False):
            try:
                keep_idx = []
                for i,atp in enumerate(atmtp):
                    #conservative: only drop true hydrogen.
                    #(some pipelines may encode atmtp in different ways; we also guard on atom name.)
                    atp_u = str(atp).strip().upper()
                    atom_u = str(atoms[i]).strip().upper() if i < len(atoms) else ""
                    if atp_u == "H" or atom_u.startswith("H"): continue
                    keep_idx.append(i)
                n_input_atoms_before = len(atmtp)
                n_removed_h = n_input_atoms_before - len(keep_idx)
                if n_removed_h > 0:
                    geometry = [geometry[i] for i in keep_idx]
                    atoms = [atoms[i] for i in keep_idx]
                    atmtp = [atmtp[i] for i in keep_idx]
                    residues = [residues[i] for i in keep_idx]
                    resnumbr = [resnumbr[i] for i in keep_idx]
                    chainid = [chainid[i] for i in keep_idx]
                    restp = [restp[i] for i in keep_idx]
                    seqidcodes = [seqidcodes[i] for i in keep_idx]
                    pdblines = [pdblines[i] for i in keep_idx]
                    self.logger.info("Stripped %d input H atoms from structure (use --keep-input-h to disable).",n_removed_h)
            except Exception as e:
                self.logger.warning("Failed to strip input hydrogens: %s",e)
        #strip *external/artificial* caps already present in the input structure (e.g., Maestro-like ACE/NMA/NME)
        #these are not experimentally resolved and should never be part of the IPA analysis.
        if not getattr(self.config,"keep_external_caps",False):
            try:
                keep_idx = []
                removed_atoms = 0
                removed_res_keys = set()
                for i in range(len(residues)):
                    rname = str(residues[i]).strip().upper()
                    key = (str(chainid[i]).strip(),int(resnumbr[i]))
                    if AA.IsExternalCap(rname):
                        removed_atoms += 1
                        removed_res_keys.add((rname,key[0],key[1]))
                        continue
                    keep_idx.append(i)
                if removed_atoms > 0:
                    geometry = [geometry[i] for i in keep_idx]
                    atoms = [atoms[i] for i in keep_idx]
                    atmtp = [atmtp[i] for i in keep_idx]
                    residues = [residues[i] for i in keep_idx]
                    resnumbr = [resnumbr[i] for i in keep_idx]
                    chainid = [chainid[i] for i in keep_idx]
                    restp = [restp[i] for i in keep_idx]
                    seqidcodes = [seqidcodes[i] for i in keep_idx]
                    pdblines = [pdblines[i] for i in keep_idx]
                    self.logger.warning("Stripped %d atoms belonging to %d external/artificial cap residues from input structure.",removed_atoms,len(removed_res_keys))
                    #extra user-facing summary: unique residue names removed
                    unique_caps = sorted({r for (r,_,_) in removed_res_keys})
                    self.logger.warning("External caps removed (unique resnames): %s",", ".join(unique_caps) if unique_caps else "<?>")
                    if getattr(self.config,"debug",False):
                        removed_list = sorted(list(removed_res_keys),key = lambda x: (x[1],x[2],x[0]))
                        for rname,ch,rn in removed_list:
                            self.logger.debug("Removed external cap residue: %s %s %s",rname,ch,rn)
            except Exception as e:
                self.logger.warning("Failed to strip external caps from structure: %s",e)
        #determine termini information
        terminal = ["N"]
        for idx in range(1,len(resnumbr)):
            if (resnumbr[idx] != resnumbr[idx - 1]) and (resnumbr[idx - 1] + 1 != resnumbr[idx]):
                terminal.append("N")
            else:
                terminal.append("0")
        for idx in range(1,len(resnumbr)):
            if terminal[idx] == "N":
                terminal[idx - 1] = "C"
        if terminal and terminal[-1] == "0":
            terminal[-1] = "C"   #the file's last atom has no "next" residue to trigger the marking above
        #convert residue name to residue id if needed
        if (self.config.active_residue_name not in ("0",None,"")) and (self.config.active_residue in ("NAN",None)):
            for ires in range(len(residues)):
                if (residues[ires] == self.config.active_residue_name) and (self.config.active_residue_chain == chainid[ires]):
                    self.config.active_residue = resnumbr[ires]
                    break
        #separate active residue from protein
        active_residue, protein_residues = self._separate_active_residue(geometry,atoms,atmtp,residues,resnumbr,chainid,restp,seqidcodes,pdblines,terminal)
        if not active_residue.atoms:
            raise InPocketAnalyzerError(f"Active residue {self.config.active_residue} in chain {self.config.active_residue_chain} not found")
        self.logger.info(f"Found active residue: {active_residue.name}")
        self.logger.info(f"Found {len(protein_residues)} protein residues")
        if AA.IsBaseAminoAcid(self.config.active_residue_name) or AA.IsExtendedAminoAcid(self.config.active_residue_name) or AA.IsPostTranslationalModification(self.config.active_residue_name):
            self.config.active_residue_is_like_aa = True
            #in this case we must fetch additional information from previous and following residues
            found_N_term_connection = False
            cannot_found_C_term_connection = False
            for atom in active_residue.atom_types:
                if atom == "OXT": cannot_found_C_term_connection = True
                elif atom == "OT2": cannot_found_C_term_connection = True
            for iresidue in range(len(protein_residues)):
                residue = protein_residues[iresidue]
                if (residue.number == active_residue.number - 1) and (residue.chain_id == active_residue.chain_id):
                    previous_residue = protein_residues[iresidue]
                    #this logic is needed because of seqid codes
                    if iresidue < len(protein_residues) - 1:
                        next_residue = protein_residues[iresidue + 1]
                        if (previous_residue.number == next_residue.number) and (previous_residue.seqidcode != next_residue.seqidcode):
                            #in this case we need to go to the last seqid code
                            for jresidue in range(iresidue + 1,len(protein_residues)):
                                next_residue = protein_residues[jresidue]
                                if previous_residue.number != next_residue.number:
                                    previous_residue = protein_residues[jresidue - 1]
                                    break
                    if AA.IsBaseAminoAcid(previous_residue.name) or AA.IsExtendedAminoAcid(previous_residue.name) or AA.IsPostTranslationalModification(previous_residue.name):
                        for iadtm in range(len(previous_residue.atoms)):
                            if (previous_residue.atom_types[iadtm] == "C") or (previous_residue.atom_types[iadtm] == "O"): 
                                active_residue.geometry.append(previous_residue.geometry[iadtm])
                                active_residue.atoms.append(previous_residue.atoms[iadtm])
                                active_residue.atom_types.append(previous_residue.atom_types[iadtm] + "p")
                                active_residue.pdb_lines.append(previous_residue.pdb_lines[iadtm])
                                found_N_term_connection = True
                    break
            active_residue.terminal = "0"
            if cannot_found_C_term_connection: active_residue.terminal += "C-term"
            if found_N_term_connection: active_residue.terminal += "N-term"
        return active_residue,protein_residues

    def _separate_active_residue(self,geometry,atoms,atmtp,residues,resnumbr,chainid,restp,seqidcodes,pdblines,terminal) -> Tuple[RES.ResidueData, List[RES.ResidueData]]:
        """Separate the active residue from the rest of the protein."""
        active_data = RES.ResidueData([],[],[],"",0,"","",[])
        protein_residues_dict = {}
        #convert active residue if needed
        if isinstance(self.config.active_residue,str):
            for i, (res,chain,num) in enumerate(zip(residues,chainid,resnumbr)):
                if ((res == self.config.active_residue) and (chain == self.config.active_residue_chain)):
                    self.config.active_residue = num
                    self.config.active_residue_name = res
                    break
        if (self.config.active_residue_name is None) or (self.config.active_residue_name == "None"):
            for i,(res,chain,num) in enumerate(zip(residues,chainid,resnumbr)):
                if ((num == self.config.active_residue) and (chain == self.config.active_residue_chain)):
                    self.config.active_residue_name = res
                    break
        #separate residues
        current_residue_idx = -1
        prev_resnum = None
        prev_chainid = None
        prev_seqcode = "0"
        for i in range(len(atoms)):
            if (resnumbr[i] == self.config.active_residue and chainid[i] == self.config.active_residue_chain):
                #active residue
                active_data.geometry.append(geometry[i])
                active_data.atoms.append(atoms[i])
                active_data.atom_types.append(atmtp[i])
                active_data.name = residues[i]
                active_data.number = resnumbr[i]
                active_data.chain_id = chainid[i]
                active_data.residue_type = restp[i]
                active_data.seqidcode = seqidcodes[i]
                active_data.pdb_lines.append(pdblines[i])
                active_data.terminal = terminal[i]
            else:
                #protein residue
                if (prev_resnum != resnumbr[i]) or (prev_seqcode != seqidcodes[i]) or (prev_chainid != chainid[i]):
                    current_residue_idx += 1
                    prev_resnum = resnumbr[i]
                    prev_seqcode = seqidcodes[i]
                    prev_chainid = chainid[i]
                if current_residue_idx not in protein_residues_dict:
                    protein_residues_dict[current_residue_idx] = RES.ResidueData([],[],[],residues[i],resnumbr[i],chainid[i],restp[i],[],terminal = terminal[i],seqidcode = seqidcodes[i])
                elif terminal[i] in ("N","C"):
                    protein_residues_dict[current_residue_idx].terminal = terminal[i]
                protein_residues_dict[current_residue_idx].geometry.append(geometry[i])
                protein_residues_dict[current_residue_idx].atoms.append(atoms[i])
                protein_residues_dict[current_residue_idx].atom_types.append(atmtp[i])
                protein_residues_dict[current_residue_idx].pdb_lines.append(pdblines[i])
        return active_data,list(protein_residues_dict.values())

    def process_ligand_states(self, active_residue: RES.ResidueData) -> List[LigandState]:
        """Process ligand protonation states and tautomers."""
        self.logger.info("Processing ligand protonation states")
        #write base ligand geometry
        ligand_name = f"{self.config.active_residue_chain}{self.config.active_residue}"
        base_file = (self.config.output_folder / self.config.geometry_folder / ligand_name)
        sdf_file_is_available = False
        try:
            GEOM.Write2XYZ(active_residue.atoms,active_residue.geometry,str(base_file))
            if self.config.active_residue_is_like_aa:
                connectivity = AA.AminoAcidConnectivity(active_residue.name,active_residue.terminal)
                local_charges = AA.GetLocalCharges(active_residue.name,active_residue.atom_types)
                output_string = GEOM.WriteSDF(active_residue.atom_types,active_residue.atoms,active_residue.geometry,connectivity,local_charges,str(base_file))
                sdf_file_is_available = True
            else: 
                if self.config.mmcif_file:
                    #get the ligand's sdf from cif
                    GEOM.WriteGemmiPDB(active_residue.geometry,
                                       active_residue.atoms,
                                       active_residue.atom_types,
                                       [active_residue.name]*len(active_residue.atoms),
                                       [active_residue.number]*len(active_residue.atoms),
                                       [active_residue.chain_id]*len(active_residue.atoms),
                                       [active_residue.residue_type]*len(active_residue.atoms),
                                       str(base_file))
                    connectivity = MONOMER.GeneralisedExtractCCDconnectivityGEMMI(self.config.pdb_file,active_residue.name)
                    output_string = GEOM.WriteSDF(active_residue.atom_types,active_residue.atoms,active_residue.geometry,connectivity,[],str(base_file))
                    sdf_file_is_available = True
                    self.logger.warning(output_string)
                else:
                    #search for ligand in the monomer library
                    connectivity = []
                    if MonomerExistsInLibrary(active_residue.name):
                        path2monomerfile = GetMonomerPath(active_residue.name)
                        connectivity = MONOMER.ExtractCCDconnectivityGEMMI(path2monomerfile)
                        output_string = GEOM.WriteSDF(active_residue.atom_types,active_residue.atoms,active_residue.geometry,connectivity,[],str(base_file))
                        sdf_file_is_available = True
                    GEOM.Write2PDB(active_residue.pdb_lines,str(base_file))
        except Exception as e:
            raise InPocketAnalyzerError(f"Failed to write ligand geometry: {e}")
        if self.config.protonation_mode == ProtonationMode.ORIGINAL and sdf_file_is_available:
            SANITY.repair_native_sdf_h_connectivity(Path(str(base_file) + ".sdf"),logger = self.logger)
        ligand_states = []
        #check whether a ligand topology is given by the user
        the_ligand_file = ""
        use_user_ligand_topology = (self.config.protonation_mode != ProtonationMode.ORIGINAL and self.config.ligand_topology not in (None,""))
        if use_user_ligand_topology:
            the_ligand_file = str(base_file) + ".sdf"
            shutil.copy(self.config.ligand_topology,the_ligand_file)
        elif sdf_file_is_available:
            the_ligand_file = str(base_file) + ".sdf"
        else:
            the_ligand_file = str(base_file) + ".pdb"
        if self.config.protonation_mode != ProtonationMode.ORIGINAL:
            #central input sanitization for ligand (backend-independent)
            the_ligand_file = SANITY.InputSanitizer.sanitize_ligand_file_for_protonation(the_ligand_file,
                                                                                         output_folder = self.config.output_folder / self.config.geometry_folder,
                                                                                         keep_input_h = self.config.keep_input_h,
                                                                                         logger = self.logger)
        elif str(the_ligand_file).lower().endswith(".sdf"):
            #original/native mode preserves the input chemistry, but we still repair missing formal charges.
            SANITY.InputSanitizer.repair_ligand_sdf_formal_charges(Path(the_ligand_file),logger = self.logger)
        if (((self.config.protonation_mode == ProtonationMode.VSEPR) or (self.config.protonation_mode == ProtonationMode.ORIGINAL)) and self.config.active_residue_is_like_aa):
            for iline in range(len(active_residue.pdb_lines)):
                line = active_residue.pdb_lines[iline].replace("HETATM","ATOM  ")
                active_residue.pdb_lines[iline] = line
            if str(the_ligand_file).lower().endswith(".sdf"):
                pdb_for_vsepr = Path(str(the_ligand_file)).with_suffix(".pdb")
                if not pdb_for_vsepr.exists():
                    pdb_lines = active_residue.pdb_lines
                    if not self.config.keep_input_h:
                        pdb_lines, n_removed = SANITY.InputSanitizer.strip_h_from_pdb_lines(pdb_lines)
                        if n_removed:
                            self.logger.debug("Stripped %d input H atoms from active-residue PDB template for VSEPR: %s",n_removed,pdb_for_vsepr.name)
                    GEOM.Write2PDB(pdb_lines,str(pdb_for_vsepr.with_suffix("")))
        try:
            if self.config.protonation_mode == ProtonationMode.ORIGINAL:
                ligand_states.append(self._prepare_original_mode_ligand_state(active_residue,base_file,str(the_ligand_file)))
            elif self.config.protonation_mode == ProtonationMode.RDKIT:
                files,charges,smiles = self._rdkit_protonation(str(the_ligand_file))
                if self.config.exclude_rdkit_bad_structures:
                    kept,rejected = SC.clean_rdkit_tautomer_smiles(smiles)
                    files_new = []
                    charges_new = []
                    smiles_new = []
                    for istate in range(len(smiles)):
                        smile = smiles[istate]
                        if smile in kept:
                            self.logger.warning(f"keeping state {istate + 1}")
                            files_new.append(files[istate])
                            charges_new.append(charges[istate])
                            smiles_new.append(smile)
                        else:
                            for smi,reason in rejected:
                                if smi == smile:
                                    self.logger.warning(f"rejected state {istate + 1} because: {reason}")
                                    break
                    files = files_new
                    charges = charges_new
                    smiles = smiles_new
                for i,(file, charge, smile) in enumerate(zip(files,charges,smiles)):
                    ligand_states.append(LigandState(Path(file),charge,i,smile))
            elif self.config.protonation_mode == ProtonationMode.OBABEL:
                new_basefilename,charge = self._obabel_protonation(str(the_ligand_file),self.config.charge,self.config.user_defined_charge)
                file = new_basefilename
                ligand_states.append(LigandState(file_path = Path(file),charge = charge,ligand_index = 0,smiles = "smile"))
            elif self.config.protonation_mode == ProtonationMode.VSEPR:
                charge,prot_file = self._vsepr_protonation(str(the_ligand_file),self.config.charge,self.config.user_defined_charge)
                ligand_states.append(LigandState(file_path = Path(prot_file),charge = charge,ligand_index = 0))
            elif self.config.protonation_mode == ProtonationMode.CCD:
                charge = self._ccd_protonation(str(base_file),self.config.charge,self.config.user_defined_charge)
                ligand_states.append(LigandState(file_path = Path(str(base_file) + ".sdf"),charge = charge,ligand_index = 0))
        except Exception as e:
            raise InPocketAnalyzerError(f"Protonation failed: {e}")
        if not ligand_states:
            raise InPocketAnalyzerError("No ligand states generated")
        #process tautomers if requested
        if self.config.tautomerise and self.config.protonation_mode != ProtonationMode.ORIGINAL:
            try:
                ligand_states = self._process_tautomers(ligand_states)
            except Exception as e:
                self.logger.warning(f"Tautomer processing failed: {e}")
        elif self.config.tautomerise and self.config.protonation_mode == ProtonationMode.ORIGINAL:
            self.logger.info("Skipping tautomer generation in original/native protonation mode")
        self._validate_prepared_ligand_states(ligand_states)
        #clean the SMILES
        for state in ligand_states:
            smile = (state.smiles.replace("C([H])","C").replace("c([H])","c").replace("C[H]","C").replace("c[H]","c").replace("C([H])","C").replace("c([H])","c").replace("C[H]","C").replace("c[H]","c"))
            state.smiles = smile
        #user selection
        if len(ligand_states) > 1:
            self.logger.warning(f"There are {len(ligand_states)} ligand chemical states generated.")
            #in non-interactive mode, default to all
            #in interactive mode (when sys.stdin.isatty()), prompt user
            if sys.stdin.isatty():
                self.logger.info(f"WARNING: list of ligand chemical states:")
                for state in ligand_states:
                    self.logger.info(state.smiles)
                codes2do = input("please select 'first', 'all', or comma-separated input list: ")
                if codes2do.lower() != "all":
                    temp_ligand_files = copy.deepcopy(ligand_states)
                    ligand_states = []
                    if codes2do.lower() == "first":
                        ligand_states.append(temp_ligand_files[0])
                    else:
                        completelist = codes2do.split(",")
                        for idl in completelist:
                            if "-" in idl:
                                newsplit = idl.split("-")
                                for idx in range(int(newsplit[0]), int(newsplit[1]) + 1):
                                    index = idx - 1
                                    if 0 <= index < len(temp_ligand_files):
                                        ligand_states.append(temp_ligand_files[index])
                            else:
                                index = int(idl) - 1
                                if 0 <= index < len(temp_ligand_files):
                                    ligand_states.append(temp_ligand_files[index])
            else:
                self.logger.info("Non-interactive mode: using all ligand states")
        #cache the *selected* ligand states as stable files in geometries/ligand_states/
        self._cache_ligand_states(ligand_states)
        self.logger.info(f"Generated {len(ligand_states)} ligand states")
        return ligand_states

    def _cache_ligand_states(self, ligand_states: List[LigandState]) -> None:
        """Stable cache:
           Copy the *selected* ligand state input files into geometries/ligand_states/.
           This folder behaves as a cache in normal and --best-debug modes:
           - it is never the "working" location
           - it must survive --best cleanup"""
        #in compact --best (no --best-debug), do NOT create the ligand_states cache folder.
        #we only want best_* outputs in the root output folder.
        if self.config.best_mode and (not self.config.best_debug):
            return
        cache_dir = self.config.output_folder / self.config.geometry_folder / "ligand_states"
        cache_dir.mkdir(parents = True,exist_ok = True)
        for st in ligand_states:
            try:
                src = Path(st.file_path)
                if not src.exists(): continue
                dst = cache_dir / src.name
                #if already cached, do nothing
                if dst.exists(): continue
                shutil.copy2(src,dst)
            except Exception as e:
                self.logger.debug(f"Failed to cache ligand state {st.file_path}: {e}")

    def _rdkit_protonation(self, base_file: str) -> Tuple[List[str], List[int], List[str]]:
            """Perform RDKit-based protonation."""
            Protoro.print_header_protoro()
            image_mode = "none"
            if getattr(self.config,"keep_png",False):
                image_mode = "highres" if getattr(self.config,"keep_png_hq",False) else "lowres"
            #only used if protoro fallback triggers; protoro will append when needed
            warning_log_path = str(self.config.output_folder / "WARNING_ligand_input.log")
            return Protoro.RDKitProtonation(base_file,self.config.ph,image_mode,logger = self.logger,warning_log_path = warning_log_path)

    def _obabel_protonation(self, base_file: str, user_charge, userdefined_charge) -> List[LigandState]:
        """Perform OpenBabel-based protonation."""
        charge = GEOM.ObabelProtonation(base_file,self.config.ph)
        new_basefilename = base_file.replace(".sdf","H.mol2").replace(".pdb","H.mol2")
        GEOM.MOL2toXYZ(new_basefilename,"ligand")
        base_file = new_basefilename.replace(".mol2",".xyz")
        charge = self._resolve_active_charge_with_user_override(charge,user_charge,userdefined_charge,"OpenBabel",base_file)
        return base_file,charge

    def _resolve_active_charge_with_user_override(self,backend_charge: int,user_charge: int,userdefined_charge: bool,backend_name: str,structure_path: str) -> int:
        """Resolve active-moiety charge when a backend charge is available."""
        raw_backend_charge = backend_charge        
        try:
            backend_charge = int(backend_charge)
        except Exception:
            self.logger.warning("Failed to parse backend-defined active-moiety charge from %s "
                                "(raw value: %r). Falling back to 0.",
                                backend_name,
                                raw_backend_charge)
            backend_charge = 0
        if not userdefined_charge:
            return backend_charge
        try: user_charge = int(user_charge)
        except Exception: user_charge = 0
        if backend_charge == user_charge: return user_charge
        if getattr(self.config,"force_user_charge",False):
            self.logger.warning("%s and user defined active-moiety charge mismatch; "
                                "using user defined charge %d (--force-user-charge enabled, backend=%d).",
                                backend_name,user_charge,backend_charge)
            return user_charge
        print(f"WARNING: {backend_name} and user defined charge mismatch")
        print("                                   user defined:",user_charge)
        print(f"                          {backend_name} defined:",backend_charge)
        return int(input("please modify the file "
                         + str(structure_path)
                         + " to the correct protonation state and redefine the charge: "))

    def _vsepr_protonation(self, base_file: str, user_charge = None, userdefined_charge: bool = False) -> Tuple[int, str]:
        """Perform VSEPR-based protonation."""
        backend_name = "VSEPR"
        input_path = str(base_file)
        if self.config.ligand_topology not in (None,""):
            self.logger.info("VSEPR ligand protonation: using prepared working copy %s derived from user topology %s.",
                             Path(input_path).name,
                             Path(str(self.config.ligand_topology)).name)
        lower = input_path.lower()
        if lower.endswith(".pdb"):
            pdb_input = input_path
            output_base = input_path[:-4]
        elif lower.endswith(".sdf"):
            output_base = str(Path(input_path).with_suffix(""))
            pdb_candidates = [str(Path(input_path).with_suffix(".pdb"))]
            if output_base.endswith("_noH"):
                pdb_candidates.append(output_base[:-4] + ".pdb")
            pdb_input = next((cand for cand in pdb_candidates if Path(cand).exists()),None)
            if pdb_input is None:
                raise InPocketAnalyzerError("VSEPR protonation requires a PDB template for SDF input; "
                                            f"tried: {', '.join(pdb_candidates)}")
        else:
            pdb_input = input_path + ".pdb"
            output_base = input_path
        if self.config.active_residue_is_like_aa:
            try:
                charge = VSEPR.ProteinProtonation(pdb_input,output_base + "H")
                structure_path = output_base + "H.xyz"
                backend_name = "ProteinProtonation (VSEPR active-residue path)"
            except Exception as e:
                self.logger.warning("ProteinProtonation failed for active residue %s in VSEPR mode (%s). "
                                    "Falling back to generic VSEPR protonation.",
                                    self.config.active_residue_name,e)
                charge = VSEPR.VSEPRProtonation(pdb_input,output_base)
                structure_path = output_base + ".xyz"
                backend_name = "VSEPR (active-residue fallback)"
        elif lower.endswith(".sdf"):
            charge = VSEPR.VSEPRSDFMOL2Protonation(input_path,input_path.replace(".sdf","_H"))
            structure_path = input_path.replace(".sdf","_H") + ".xyz"
        else:
            charge = VSEPR.VSEPRProtonation(pdb_input,output_base)
            structure_path = output_base + ".xyz"
        charge = self._resolve_active_charge_with_user_override(charge,user_charge,userdefined_charge,backend_name,structure_path)
        return charge,structure_path

    def _calculate_ligand_charge_original_mode(self,ligand_file: Optional[str] = None) -> Optional[int]:
        """Infer ligand charge for original/native mode using VSEPR as a charge probe."""
        if not ligand_file: return None
        ligand_path = Path(ligand_file)
        probe_input_xyz = (ligand_path
                           if ligand_path.suffix.lower() == ".xyz"
                           else ligand_path.with_suffix(".xyz"))
        if not probe_input_xyz.exists():
            self.logger.info("Original/native ligand VSEPR charge probe skipped for %s: XYZ input not found (%s).",
                             ligand_path.name,
                             probe_input_xyz.name)
            return None
        if not getattr(self.config,"keep_input_h",False):
            self.logger.warning("Original/native ligand VSEPR charge probe for %s is running with keep_input_h=False; "
                                "input hydrogens may already have been stripped before charge inference.",
                                ligand_path.name)
        probe_base = probe_input_xyz.with_name(probe_input_xyz.stem + "_vsepr_chargeprobe")
        try:
            adjusted_charge,probe_meta = VSEPR.VSEPRChargeProbeXYZ(str(probe_input_xyz),str(probe_base))
        except Exception as e:
            self.logger.warning("Original/native ligand VSEPR charge probe failed for %s: %s",ligand_path.name,e)
            return None
        if adjusted_charge is None:
            if probe_meta.get("heavy_atoms_unchanged") is False:
                self.logger.warning("Original/native ligand VSEPR charge probe rejected for %s: heavy-atom framework changed "
                                    "(input heavy = %s output heavy = %s). Falling back to SDF charge.",
                                    ligand_path.name,
                                    probe_meta.get("input_heavy_atoms","?"),
                                    probe_meta.get("output_heavy_atoms","?"))
            else:
                self.logger.warning("Original/native ligand VSEPR charge probe did not return a usable charge for %s. "
                                    "stdout = %r stderr = %r",
                                    ligand_path.name,
                                    probe_meta.get("stdout",""),
                                    probe_meta.get("stderr",""))
            return None
        raw_charge = int(probe_meta.get("charge_vsepr",adjusted_charge))
        delta_h = int(probe_meta.get("delta_h",0))
        input_h = int(probe_meta.get("input_h",0))
        output_h = int(probe_meta.get("output_h",0))
        if delta_h == 0:
            self.logger.info("Original/native ligand VSEPR charge probe for %s: raw_charge = %+d input_H = %d output_H = %d "
                             "(no H correction needed).",ligand_path.name,raw_charge,input_h,output_h)
        else:
            self.logger.warning("Original/native ligand VSEPR charge probe for %s changed proton count: "
                                "raw_charge=%+d input_H=%d output_H=%d delta_H=%+d -> adjusted_input_charge=%+d",
                                ligand_path.name,raw_charge,input_h,output_h,delta_h,int(adjusted_charge))
        return int(adjusted_charge)

    def _resolve_original_mode_ligand_charge(self,ligand_file: Optional[str] = None) -> int:
        """Resolve the active-moiety charge for original/native mode.
           Current fallback order:
             1) user-defined charge (--charge / --force-user-charge workflow)
             2) VSEPR ligand charge probe (charge-only, with H-count correction)
             3) total formal charge already encoded in the ligand SDF
             4) 0"""
        if getattr(self.config,"user_defined_charge",False):
            try:
                user_charge = int(self.config.charge)
                self.logger.info("Original/native ligand charge: using user-defined charge %+d for %s.",
                                 user_charge,
                                 Path(ligand_file).name if ligand_file else "<unknown>")
                return user_charge
            except Exception:
                self.logger.warning("Failed to parse user-defined ligand charge %r in original/native mode; using 0.",getattr(self.config,"charge",None))
                return 0
        auto_charge = self._calculate_ligand_charge_original_mode(ligand_file)
        if auto_charge is not None:
            if ligand_file and str(ligand_file).lower().endswith(".sdf"):
                sdf_charge = SANITY.InputSanitizer.get_sdf_total_formal_charge(Path(ligand_file),logger = self.logger)
                if (sdf_charge is not None) and (int(sdf_charge) != int(auto_charge)):
                    self.logger.warning("Original/native ligand charge mismatch for %s: VSEPR probe = %+d SDF formal charge = %+d. "
                                        "Using VSEPR probe charge.",Path(ligand_file).name,int(auto_charge),int(sdf_charge))
            return int(auto_charge)
        if ligand_file and str(ligand_file).lower().endswith(".sdf"):
            sdf_charge = SANITY.InputSanitizer.get_sdf_total_formal_charge(Path(ligand_file),logger = self.logger)
            if sdf_charge is not None:
                self.logger.info("Original/native ligand charge fallback for %s: using SDF formal charge %+d.",
                                 Path(ligand_file).name,int(sdf_charge))
                return int(sdf_charge)
        self.logger.warning("Original/native ligand charge fallback for %s: no user/VSEPR/SDF charge available, using 0.",
                            Path(ligand_file).name if ligand_file else "<unknown>")
        return 0

    def _prepare_original_mode_ligand_state(self,active_residue: RES.ResidueData,base_file: Path,the_ligand_file: str) -> LigandState:
        """Prepare exactly one ligand state while preserving the input chemistry."""
        if self.config.active_residue_is_like_aa:
            pdb_path = Path(str(base_file) + ".pdb")
            if not pdb_path.exists():
                GEOM.Write2PDB(active_residue.pdb_lines,str(base_file))
            #keep the original heavy atoms/H atoms, then complete only the
            #newly exposed termini introduced by the residue capping context.
            #ProteinProtonation is the preferred path for amino-acid-like residues, but PTM/non-canonical names such as HIC can fail there.
            #in original/native mode if the template-driven protonator fails to interpret the residue name, we fall back to the generic VSEPR path.
            try:
                backend_charge = VSEPR.ProteinProtonation(str(pdb_path),str(base_file) + "H")
                prot_file = Path(str(base_file) + "H.xyz")
                backend_name = "ProteinProtonation (original/native)"
            except Exception as e:
                self.logger.warning("ProteinProtonation failed for active residue %s in original/native mode (%s). "
                                    "Falling back to generic VSEPR completion for cap/terminal hydrogens.",active_residue.name,e)
                backend_charge = VSEPR.VSEPRProtonation(str(pdb_path),str(base_file))
                prot_file = Path(str(base_file) + ".xyz")
                backend_name = "VSEPR (original/native fallback)"
            charge = self._resolve_active_charge_with_user_override(backend_charge,self.config.charge,self.config.user_defined_charge,backend_name,str(prot_file))
            return LigandState(file_path = prot_file,charge = charge,ligand_index = 0)
        charge = self._resolve_original_mode_ligand_charge(the_ligand_file)
        ligand_state_file = Path(the_ligand_file)
        if ligand_state_file.suffix.lower() == ".pdb":
            ligand_state_file = Path(str(base_file) + ".xyz")
        return LigandState(file_path = ligand_state_file,charge = charge,ligand_index = 0)

    def _ccd_protonation(self, base_file: str, user_charge: int, userdefined_charge: bool) -> int:
        """Perform CCD-based protonation."""
        if not self.config.ccd_file:
            raise InPocketAnalyzerError("CCD file required for CCD protonation")
        charge = VSEPR.CCDProtonation(base_file + ".pdb",base_file + ".pdb",self.config.ccd_file,base_file)
        charge = self._resolve_active_charge_with_user_override(charge,user_charge,userdefined_charge,"CCD",base_file + ".sdf")
        return charge

    def _process_tautomers(self, ligand_states: List[LigandState]) -> List[LigandState]:
        """Process tautomers for each ligand state.
        Behaviour:
        - Calls process_tautomers() which generates raw tautomers in a temporary folder.
        - Filters them by RMSD using _filter_unique_tautomers().
        - Writes curated tautomers to geometries/ligand_states/.
        - Removes the temporary <base_name>_tautomers directory."""
        unique_states: List[LigandState] = []
        geom_root = self.config.output_folder / self.config.geometry_folder
        for state in ligand_states:
            try:
                #generate raw tautomers into a temporary <base>_tautomers folder
                image_mode = "none"
                if getattr(self.config,"keep_png",False):
                    image_mode = "highres" if getattr(self.config,"keep_png_hq",False) else "lowres"
                tautomers,smiles = Tautoro.process_tautomers(str(state.file_path),str(geom_root),self.config.max_tautomers,image_mode = image_mode)
                if self.config.exclude_rdkit_bad_structures:
                    kept,rejected = SC.clean_rdkit_tautomer_smiles(smiles)
                    tautomers_new = []
                    smiles_new = []
                    for istate in range(len(smiles)):
                        smile = smiles[istate]
                        if smile in kept:
                            self.logger.warning(f"keeping state {istate + 1}")
                            tautomers_new.append(tautomers[istate])
                            smiles_new.append(smile)
                        else:
                            for smi,reason in rejected:
                                if smi == smile:
                                    self.logger.warning(f"rejected state {istate + 1} because: {reason}")
                                    break
                    tautomers = tautomers_new
                    smiles = smiles_new
                #no tautomers: keep original state
                if not tautomers:
                    unique_states.append(state)
                    continue
                #temporary tautomer directory created by process_tautomers
                tautomer_dir = Path(tautomers[0]).parent
                #filter unique tautomers based on RMSD
                filtered_tautomers,filtered_smiles = self._filter_unique_tautomers([str(f) for f in tautomers],smiles)
                #if filtering killed everything, keep the original state (but still try to clean tmp dir)
                if not filtered_tautomers:
                    unique_states.append(state)
                    if tautomer_dir.exists():
                        try:
                            shutil.rmtree(tautomer_dir)
                        except Exception as e:
                            self.logger.debug(f"Failed to remove temporary tautomer dir {tautomer_dir}: {e}")
                    continue
                #target directory for curated ligand states.
                #normal + --best-debug: stable cache in geometries/ligand_states/
                #compact --best: write to a temp folder so geometries/ligand_states/ is never created.
                if self.config.best_mode and (not self.config.best_debug): target_dir = geom_root / "_ligand_states_tmp"
                else: target_dir = geom_root / "ligand_states"
                target_dir.mkdir(parents = True,exist_ok = True)
                base_name = Path(state.file_path).stem  # e.g., "Z9999_prot_1"
                added_any = False
                for i,taut_file in enumerate(filtered_tautomers):
                    expected_name = f"{base_name}_tautomer_{i+1}.sdf"
                    new_path = target_dir / expected_name
                    try:
                        shutil.copy2(taut_file,new_path)
                        smile = filtered_smiles[i] if i < len(filtered_smiles) else ""
                        #keep positional arguments to match existing LigandState signature
                        unique_states.append(LigandState(Path(new_path),state.charge,0,smile))
                        added_any = True
                        self.logger.debug(f"Created tautomer: {new_path}")
                    except Exception as e:
                        self.logger.warning(f"Failed to process tautomer {taut_file}: {e}")
                #if nothing usable ended up being copied, fall back to original state
                if not added_any: unique_states.append(state)
                #if requested, preserve the RDKit PNG grid generated by Tautoro before deleting the temp dir
                if getattr(self.config, "keep_png", False) and tautomer_dir.exists():
                    try:
                        base_name_for_png = Path(state.file_path).stem
                        png_name = (f"{base_name_for_png}_tautomers_highres.png"
                                    if getattr(self.config, "keep_png_hq", False)
                                    else f"{base_name_for_png}_tautomers.png")
                        src_png = tautomer_dir / png_name
                        if src_png.exists():
                            png_dir = self.config.output_folder / "png"
                            png_dir.mkdir(parents = True,exist_ok = True)
                            shutil.move(str(src_png),str(png_dir / src_png.name))
                    except Exception as e:
                        self.logger.debug(f"Failed to preserve tautomer PNG from {tautomer_dir}: {e}")
                #clean up the temporary tautomer directory
                if tautomer_dir.exists():
                    try: shutil.rmtree(tautomer_dir)
                    except Exception as e: self.logger.debug(f"Failed to remove temporary tautomer dir {tautomer_dir}: {e}")
            except Exception as e:
                self.logger.warning(f"Failed to process tautomers for {state.file_path}: {e}")
                unique_states.append(state)
        return unique_states

    def _filter_unique_tautomers(self, tautomers: List[str], smiles: List[str]) -> Tuple[List[str], List[str]]:
        """Filter tautomers to keep only unique structures based on RMSD."""
        unique_tautomers = []
        unique_smiles = []
        for taut1,smile1 in zip(tautomers,smiles):
            is_unique = True
            for taut2 in unique_tautomers:
                try:
                    rmsd = RMSD.RMSDcalculation(taut1,taut2)
                    if rmsd < self.config.rmsd_tautomer_threshold:
                        is_unique = False
                        break
                except Exception as e:
                    self.logger.warning(f"RMSD calculation failed: {e}")
                    continue
            if is_unique:
                unique_tautomers.append(taut1)
                unique_smiles.append(smile1)
        return unique_tautomers,unique_smiles

    def _select_ligand_states(self, ligand_states: List[LigandState]) -> List[LigandState]:
        """Handle user selection of ligand states."""
        if len(ligand_states) <= 1 or not sys.stdin.isatty():
            return ligand_states
        self.logger.info(f"WARNING: there are {len(ligand_states)} ligand protonation states generated.")
        selection = input("Please select 'first', 'all', or input list (e.g., 1,3-5): ")
        if selection.lower() == "all":
            return ligand_states
        elif selection.lower() == "first":
            return [ligand_states[0]]
        else:
            return self._parse_selection(selection, ligand_states)

    def _parse_selection(self, selection: str, ligand_states: List[LigandState]) -> List[LigandState]:
        """Parse user selection string."""
        selected = []
        for part in selection.split(","):
            if "-" in part:
                start,end = map(int,part.split("-"))
                for i in range(start - 1,min(end,len(ligand_states))):
                    if 0 <= i < len(ligand_states): selected.append(ligand_states[i])
            else:
                idx = int(part) - 1
                if 0 <= idx < len(ligand_states):
                    selected.append(ligand_states[idx])
        return selected

    def find_neighbour_residues(self,active_residue: RES.ResidueData,protein_residues: List[RES.ResidueData],ligand_states: List[LigandState]) -> List[RES.ResidueData]:
        """Find residues within critical distance of the active residue."""
        self.logger.info(f"Finding neighbor residues within {self.config.critical_distance} Å")
        neighbors = []
        r2_critical = self.config.critical_distance**2
        covalent_threshold2 = self.config.covalent_threshold**2  #legacy fallback for ±1 residues when peptide-bond atoms cannot be identified
        #covalent-radius based threshold (more specific than a single fixed distance).
        #tune scale if needed (e.g. 1.2–1.35 depending on input quality).

        def _covalent_pair_threshold2(ei: str, ej: str) -> float:
            ri = AD._covalent_radius_angstrom(ei)
            rj = AD._covalent_radius_angstrom(ej)
            if ri is None or rj is None: thr = self.config.covalent_threshold
            else: thr = UC.covalent_scale*(ri + rj)
            return thr*thr
        
        def _atom_coord_by_type(residue: RES.ResidueData, atom_type: str):
            """Best-effort coordinate lookup by PDB-style atom_type (e.g. 'N', 'C')."""
            wanted = str(atom_type).strip().upper()
            atypes = getattr(residue,"atom_types",[]) or []
            geom = getattr(residue,"geometry",[]) or []
            if not atypes or not geom or (len(atypes) != len(geom)): return None
            for i,nm in enumerate(atypes):
                if str(nm).strip().upper() == wanted:
                    x,y,z = geom[i]
                    return (x,y,z)
            return None

        def _active_state_coord(atom_type: str, elem: str, max_map_dist: float = 1.5):
            """Map an atom from active_residue.geometry (by atom_type) to the active state
            geometry (actres_geometry) using nearest-by-coordinate matching on element."""
            ref = _atom_coord_by_type(active_residue,atom_type)
            if ref is None: return None
            idx = CONNECT._nearest_index_by_element_and_coord(actres_atoms,actres_geometry,str(elem).strip().upper(),ref,max_dist = max_map_dist)
            if idx is None:
                return ref
            x,y,z = actres_geometry[idx]
            return (x,y,z)
        
        #read active residue geometry (use first ligand state)
        try:
            actres_atoms, actres_geometry = GEOM.ReadGeometry(str(ligand_states[0].file_path))
        except Exception as e:
            self.logger.warning(f"Failed to read active residue geometry: {e}")
            actres_atoms = active_residue.atoms
            actres_geometry = active_residue.geometry
        #use a set for O(1) membership and to avoid duplicates
        #include seqidcode when present to avoid collisions on insertion codes
        residues2skip = set()
        #number of residues around the active residue to skip when checking for covalent linkages.
        covalent_skip_window = 1       #default 1
        self.logger.info("List of residues near the ligand and their distances:")
        for protein_res in protein_residues:
            #skip if already in skip list
            res_key = [protein_res.chain_id,protein_res.number]
            res_key = (protein_res.chain_id,protein_res.number,getattr(protein_res,"seqidcode","0"))
            if res_key in residues2skip: continue
            #skip water if requested
            if self.config.skip_water and protein_res.name.upper() in ["HOH","SOL","WAT","H2O","TP3","TIP"]:
                continue
            #check for covalent linkage
            if (protein_res.chain_id == self.config.active_residue_chain
                and abs(protein_res.number - self.config.active_residue) <= covalent_skip_window
                and protein_res.number != self.config.active_residue):
                #(+/-1) covalent skip:
                #for peptide-like ligands, only check the peptide-bond atoms (C-N).
                #fall back to legacy all-atom covalent check only if required atoms are not available.
                did_specific_check = False
                if self.config.active_residue_is_like_aa:
                    try:
                        aN = _active_state_coord("N","N")
                        aC = _active_state_coord("C","C")
                        if protein_res.number == (self.config.active_residue - 1):
                            #previous residue: C(prev) -- N(active)
                            pC = _atom_coord_by_type(protein_res,"C")
                            if pC is not None and aN is not None:
                                did_specific_check = True
                                d2 = sum((float(a) - float(p))**2 for a,p in zip(aN,pC))
                                if d2 < _covalent_pair_threshold2("N","C"):
                                    residues2skip.add(res_key)
                                    self.logger.info(f"Skipping {protein_res.name} {protein_res.chain_id} {protein_res.number} {np.sqrt(d2):.2f} Å (peptide bond)")
                                    continue
                        elif protein_res.number == (self.config.active_residue + 1):
                            #next residue: C(active) -- N(next)
                            pN = _atom_coord_by_type(protein_res,"N")
                            if pN is not None and aC is not None:
                                did_specific_check = True
                                d2 = sum((float(a) - float(p))**2 for a,p in zip(aC,pN))
                                if d2 < _covalent_pair_threshold2("C","N"):
                                    residues2skip.add(res_key)
                                    self.logger.info(f"Skipping {protein_res.name} {protein_res.chain_id} {protein_res.number} {np.sqrt(d2):.2f} Å (peptide bond)")
                                    continue
                    except Exception:
                        did_specific_check = False
                if not did_specific_check:
                    #legacy behaviour (all heavy atoms, fixed threshold)
                    mindist2 = float("inf")
                    for p_geom in protein_res.geometry:
                        for a_idx,a_geom in enumerate(actres_geometry):
                            if actres_atoms[a_idx] == "H": continue  #skip protons
                            dist2 = sum((a - p)**2 for a,p in zip(a_geom,p_geom))
                            mindist2 = min(mindist2,dist2)
                            if mindist2 < covalent_threshold2: break
                        if mindist2 < covalent_threshold2: break
                    if mindist2 < covalent_threshold2:
                        residues2skip.add(res_key)
                        self.logger.warning(f"Skipping {protein_res.name} {protein_res.chain_id} {protein_res.number} {np.sqrt(mindist2):.2f} Å (within covalent range)")
                        continue
            #check distance
            mindist2 = float("inf")
            covalent_hit = False
            covalent_d2 = float("inf")
            p_atoms = list(getattr(protein_res,"atoms",[]) or [])
            for a_idx,a_geom in enumerate(actres_geometry):
                ai = str(actres_atoms[a_idx]).strip().upper()
                if ai == "H": continue  #skip protons
                for p_idx,p_geom in enumerate(protein_res.geometry):
                    dist2 = sum((a - p)**2 for a,p in zip(a_geom,p_geom))
                    if dist2 < mindist2:
                        mindist2 = dist2
                    #second-stage filter: exclude residues covalently bound to the ligand
                    if not covalent_hit:
                        pj = str(p_atoms[p_idx]).strip().upper() if p_idx < len(p_atoms) else ""
                        if pj and pj != "H":
                            if dist2 < _covalent_pair_threshold2(ai,pj):
                                if ION.IsMetal(pj)[0]: 
                                    self.logger.warning(f"Metal-containing system found for {protein_res.name} {protein_res.chain_id} {protein_res.number}: interaction energies with isolated metals are unreliable")
                                    continue
                                else:
                                    covalent_hit = True
                                    covalent_d2 = dist2
                    else:
                        #track best covalent distance for logging
                        if dist2 < covalent_d2:
                            pj = str(p_atoms[p_idx]).strip().upper() if p_idx < len(p_atoms) else ""
                            if pj and pj != "H" and dist2 < _covalent_pair_threshold2(ai,pj):
                                covalent_d2 = dist2
            if mindist2 < r2_critical:
                if covalent_hit:
                    residues2skip.add(res_key)
                    self.logger.warning(f"Skipping {protein_res.name} {protein_res.chain_id} {protein_res.number} {np.sqrt(covalent_d2):.2f} Å (within covalent range)")
                    continue
                neighbors.append(protein_res)
                distance = np.sqrt(mindist2)
                self.logger.info(f"    {protein_res.name} {protein_res.chain_id} {protein_res.number} {distance:.2f} Å")
        self.logger.info(f"Found {len(neighbors)} neighbor residues")
        #create caps for neighbors
        if self.config.cap_residues:
            self._create_caps_for_neighbors(neighbors,protein_residues)
        if self.config.flip_residues and self.config.protonation_mode != ProtonationMode.ORIGINAL:
            flip_residues = []
            for residue in neighbors:
                should_flip = ((residue.name == "ASN" and self.config.flip_ASN) or (residue.name == "GLN" and self.config.flip_GLN) or (residue.name == "HIS" and self.config.flip_HIS))
                if should_flip:
                    flipped_residue = self._flip_residue(residue)
                    flipped_residue.name = residue.name + "_flip"        #mark as flipped
                    flipped_residue.residue_id += "f"                    #mark as flipped
                    flip_residues.append(flipped_residue)
            for residue in flip_residues:
                neighbors.append(residue)
        return neighbors

    def _create_caps_for_neighbors(self, neighbours: List[RES.ResidueData], all_protein_residues: List[RES.ResidueData]) -> None:
        """Create caps for neighbour residues using full protein context."""
        index = 0
        for neighbour in neighbours:
            #skip if not an amino acid or is N-terminal
            not_aa = (not AA.IsBaseAminoAcid(neighbour.name)) and (not AA.IsExtendedAminoAcid(neighbour.name)) and (not AA.IsPostTranslationalModification(neighbour.name))
            not_na = (not NA.IsBaseNucleicAcid(neighbour.name)) and (not NA.IsExtendedNucleicAcid(neighbour.name))
            if not_aa and not_na:
                continue
            if neighbour.terminal == "N":
                continue
            #find the PREVIOUS residue in the sequence
            idx_this_res = -1
            for ires in range(index,len(all_protein_residues)):
                residue = all_protein_residues[ires]
                if (residue.name == neighbour.name) and (residue.chain_id == neighbour.chain_id) and (residue.number == neighbour.number) and (residue.seqidcode == neighbour.seqidcode):
                    idx_this_res = ires
                    index = ires
                    break
            if (idx_this_res < 0):
                self.logger.warning("Failed to find the residue %s %s %s in the container arrays",neighbour.chain_id,neighbour.name,neighbour.number)
                continue
            if (idx_this_res < 1):              #this means that the termini were not properly defined
                self.logger.warning("Terminal-status of residue with a bug.")
                continue
            prev_residue = all_protein_residues[idx_this_res - 1]
            #amino acids
            if AA.IsBaseAminoAcid(neighbour.name) or AA.IsExtendedAminoAcid(neighbour.name) or AA.IsPostTranslationalModification(neighbour.name):
                if (not AA.IsBaseAminoAcid(prev_residue.name)) and (not AA.IsExtendedAminoAcid(prev_residue.name)) and (not AA.IsPostTranslationalModification(prev_residue.name)):
                    continue
                #add C and O atoms from previous residue as caps
                for j,atp in enumerate(prev_residue.atom_types):
                    if atp == "C":
                        #add to main geometry lists
                        neighbour.geometry.append(prev_residue.geometry[j].copy())
                        neighbour.atoms.append("C")
                        neighbour.atom_types.append("C#")
                        neighbour.pdb_lines.append(prev_residue.pdb_lines[j])
                        #add to cap lists
                        neighbour.cap_geometry.append(prev_residue.geometry[j].copy())
                        neighbour.cap_atoms.append("C")
                    elif atp == "O":
                        neighbour.geometry.append(prev_residue.geometry[j].copy())
                        neighbour.atoms.append("O")
                        neighbour.atom_types.append("O#")
                        neighbour.pdb_lines.append(prev_residue.pdb_lines[j])
                        neighbour.cap_geometry.append(prev_residue.geometry[j].copy())
                        neighbour.cap_atoms.append("O")
            elif NA.IsBaseNucleicAcid(neighbour.name) or NA.IsExtendedNucleicAcid(neighbour.name):
                if (not NA.IsBaseNucleicAcid(prev_residue.name)) and (not NA.IsExtendedNucleicAcid(prev_residue.name)):
                    continue
                #add C and O atoms from previous residue as caps
                for j,atp in enumerate(prev_residue.atom_types):
                    if atp == "C3'":
                        #add to main geometry lists
                        neighbour.geometry.append(prev_residue.geometry[j].copy())
                        neighbour.atoms.append("C")
                        neighbour.atom_types.append("C#")
                        neighbour.pdb_lines.append(prev_residue.pdb_lines[j])
                        #add to cap lists
                        neighbour.cap_geometry.append(prev_residue.geometry[j].copy())
                        neighbour.cap_atoms.append("C")
                    elif atp == "O3'":
                        neighbour.geometry.append(prev_residue.geometry[j].copy())
                        neighbour.atoms.append("O")
                        neighbour.atom_types.append("O#")
                        neighbour.pdb_lines.append(prev_residue.pdb_lines[j])
                        neighbour.cap_geometry.append(prev_residue.geometry[j].copy())
                        neighbour.cap_atoms.append("O")
        self.logger.info(f"Added caps to {sum(1 for r in neighbours if r.cap_atoms)} neighbour residues")

    def calculate_deformation_energies(self, ligand_states: List[LigandState]) -> None:
        """Calculate deformation energies for all ligand states."""
        self.logger.info("Calculating ligand deformation energies")
        #acid_base_correction
        max_ligand_charge = -10**1000
        calculations = []
        ip_calculation_mode = CalculationType.IN_POCKET
        opt_calculation_mode = CalculationType.FULL_OPT
        if self.config.skip_deformation_energy:
            ip_calculation_mode = CalculationType.BASE_ENERGY
            opt_calculation_mode = CalculationType.BASE_ENERGY
        #prepare all calculations
        for i,state in enumerate(ligand_states):
            #add index
            ligand_states[i].ligand_index = i
            base_file = str(state.file_path)
            #base energy calculation
            base_calc = (base_file,
                         state.charge,
                         CalculationType.BASE_ENERGY,
                         base_file.replace(".xyz","_base").replace(".sdf","_base"))
            calculations.append(("base",i,base_calc))
            #in-pocket calculation
            ip_file = base_file.replace(".xyz","_ip").replace(".sdf","_ip")
            ip_calc = (base_file,state.charge,ip_calculation_mode,ip_file)
            calculations.append(("ip",i,ip_calc))
            #full optimization (from IP structure)
            opt_file = base_file.replace(".xyz","_opt").replace(".sdf","_opt")
            opt_calc = (f"{ip_file}.xyz",
                        state.charge,
                        opt_calculation_mode,
                        opt_file)
            calculations.append(("opt",i,opt_calc))
            #update information for acid-base corrections
            if state.charge > max_ligand_charge: max_ligand_charge = state.charge
        for state in ligand_states:
            state.acid_base_correction = PKA.AcidBaseCorrection(max_ligand_charge - state.charge)
        #run calculations
        if self.config.parallel_processing: self._calculate_deformation_parallel(ligand_states,calculations)
        else: self._calculate_deformation_sequential(ligand_states,calculations)
        #STRICT sanity check: ligand IP energy must exist.
        #missing (None) or ~0.0 Hartree is treated as a hard failure because it indicates a broken QM run or inconsistent input/charge.
        bad_states = []
        for st in ligand_states:
            eip = st.energy_ip
            ok = True
            if eip is None: ok = False
            else:
                try:
                    v = float(eip)
                    # Total QM energies should never be ~0 Hartree; treat this as a sentinel.
                    if (not math.isfinite(v)) or abs(v) < 1e-8: ok = False
                except (TypeError, ValueError):
                    ok = False
            if not ok: bad_states.append(st)
        if bad_states:
            for st in bad_states:
                idx = int(getattr(st,"ligand_index",-1))
                #idx is 0-based internally; outputs are 1-based
                out_n = idx + 1 if idx >= 0 else "?"
                self.logger.error("FATAL: Missing/zero ligand IP energy for ligand state %s "
                                  "(charge = %s, file = %s). This indicates a failed QM job or invalid input. "
                                  "Please inspect %s/%s/ligand_ip_%s.out and related QM outputs.",
                                  out_n,
                                  getattr(st,"charge","?"),
                                  getattr(st,"file_path","?"),
                                  self.config.output_folder,
                                  self.config.output_folder_qm,
                                  out_n)
            raise InPocketAnalyzerError("Missing/zero ligand IP energy for one or more ligand states. Aborting analysis.")
        #organize files
        self._organize_ligand_files(ligand_states)
        #report results
        self._report_deformation_results(ligand_states)

    def _calculate_deformation_sequential(self,ligand_states: List[LigandState],calculations: List[Tuple[str, int, Tuple]]) -> None:
        """Calculate deformation energies sequentially and write ligand_ip/opt .out."""
        #group calculations by state and type
        calc_by_state: Dict[int,Dict[str,Tuple]] = {}
        for calc_type,state_idx, calc_params in calculations:
            if state_idx not in calc_by_state:
                calc_by_state[state_idx] = {}
            calc_by_state[state_idx][calc_type] = calc_params
        #process each state completely (base -> ip -> opt)
        for state_idx in sorted(calc_by_state.keys()):
            state = ligand_states[state_idx]
            state_calcs = calc_by_state[state_idx]
            #1) base energy calculation (distortion reference)
            if "base" in state_calcs:
                (energy,rmsd,totalpol,atmlist,charges,polarisabilities,converged) = self.qm_calc.calculate_energy(*state_calcs["base"])
                if energy is not None:
                    state.energies["base"] = energy
                    state.rmsd["base"] = rmsd
                    state.converged["base"] = converged
            #2) IP calculation
            if "ip" in state_calcs:
                (energy,rmsd,totalpol,atmlist,charges,polarisabilities,converged) = self.qm_calc.calculate_energy(*state_calcs["ip"])
                if energy is not None:
                    state.energies["ip"] = energy
                    state.rmsd["ip"] = rmsd
                    state.converged["ip"] = converged
                    #write ligand_ip_N.out as in the serial code
                    ip_out = (self.config.output_folder / self.config.output_folder_qm / f"ligand_ip_{state_idx + 1}.out")
                    QM.Write2File(energy,rmsd,totalpol,atmlist,charges,polarisabilities,converged,ip_out)
                    if converged:
                        base_file = str(state.file_path)
                        ip_file = state_calcs["ip"][3]
                        GEOM.ReCenterAndWrite(base_file, ip_file + ".xyz", ip_file)
                        if self.config.qm_params.calculate_density:
                            density_file = self._density_prefix(ip_file)
                            QM.GetDensity(ip_file + ".xyz",str(state.charge),str(self.config.qm_params.electronic_temperature),str(int(self.config.qm_params.use_solvation)),self.config.qm_params.solvent_name,density_file)
            #3) OPT calculation (uses IP structure as input)
            if "opt" in state_calcs:
                (energy,rmsd,totalpol,atmlist,charges,polarisabilities,converged) = self.qm_calc.calculate_energy(*state_calcs["opt"])
                if energy is not None:
                    state.energies["opt"] = energy
                    state.rmsd["opt"] = rmsd
                    state.converged["opt"] = converged
                    #write ligand_opt_N.out as in the serial code
                    opt_out = (self.config.output_folder / self.config.output_folder_qm / f"ligand_opt_{state_idx + 1}.out")
                    QM.Write2File(energy,rmsd,totalpol,atmlist,charges,polarisabilities,converged,opt_out)

    def _calculate_deformation_parallel(self, ligand_states: List[LigandState], calculations: List[Tuple[str, int, Tuple]]) -> None:
        """Calculate deformation energies in parallel and write ligand_ip/opt .out."""
        #stage 1: Base and IP calculations (can run in parallel)
        stage1_calcs: List[Tuple[str,int,Tuple]] = []
        #stage 2: OPT calculations
        stage2_calcs: List[Tuple[str,int,Tuple]] = []
        for calc_type,state_idx,calc_params in calculations:
            if calc_type in ["base","ip"]:
                stage1_calcs.append((calc_type,state_idx,calc_params))
            else:
                stage2_calcs.append((calc_type,state_idx,calc_params))
        #stage 1 (base + ip)
        if stage1_calcs:
            with ProcessPoolExecutor(max_workers = self.config.max_workers) as executor:
                future_to_calc = {executor.submit(self.qm_calc.calculate_energy,*calc[2]): calc
                                  for calc in stage1_calcs}
                for future in as_completed(future_to_calc):
                    calc_type,state_idx,calc_params = future_to_calc[future]
                    state = ligand_states[state_idx]
                    try:
                        (energy,rmsd,totalpol,atmlist,charges,polarisabilities,converged) = future.result()
                        if energy is not None:
                            state.energies[calc_type] = energy
                            state.rmsd[calc_type] = rmsd
                            state.converged[calc_type] = converged
                            if calc_type == "ip":
                                #write ligand_ip_N.out
                                ip_out = (self.config.output_folder / self.config.output_folder_qm / f"ligand_ip_{state_idx + 1}.out")
                                QM.Write2File(energy,rmsd,totalpol,atmlist,charges,polarisabilities,converged,ip_out)
                                #track minimum IP energy for acid-base correction
                                if energy < getattr(self,"minimum_energy",float("inf")):
                                    self.minimum_energy = energy
                    except Exception as e:
                        self.logger.error(f"Stage 1 calculation failed: {e}")
        #after stage 1, re-center IP structures
        for calc_type,state_idx,calc_params in stage1_calcs:
            if calc_type == "ip":
                state = ligand_states[state_idx]
                if state.converged.get("ip",False):
                    base_file = str(state.file_path)
                    ip_file = calc_params[3]
                    try:
                        GEOM.ReCenterAndWrite(base_file,ip_file + ".xyz",ip_file)
                        #calculate the density?
                        if self.config.qm_params.calculate_density:
                            density_file = self._density_prefix(ip_file)
                            QM.GetDensity(ip_file + ".xyz",str(state.charge),str(self.config.qm_params.electronic_temperature),str(int(self.config.qm_params.use_solvation)),self.config.qm_params.solvent_name,density_file)
                    except Exception as e:
                        self.logger.warning("Failed to recenter/write or compute density for %s: %s",ip_file,e)
        #stage 2: OPT calculations
        if stage2_calcs:
            with ProcessPoolExecutor(max_workers = self.config.max_workers) as executor:
                future_to_calc = {executor.submit(self.qm_calc.calculate_energy,*calc[2]): calc
                                  for calc in stage2_calcs}
                for future in as_completed(future_to_calc):
                    calc_type,state_idx,calc_params = future_to_calc[future]
                    state = ligand_states[state_idx]
                    try:
                        (energy,rmsd,totalpol,atmlist,charges,polarisabilities,converged) = future.result()
                        if energy is not None:
                            state.energies[calc_type] = energy
                            state.rmsd[calc_type] = rmsd
                            state.converged[calc_type] = converged
                            if calc_type == "opt":
                                #write ligand_opt_N.out
                                opt_out = (self.config.output_folder / self.config.output_folder_qm / f"ligand_opt_{state_idx + 1}.out")
                                QM.Write2File(energy,rmsd,totalpol,atmlist,charges,polarisabilities,converged,opt_out)
                    except Exception as e:
                        self.logger.error(f"Stage 2 calculation failed: {e}")

    def _organize_ligand_files(self,ligand_states: List[LigandState]) -> None:
        """Organize ligand files into appropriate directories."""
        for i,state in enumerate(ligand_states):
            base_file = str(state.file_path)
            original_name = Path(base_file).stem  #get the original name without extension
            #move IP files
            ip_file = base_file.replace(".xyz","_ip.xyz").replace(".sdf","_ip.xyz")
            if Path(ip_file).exists():
                #use original name pattern instead of generic ligand_ip_{i+1}
                ip_filename = f"{original_name}_ip.xyz"
                self.file_manager.move_file(Path(ip_file),"ip",ip_filename)
            #move OPT files
            opt_file = base_file.replace(".xyz","_opt.xyz").replace(".sdf","_opt.xyz")
            if Path(opt_file).exists():
                #use original name pattern instead of generic ligand_opt_{i+1}
                opt_filename = f"{original_name}_opt.xyz"
                self.file_manager.move_file(Path(opt_file),"opt",opt_filename)
            #move original ligand files produced by protonation
            #IMPORTANT: geometries/ligand_states/ is a stable cache.
            #we keep those files in place, but we still move derived *_ip.xyz / *_opt.xyz outputs.
            base_path = Path(base_file)
            if base_path.parent.name == "ligand_states":
                #cached ligand state: keep it where it is
                continue
            if base_file.endswith(".xyz"):
                self.file_manager.move_file(base_path,"start_xyz")
            elif base_file.endswith(".sdf"):
                self.file_manager.move_file(base_path,"start_sdf")
        #move base ligand PDB file (active residue entry)
        ligand_pdb = self.file_manager.get_path("geometries",f"{self.config.active_residue_chain}{self.config.active_residue}.pdb")
        if ligand_pdb.exists():
            self.file_manager.move_file(ligand_pdb,"start_pdb")

    def _final_geometry_cleanup(self) -> None:
        """Geometry cleanup
           Rules (applied to files in the *top-level* geometries folder):
           - *_ip.xyz     -> geometries/ip/
           - *_opt.xyz    -> geometries/opt/
           - other .xyz   -> geometries/start_xyz/
           - .sdf         -> geometries/start_sdf/
           - .smi         -> geometries/smiles/
           - .pdb         -> geometries/start_pdb/
           - .png         -> deleted"""
        geom_root = self.file_manager.dirs.get("geometries")
        if geom_root is None:
            return
        try:
            entries = list(geom_root.iterdir())
        except Exception as e:
            self.logger.warning("Could not list geometry folder %s during final cleanup: %s",geom_root,e)
            return
        for item in entries:
            #only handle top-level files; skip subdirectories (ip/, opt/, etc.)
            if not item.is_file(): continue
            name = item.name
            suffix = item.suffix.lower()
            try:
                if suffix == ".xyz":
                    #decide the category based on the filename
                    if name.endswith("_ip.xyz"): dest_cat = "ip"
                    elif name.endswith("_opt.xyz"): dest_cat = "opt"
                    else: dest_cat = "start_xyz"
                    self.file_manager.move_file(item,dest_cat)
                elif suffix == ".sdf":
                    self.file_manager.move_file(item,"start_sdf")
                elif suffix == ".smi":
                    self.file_manager.move_file(item,"smiles")
                elif suffix == ".pdb":
                    self.file_manager.move_file(item,"start_pdb")
                elif suffix == ".png":
                    #keep or delete PNG previews (residue_prots.png, ligand_prots.png, etc.)
                    if getattr(self.config,"keep_png",False):
                        #place png/ next to geometries/ and output/ (i.e., inside output_folder root)
                        png_dir = self.config.output_folder / "png"
                        png_dir.mkdir(parents = True,exist_ok = True)
                        dest = png_dir / item.name
                        #avoid overwriting: add incremental suffix if needed
                        if dest.exists():
                            stem = dest.stem
                            suf = dest.suffix
                            n = 1
                            while True:
                                cand = png_dir / f"{stem}_{n}{suf}"
                                if not cand.exists():
                                    dest = cand
                                    break
                                n += 1
                        shutil.move(str(item),str(dest))
                    else:
                        item.unlink()
            except Exception as e:
                self.logger.warning("Failed to clean geometry file %s: %s",item,e)

    def _write_residue_csv_for_ligand(self, ligand_name: str, ligand_charge: int, residues: List[RES.ResidueData], target_file: Path) -> None:
        """Write a CSV-like text file describing one ligand + all residue states.
           Format:
               path-to-geometries: <ip_dir>/
               extension: .xyz
               ligand: <ligand_name> <ligand_charge>
               r:RES;CHAIN;RESNUM;STATE;RES_FILE_BASENAME
           This is intended to be consumed by downstream IPA workflows."""
        ip_dir = self.file_manager.dirs.get("ip")
        if ip_dir is None:
            self.logger.warning("IP directory not found; cannot write residue CSV '%s'",target_file)
            return
        #IP geometries are written as .xyz
        extension = ".xyz"
        try:
            target_file.parent.mkdir(parents = True,exist_ok = True)
            with target_file.open("w") as f:
                #header lines
                f.write(f"path-to-geometries: {ip_dir}/\n")
                f.write(f"extension: {extension}\n")
                f.write(f"ligand: {ligand_name} {ligand_charge}\n")
                #one line per residue/protonation state
                for res in residues:
                    name = getattr(res,"name","")
                    residue_id = getattr(res,"residue_id","")
                    #NOTE: we no longer skip flipped residues here.
                    #flip information is encoded in name/residue_id/ip_basename
                    #(e.g. HIS_flip_A179_prot_3_ip).
                    resname = (name or "")[:3]
                    chain = getattr(res,"chain_id","")
                    number = getattr(res,"number",0)
                    state = getattr(res,"protonation_state","0")
                    prot_file = getattr(res,"protonated_file",None)
                    #derive the IP basename:
                    #  if we know the protonated geometry filename, use its stem:
                    #      LYS_A78.xyz            -> LYS_A78_ip
                    #      LYS_A78_prot_2.xyz     -> LYS_A78_prot_2_ip
                    #      HIS_flip_A179_prot_3   -> HIS_flip_A179_prot_3_ip
                    #  otherwise, fall back to residue_id-based name.
                    if prot_file:
                        stem = Path(prot_file).stem
                        if stem.endswith("_ip"): ip_basename = stem
                        else: ip_basename = f"{stem}_ip"
                    else:
                        if residue_id: ip_basename = f"{residue_id}_ip"
                        else: ip_basename = f"{resname}_{chain}{number}_ip"
                    f.write(f"r:{resname};{chain};{number};{state};{ip_basename}\n")
        except Exception as e:
            self.logger.error("Failed to write residue CSV '%s': %s",target_file,e)

    def _write_best_unique_residue_csv( self, ligand_name: str, ligand_charge: int, residues: List[RES.ResidueData], best_results: List[BindingResult], target_file: Path) -> None:
        """Write best_list_residues_unique.csv: one entry per physical residue (RES, chain, resnum),
           choosing the *most favorable* variant among protonation/flip states for the best ligand.
           Selection criterion: minimal BindingResult.binding_energy (most negative, kcal/mol)."""
        #map residue_id -> ResidueData for reliable ip_basename generation
        res_by_id: Dict[str,RES.ResidueData] = {r.residue_id: r for r in residues}
        #pick best variant per (resname,chain,number)
        best_variant_by_key: Dict[Tuple[str,str,int],BindingResult] = {}
        for br in best_results:
            rid = getattr(br,"residue_name","")
            if not rid: continue
            r = res_by_id.get(rid)
            if r is None:
                #if we cannot resolve the residue object, skip (safer than guessing filenames)
                continue
            resname = (getattr(r,"name","") or "")[:3]
            chain = getattr(r,"chain_id","")
            number = int(getattr(r,"number",0))
            key = (resname,chain,number)
            prev = best_variant_by_key.get(key)
            if (prev is None) or (br.binding_energy < prev.binding_energy):
                best_variant_by_key[key] = br
        #convert selected BindingResults -> selected ResidueData list (stable order)
        selected_residues: List[RES.ResidueData] = []
        for key in sorted(best_variant_by_key.keys(),key = lambda x: (x[1],x[2],x[0])):
            br = best_variant_by_key[key]
            r = res_by_id.get(br.residue_name)
            if r is not None: selected_residues.append(r)
        #reuse the standard writer to guarantee consistent formatting + filename logic
        self._write_residue_csv_for_ligand(ligand_name = ligand_name,ligand_charge = ligand_charge,residues = selected_residues,target_file = target_file)

    def _collect_ip_basenames_from_residue_csv(self, csv_path: Path) -> List[str]:
        """Parse a residue CSV (best_list_residues*.csv) and return the list of IP basenames
           that must exist in geometries/ip as <basename>.xyz.
           We collect:
             - ligand basename from: 'ligand: <name> <charge>'
             - residue basenames from: 'r:...;...;...;...;<basename>'"""
        basenames: List[str] = []
        if not csv_path.exists(): return basenames
        try:
            for line in csv_path.read_text().splitlines():
                line = line.strip()
                if not line: continue
                if line.startswith("ligand:"):
                    #ligand: Z9999_prot_1_tautomer_1_ip -1
                    parts = line.split()
                    if len(parts) >= 2: basenames.append(parts[1])
                elif line.startswith("r:"):
                    #r:MET;A;69;0;MET_A69_ip
                    #basename is after last ';'
                    base = line.split(";")[-1].strip()
                    if base:
                        basenames.append(base)
        except Exception as e:
            self.logger.warning("Failed to parse residue CSV %s: %s",csv_path,e)
        #unique, stable order
        seen = set()
        out: List[str] = []
        for b in basenames:
            if b not in seen:
                seen.add(b)
                out.append(b)
        return out

    def _protect_and_restore_best_ip_xyz(self) -> None:
        """In compact --best mode, cleanup_outputs() removes geometries/ip,
           but the fragment analysis tools need the IP xyz files referenced in best_list_residues*.csv.
           This function:
             - reads best_list_residues_unique.csv if present, else best_list_residues.csv
             - copies required geometries/ip/<basename>.xyz to a temp folder
             - after cleanup, restores those xyz files back into geometries/ip/"""
        #only relevant for compact best
        if (not self.config.best_mode) or self.config.best_debug: return
        geom_root = self.config.output_folder / self.config.geometry_folder
        ip_dir = geom_root / "ip"
        best_csv = self.config.output_folder / "best_list_residues.csv"
        best_unique_csv = self.config.output_folder / "best_list_residues_unique.csv"
        chosen_csv = best_unique_csv if best_unique_csv.exists() else best_csv
        needed_basenames = self._collect_ip_basenames_from_residue_csv(chosen_csv)
        if not needed_basenames: return
        tmp_keep = self.config.output_folder / "_ip_keep_tmp"
        try:
            if tmp_keep.exists():
                shutil.rmtree(tmp_keep,ignore_errors = True)
            tmp_keep.mkdir(parents = True,exist_ok = True)
            #copy out required xyz files (if they exist)
            for base in needed_basenames:
                src = ip_dir / f"{base}.xyz"
                if src.exists():shutil.copy2(src,tmp_keep / src.name)
                else: self.logger.warning("Requested IP xyz missing before cleanup: %s",src)
        except Exception as e:
            self.logger.warning("Failed to protect IP xyz before cleanup: %s",e)
            return
        #cleanup happens outside this function
        #we restore in a second step after cleanup_outputs() is called
        #to do that cleanly, we store tmp_keep path on self and restore later
        self._ip_keep_tmp = tmp_keep

    def _compute_total_interaction_energy_and_list_res(self, results: List[BindingResult]) -> Tuple[float, List[str]]:
        """Reproduce the same residue-selection logic used in write_results()
           to compute the "IPA Interaction Energy" and the residue list contributing to it."""
        #when residues are NOT capped (--no-cap), multiple protonation variants can share the same *net charge* because a backbone N proton can compensate a side-chain deprotonation (or vice-versa)
        #in that regime we must select exactly ONE "standard" variant per residue
        if not getattr(self.config,"cap_residues",True):
            return self._compute_total_interaction_energy_no_cap(results)
        best_by_residue = {}
        for result in results:
            #parse residue info properly to get protonation state
            residue_parts = result.residue_name.replace("_flip","").split("_")
            base_name = residue_parts[0]
            chain_res = residue_parts[1] if len(residue_parts) > 1 else "A1"
            #extract chain and residue number
            chain_id = "".join([c for c in chain_res if not c.isdigit()]).replace("f","")
            res_num = "".join([c for c in chain_res if c.isdigit()])
            state_raw = result.residue_state
            charge_int = FORMAT.Parse_Charge_State(state_raw)
            proton_state = FORMAT.Format_Charge_State(charge_int)
            #flipped?
            is_flipped = result.residue_name.endswith("f") or ("_flip" in result.residue_name)
            if is_flipped and base_name not in ("THR","SER","TYR"):
                continue
            take = False
            if base_name in ("ARG","LYS"):
                take = (charge_int == 1)
            elif base_name in ("TYR","HIS"):
                take = (charge_int == 0)
            elif base_name in ("ASP","GLU"):
                take = (charge_int == -1)
            else:
                take = True
            if not take:
                continue
            key = f"{base_name}_{chain_id}{res_num}"
            if (key not in best_by_residue or result.binding_energy < best_by_residue[key]["energy"]):
                best_by_residue[key] = {"energy": result.binding_energy,
                                        "picked": f"{base_name}_{chain_id}{res_num}_{proton_state}"}
        total_interaction = sum(v["energy"] for v in best_by_residue.values())
        picked = [v["picked"] for v in best_by_residue.values()]
        return total_interaction,picked

    def _standard_charge_target(self, base_name: str) -> int:
        """Target charge for the single "standard" variant used in totals when --no-cap."""
        b = str(base_name or "").strip().upper()
        if b in ("ARG", "LYS"): return 1
        if b in ("ASP", "GLU"): return -1
        #TYR and HIS are neutral in our "standard" set
        return 0

    def _processed_residue_variant(self, residue_id: str):
        """Return processed ResidueData for a residue_id (best-effort)."""
        m = getattr(self,"_processed_residues_by_id",None) or {}
        return m.get(residue_id)

    def _no_cap_atom_neighbor_count(self, residue, atom_type_candidates: Tuple[str, ...], elem: str, max_map_dist: float = 1.2, neigh_scale: float = 1.3) -> Optional[int]:
        """Map an atom from the *original* residue by atom_type name(s) to the processed geometry
           (via nearest-by-coordinate matching of element), then count covalent neighbours."""
        #original
        atypes = list(getattr(residue,"atom_types",[]) or [])
        geom0 = list(getattr(residue,"geometry",[]) or [])
        if not atypes or not geom0: 
            return None
        wanted = {s for s in (str(x).strip().upper() for x in (atom_type_candidates or ())) if s}
        orig_idx = None
        for i,nm in enumerate(atypes):
            if str(nm).strip().upper() in wanted:
                orig_idx = i
                break
        if orig_idx is None: 
            return None
        ref = geom0[orig_idx]
        if not (isinstance(ref,(list,tuple)) and len(ref) == 3):
            return None
        #processed
        pe = list(getattr(residue,"processed_atoms",[]) or [])
        pg = list(getattr(residue,"processed_geometry",[]) or [])
        if not pe or not pg or len(pe) != len(pg):
            return None
        proc_idx = CONNECT._nearest_index_by_element_and_coord(pe,pg,elem,ref,max_dist = max_map_dist)
        if proc_idx is None:
            return None
        return CONNECT._count_covalent_neighbors(pe,pg,proc_idx,scale = neigh_scale)

    def _no_cap_variant_penalty(self, base_name: str, residue_variant) -> int:
        """Penalty score: lower is better. Non-checkable variants are penalized but not discarded."""
        b = str(base_name or "").strip().upper()
        penalty = 0
        #backbone-N check requested for HIS/LYS/TYR in --no-cap mode
        if b in ("HIS","LYS","TYR"):
            n_nb = self._no_cap_atom_neighbor_count(residue_variant,("N",),"N")
            if n_nb is None:
                penalty += 50
            elif n_nb != 3:
                penalty += 100
        #HIS tautomer check: want HID/HSD (proton on ND)
        if b == "HIS":
            nd_nb = self._no_cap_atom_neighbor_count(residue_variant,("ND1","ND"),"N")
            if nd_nb is None:
                penalty += 50
            elif nd_nb != 3:
                penalty += 100
        return penalty

    def _compute_total_interaction_energy_no_cap(self, results: List[BindingResult]) -> Tuple[float, List[str]]:
        """Compute total interaction by selecting exactly one "standard" variant per residue.
           Rules:
             - flipped variants are excluded entirely from the total
             - pick a single variant per (resname, chain, resnum)
             - default standard charge is 0; ARG/LYS: +1; ASP/GLU: -1
             - for HIS/LYS/TYR disambiguate further using neighbour-count checks on the backbone N
               (and ND1 for HIS)"""
        total_interaction = 0.0
        picked: List[str] = []
        #group candidates by residue identity (ignore protonation index)
        grouped: Dict[Tuple[str,str,str],List[BindingResult]] = {}
        for result in results or []:
            #flipped?
            is_flipped = result.residue_name.endswith("f") or ("_flip" in result.residue_name)
            if is_flipped:
                continue
            residue_parts = result.residue_name.replace("_flip","").split("_")
            base_name = residue_parts[0] if residue_parts else ""
            chain_res = residue_parts[1] if len(residue_parts) > 1 else "A1"
            chain_id = "".join([c for c in chain_res if not c.isdigit()]).replace("f","")
            res_num = "".join([c for c in chain_res if c.isdigit()])
            key = (str(base_name).strip().upper(),chain_id,res_num)
            grouped.setdefault(key,[]).append(result)
        #select exactly one candidate per residue
        for (base_u,chain_id,res_num),cand_list in grouped.items():
            target_q = self._standard_charge_target(base_u)
            best: Optional[BindingResult] = None
            best_pen = float("inf")
            best_e = float("inf")
            best_q = 0
            #1) primary: exact target charge
            for r in cand_list:
                q = FORMAT.Parse_Charge_State(getattr(r,"residue_state",None))
                if q != target_q:
                    continue
                rv = self._processed_residue_variant(getattr(r,"residue_name",""))
                pen = self._no_cap_variant_penalty(base_u,rv) if rv is not None else 50
                e = float(getattr(r,"binding_energy",0.0) or 0.0)
                if (pen < best_pen) or (pen == best_pen and e < best_e):
                    best = r
                    best_pen = pen
                    best_e = e
                    best_q = q
            #2) fallback: closest charge (keeps totals robust if a target charge is missing)
            if best is None:
                for r in cand_list:
                    q = FORMAT.Parse_Charge_State(getattr(r,"residue_state",None))
                    delta = abs(int(q) - int(target_q))
                    rv = self._processed_residue_variant(getattr(r,"residue_name",""))
                    pen = (1000 + 10*delta) + (self._no_cap_variant_penalty(base_u,rv) if rv is not None else 50)
                    e = float(getattr(r,"binding_energy",0.0) or 0.0)
                    if (pen < best_pen) or (pen == best_pen and e < best_e):
                        best = r
                        best_pen = pen
                        best_e = e
                        best_q = q
            if best is None:
                continue
            total_interaction += float(getattr(best,"binding_energy",0.0) or 0.0)
            picked.append(f"{base_u}_{chain_id}{res_num}_{FORMAT.Format_Charge_State(best_q)}")
        return total_interaction, picked

    def _swap_sdf_coordinates_and_write(self, template_sdf: Path, xyz_file: Path, out_sdf: Path, props: Dict[str, Union[str, float, int]]) -> None:
        """Create an SDF by taking bonding/connectivity from template_sdf and
           replacing the molblock coordinates with coordinates from xyz_file.
           Then append SDF property blocks (props) before $$$$."""
        SANITY.InputSanitizer.write_template_sdf_with_xyz(template_sdf = template_sdf,xyz_file = xyz_file,out_sdf = out_sdf,props = props,logger = self.logger)

    def _write_ligand_states_ip_sdfs(self, ligand_states: List[LigandState], binding_results: List[List[BindingResult]]) -> None:
        """Create geometries/ligand_states_ip/*.sdf:
           - connectivity from geometries/ligand_states/<state>.sdf (cache) if present
           - coordinates from geometries/ip/<state_stem>_ip.xyz
           - append SMILES / IPA_EBIND / LIST_RES props"""
        #only for normal + --best-debug, not for --best
        create_all = (not self.config.best_mode) or self.config.best_debug
        if not create_all:
            return
        geom_root = self.config.output_folder / self.config.geometry_folder
        cache_dir = geom_root / "ligand_states"
        ip_dir = geom_root / "ip"
        out_dir = geom_root / "ligand_states_ip"
        out_dir.mkdir(parents = True,exist_ok = True)
        #need ligand_min_energy to reproduce IPA_EBIND formula used in stats.txt
        ligand_min_energy = float("inf")
        for st in ligand_states:
            if st.energy_ip is not None:
                ligand_min_energy = min(ligand_min_energy,st.energy_ip)
        for lig_idx,state in enumerate(ligand_states):
            stem = Path(state.file_path).stem
            #template sdf: prefer cached ligand_states/<name>.sdf
            template = cache_dir / f"{stem}.sdf"
            if not template.exists():
                #fallback: only if original is an SDF
                if Path(state.file_path).suffix.lower() == ".sdf":
                    template = Path(state.file_path)
                else:
                    native_candidates: List[Path] = []
                    start_sdf_dir = geom_root / "start_sdf"
                    #original/native amino-acid-like ligands often end up as <stem>H.xyz
                    #while the preserved template remains <stem>.sdf.
                    if stem.endswith("H"):
                        native_candidates.extend([cache_dir / f"{stem[:-1]}.sdf",geom_root / f"{stem[:-1]}.sdf",start_sdf_dir / f"{stem[:-1]}.sdf"])
                    native_candidates.extend([cache_dir / f"{stem}.sdf",geom_root / f"{stem}.sdf",start_sdf_dir / f"{stem}.sdf"])
                    template = next((cand for cand in native_candidates if cand.exists()), None)
                    if template is None:
                        self.logger.debug("Skipping ligand_states_ip for %s: no SDF template available",stem)
                        continue
                    else:
                        self.logger.info("ligand_states_ip for %s: using fallback SDF template %s.",stem,template.name)
            xyz_ip = ip_dir / f"{stem}_ip.xyz"
            if not xyz_ip.exists():
                #if IP did not converge / no file produced, skip
                self.logger.debug("Skipping ligand_states_ip for %s: missing IP xyz (%s)",stem,xyz_ip)
                continue
            #interaction energy + residue list (same logic as write_results)
            total_interaction,picked_res = self._compute_total_interaction_energy_and_list_res(binding_results[lig_idx])
            #IPA_EBIND (same definition as stats.txt)
            ipa_ebind = None
            if state.deformation_energy is not None and state.energy_ip is not None:
                ipa_ebind = (state.deformation_energy + (state.energy_ip - ligand_min_energy + state.acid_base_correction)*self.Eh2kcalmol + total_interaction)
            #match best_lig.sdf conventions:
            #  IPA_EBIND printed with 3 decimals
            #  LIST_RES as chain+resnum list (e.g. "A69,A124,...")
            list_res_simple: List[str] = []
            for r in picked_res:
                #picked_res entries look like "MET_A69_0" etc.
                parts = r.split("_")
                if len(parts) >= 2:
                    chain_res = parts[1].replace("f", "")  #e.g. "A69"
                    list_res_simple.append(chain_res)
            #match best_lig.sdf ordering: lexicographic sort (string sort)
            #example: "A100" comes before "A69"
            list_res_sorted = sorted(list_res_simple)
            props = {"SMILES": (state.smiles or ""),"IPA_EBIND": (f"{ipa_ebind:.3f}" if ipa_ebind is not None else "NA"),"LIST_RES": ",".join(list_res_sorted)}
            out_sdf = out_dir / f"{stem}_ip.sdf"
            try:
                self._swap_sdf_coordinates_and_write(template,xyz_ip,out_sdf,props)
            except Exception as e:
                self.logger.warning("Failed to write ligand_states_ip SDF for %s: %s",stem,e)

    def _report_deformation_results(self, ligand_states: List[LigandState]) -> None:
        """Report deformation energy results."""
        for i,state in enumerate(ligand_states):
            self.logger.info(f"Ligand state {i + 1} results:")
            if state.deformation_energy is not None:
                self.logger.info(f"  Deformation Energy: {state.deformation_energy:.2f} kcal/mol")
            if state.distortion_energy is not None:
                self.logger.info(f"  Distortion Energy: {state.distortion_energy:.2f} kcal/mol")
            #note: no relaxation_energy here; LigandState doesn't define it.
            #any "relaxation" contributions are folded into the IPA Binding Energy and related quantities reported later in write_results().
            
    def _validate_prepared_ligand_states(self, ligand_states: List[LigandState]) -> None:
        """Validate ligand structures after protonation/tautomerisation and before QM."""
        for state in (ligand_states or []):
            path = Path(state.file_path)
            fix_sdf_h = (self.config.protonation_mode != ProtonationMode.ORIGINAL)
            try:
                if path.suffix.lower() == ".sdf":
                    SANITY.InputSanitizer.repair_prepared_sdf_h_connectivity(path,logger = self.logger)
                SANITY.InputSanitizer.validate_prepared_structure(path,logger = self.logger,fix_sdf_h = fix_sdf_h)
            except Exception as e:
                self.logger.error("Prepared ligand structure validation failed for %s: %s",path.name,e)
                raise InPocketAnalyzerError(f"Prepared ligand structure validation failed for {path.name}: {e}")

    def _validate_prepared_residues(self, residues: List[RES.ResidueData]) -> None:
        """Validate processed residue structures before QM pair calculations."""
        for residue in (residues or []):
            prot_file = getattr(residue,"protonated_file",None)
            if not prot_file:
                continue
            path = Path(str(prot_file))
            try:
                if path.suffix.lower() == ".sdf":
                    SANITY.InputSanitizer.repair_prepared_sdf_h_connectivity(path,logger = self.logger)
                SANITY.InputSanitizer.validate_prepared_structure(path,logger = self.logger,fix_sdf_h = (self.config.protonation_mode != ProtonationMode.ORIGINAL))
            except Exception as e:
                self.logger.error("Prepared residue structure validation failed for %s (%s): %s",getattr(residue,"residue_id",path.name),path.name,e)
                raise InPocketAnalyzerError(f"Prepared residue structure validation failed for {getattr(residue,'residue_id',path.name)}: {e}")

    def process_neighbor_residues(self, neighbor_residues: List[RES.ResidueData]) -> List[RES.ResidueData]:
        """Process neighbor residues with protonation."""
        processed = self.residue_processor.process_residues(neighbor_residues)
        self._validate_prepared_residues(processed)
        return processed

    def calculate_binding_energies(self, ligand_states: List[LigandState], processed_residues: List[RES.ResidueData]) -> List[List[BindingResult]]:
        """Calculate binding energies between ligands and residues."""
        return self.binding_calculator.calculate_binding_energies(ligand_states,processed_residues)

    def _flip_residue(self, residue: RES.ResidueData) -> RES.ResidueData:
        """Create flipped version of residue."""
        flipped = copy.deepcopy(residue)
        if residue.name == "ASN":
            #swap OD1 with ND2 or OD2 with ND1
            for i,atp in enumerate(residue.atom_types):
                if atp == "OD1":
                    flipped.atom_types[i] = "ND2"
                    flipped.atoms[i] = "N"
                    coordinate_to_find = str("{:.3f}".format(flipped.geometry[i][2]))
                    string_to_replace = flipped.pdb_lines[i][flipped.pdb_lines[i].rfind(coordinate_to_find):]
                    new_string = string_to_replace.replace("O","N")
                    flipped.pdb_lines[i] = (flipped.pdb_lines[i].replace("OD1","ND2").replace(string_to_replace,new_string))
                elif atp == "OD2":
                    flipped.atom_types[i] = "ND1"
                    flipped.atoms[i] = "N"
                    coordinate_to_find = str("{:.3f}".format(flipped.geometry[i][2]))
                    string_to_replace = flipped.pdb_lines[i][flipped.pdb_lines[i].rfind(coordinate_to_find):]
                    new_string = string_to_replace.replace("O","N")
                    flipped.pdb_lines[i] = (flipped.pdb_lines[i].replace("OD2","ND1").replace(string_to_replace,new_string))
                elif atp == "ND1":
                    flipped.atom_types[i] = "OD2"
                    flipped.atoms[i] = "O"
                    coordinate_to_find = str("{:.3f}".format(flipped.geometry[i][2]))
                    string_to_replace = flipped.pdb_lines[i][flipped.pdb_lines[i].rfind(coordinate_to_find):]
                    new_string = string_to_replace.replace("N","O")
                    flipped.pdb_lines[i] = (flipped.pdb_lines[i].replace("ND1","OD2").replace(string_to_replace,new_string))
                elif atp == "ND2":
                    flipped.atom_types[i] = "OD1"
                    flipped.atoms[i] = "O"
                    coordinate_to_find = str("{:.3f}".format(flipped.geometry[i][2]))
                    string_to_replace = flipped.pdb_lines[i][flipped.pdb_lines[i].rfind(coordinate_to_find):]
                    new_string = string_to_replace.replace("N", "O")
                    flipped.pdb_lines[i] = (flipped.pdb_lines[i].replace("ND2","OD1").replace(string_to_replace,new_string))
        elif residue.name == "GLN":
            #swap OE1 with NE2 or OE2 with NE1
            for i,atp in enumerate(residue.atom_types):
                if atp == "OE1":
                    flipped.atom_types[i] = "NE2"
                    flipped.atoms[i] = "N"
                    coordinate_to_find = str("{:.3f}".format(flipped.geometry[i][2]))
                    string_to_replace = flipped.pdb_lines[i][flipped.pdb_lines[i].rfind(coordinate_to_find):]
                    new_string = string_to_replace.replace("O","N")
                    flipped.pdb_lines[i] = (flipped.pdb_lines[i].replace("OE1","NE2").replace(string_to_replace,new_string))
                elif atp == "OE2":
                    flipped.atom_types[i] = "NE1"
                    flipped.atoms[i] = "N"
                    coordinate_to_find = str("{:.3f}".format(flipped.geometry[i][2]))
                    string_to_replace = flipped.pdb_lines[i][flipped.pdb_lines[i].rfind(coordinate_to_find):]
                    new_string = string_to_replace.replace("O","N")
                    flipped.pdb_lines[i] = (flipped.pdb_lines[i].replace("OE2","NE1").replace(string_to_replace,new_string))
                elif atp == "NE1":
                    flipped.atom_types[i] = "OE2"
                    flipped.atoms[i] = "O"
                    coordinate_to_find = str("{:.3f}".format(flipped.geometry[i][2]))
                    string_to_replace = flipped.pdb_lines[i][flipped.pdb_lines[i].rfind(coordinate_to_find):]
                    new_string = string_to_replace.replace("N","O")
                    flipped.pdb_lines[i] = (flipped.pdb_lines[i].replace("NE1","OE2").replace(string_to_replace,new_string))
                elif atp == "NE2":
                    flipped.atom_types[i] = "OE1"
                    flipped.atoms[i] = "O"
                    coordinate_to_find = str("{:.3f}".format(flipped.geometry[i][2]))
                    string_to_replace = flipped.pdb_lines[i][
                        flipped.pdb_lines[i].rfind(coordinate_to_find) :
                    ]
                    new_string = string_to_replace.replace("N","O")
                    flipped.pdb_lines[i] = (flipped.pdb_lines[i].replace("NE2","OE1").replace(string_to_replace,new_string))
        elif residue.name == "HIS":
            #complex HIS flipping, 4 atoms swapped
            for i,atp in enumerate(residue.atom_types):
                if atp == "ND1":
                    flipped.atom_types[i] = "CD2"
                    flipped.atoms[i] = "C"
                    coordinate_to_find = str("{:.3f}".format(flipped.geometry[i][2]))
                    string_to_replace = flipped.pdb_lines[i][flipped.pdb_lines[i].rfind(coordinate_to_find):]
                    new_string = string_to_replace.replace("N","C")
                    flipped.pdb_lines[i] = (flipped.pdb_lines[i].replace("ND1","CD2").replace(string_to_replace,new_string))
                elif atp == "ND2":
                    flipped.atom_types[i] = "CD1"
                    flipped.atoms[i] = "C"
                    coordinate_to_find = str("{:.3f}".format(flipped.geometry[i][2]))
                    string_to_replace = flipped.pdb_lines[i][flipped.pdb_lines[i].rfind(coordinate_to_find):]
                    new_string = string_to_replace.replace("N","C")
                    flipped.pdb_lines[i] = (flipped.pdb_lines[i].replace("ND2","CD1").replace(string_to_replace,new_string))
                elif atp == "CD1":
                    flipped.atom_types[i] = "ND2"
                    flipped.atoms[i] = "N"
                    coordinate_to_find = str("{:.3f}".format(flipped.geometry[i][2]))
                    string_to_replace = flipped.pdb_lines[i][flipped.pdb_lines[i].rfind(coordinate_to_find):]
                    new_string = string_to_replace.replace("C","N")
                    flipped.pdb_lines[i] = (flipped.pdb_lines[i].replace("CD1","ND2").replace(string_to_replace,new_string))
                elif atp == "CD2":
                    flipped.atom_types[i] = "ND1"
                    flipped.atoms[i] = "N"
                    coordinate_to_find = str("{:.3f}".format(flipped.geometry[i][2]))
                    string_to_replace = flipped.pdb_lines[i][flipped.pdb_lines[i].rfind(coordinate_to_find):]
                    new_string = string_to_replace.replace("C","N")
                    flipped.pdb_lines[i] = (flipped.pdb_lines[i].replace("CD2","ND1").replace(string_to_replace,new_string))
                elif atp == "NE1":
                    flipped.atom_types[i] = "CE2"
                    flipped.atoms[i] = "C"
                    coordinate_to_find = str("{:.3f}".format(flipped.geometry[i][2]))
                    string_to_replace = flipped.pdb_lines[i][flipped.pdb_lines[i].rfind(coordinate_to_find):]
                    new_string = string_to_replace.replace("N","C")
                    flipped.pdb_lines[i] = (flipped.pdb_lines[i].replace("NE1","CE2").replace(string_to_replace,new_string))
                elif atp == "NE2":
                    flipped.atom_types[i] = "CE1"
                    flipped.atoms[i] = "C"
                    coordinate_to_find = str("{:.3f}".format(flipped.geometry[i][2]))
                    string_to_replace = flipped.pdb_lines[i][flipped.pdb_lines[i].rfind(coordinate_to_find):]
                    new_string = string_to_replace.replace("N","C")
                    flipped.pdb_lines[i] = (flipped.pdb_lines[i].replace("NE2","CE1").replace(string_to_replace,new_string))
                elif atp == "CE1":
                    flipped.atom_types[i] = "NE2"
                    flipped.atoms[i] = "N"
                    coordinate_to_find = str("{:.3f}".format(flipped.geometry[i][2]))
                    string_to_replace = flipped.pdb_lines[i][flipped.pdb_lines[i].rfind(coordinate_to_find):]
                    new_string = string_to_replace.replace("N","C")
                    flipped.pdb_lines[i] = (flipped.pdb_lines[i].replace("CE1","NE2").replace(string_to_replace,new_string))
                elif atp == "CE2":
                    flipped.atom_types[i] = "NE1"
                    flipped.atoms[i] = "N"
                    coordinate_to_find = str("{:.3f}".format(flipped.geometry[i][2]))
                    string_to_replace = flipped.pdb_lines[i][flipped.pdb_lines[i].rfind(coordinate_to_find):]
                    new_string = string_to_replace.replace("N","C")
                    flipped.pdb_lines[i] = (flipped.pdb_lines[i].replace("CE2","NE1").replace(string_to_replace,new_string))
        return flipped

    def write_results(self, ligand_states: List[LigandState], binding_results: List[List[BindingResult]]) -> None:
        """Write analysis results to files."""
        self.logger.info("Writing results")
        output_dir = self.config.output_folder
        ligand_int_energies = [0.0]*len(ligand_states)
        ligand_min_energy = float("inf")
        #write interaction energies
        with open(output_dir / "interaction_energies.csv", "w") as f:
            for lig_idx,(ligand,results) in enumerate(zip(ligand_states,binding_results)):
                if (self.config.protonation_mode == ProtonationMode.RDKIT and ligand.smiles):
                    f.write(f"{ligand.smiles}:\n")
                else:
                    f.write(f"Ligand state {lig_idx + 1}/{len(ligand_states)}: \n")
                #STRICT: energy_ip must be present and valid
                _eip = ligand.energy_ip
                if _eip is None:
                    self.logger.error("FATAL: ligand energy_ip is missing for ligand state %d. "
                                      "QM job likely failed; check input and QM outputs.",
                                      lig_idx + 1)
                    raise InPocketAnalyzerError("Missing ligand IP energy (energy_ip = None).")
                if (not math.isfinite(float(_eip))) or abs(float(_eip)) < 1e-8:
                    self.logger.error("FATAL: ligand energy_ip is invalid/zero for ligand state %d: %s. "
                                      "QM job likely failed; check input and QM outputs.",
                                      lig_idx + 1,
                                      _eip)
                    raise InPocketAnalyzerError("Invalid/zero ligand IP energy.")
                ligand_min_energy = min(ligand_min_energy,float(_eip))
                if self.config.edda:
                    f.write("res;chain;resn;charge;res+resn;flip;clash;Eint(kcal/mol);CapCorr(kcal/mol)\n")
                else:
                    f.write("res;chain;resn;charge;res+resn;flip;Eint(kcal/mol);CapCorr(kcal/mol)\n")
                total_interaction = 0.0
                for res_j,result in enumerate(results):
                    #parse residue info properly to get protonation state
                    residue_parts = result.residue_name.replace("_flip","").split("_")
                    base_name = residue_parts[0]  #e.g., "TYR"
                    chain_res = (residue_parts[1] if len(residue_parts) > 1 else "A1")  #e.g., "B253"
                    #extract chain and residue number
                    chain_id = "".join([c for c in chain_res if not c.isdigit()]).replace("f","")
                    res_num = "".join([c for c in chain_res if c.isdigit()])
                    #get protonation state
                    #residue charge/state (supports legacy '+'/'--' and new '+1'/'-2' encodings)
                    _state_raw = result.residue_state
                    _charge_int = FORMAT.Parse_Charge_State(_state_raw)
                    proton_state = FORMAT.Format_Charge_State(_charge_int)
                    #check for flipped residues
                    is_flipped = result.residue_name.endswith("f") or ("_flip" in result.residue_name)
                    clash_flag = False
                    if self.config.edda:
                        cl = getattr(result,"clash",None)
                        if cl is not None and getattr(cl,"short_contact_indices",None):
                            scores = list(getattr(cl,"short_contact_indices",[]) or [])
                            labels = list(getattr(cl,"atom_labels",[]) or [])
                            rid = str(getattr(cl,"residue_id","") or "")
                            res_scores = []
                            #try label-based selection first ONLY if labels look residue-prefixed.
                            used_label_split = False
                            if rid and labels and len(labels) == len(scores):
                                for lab,v in zip(labels,scores):
                                    if isinstance(lab,str) and lab.startswith(rid + ":"):
                                        try:
                                            res_scores.append(float(v))
                                            used_label_split = True
                                        except Exception:
                                            pass
                            #if label-based produced nothing (or labels aren't prefixed), fall back to split_index.
                            if not used_label_split:
                                split = int(getattr(cl,"split_index",0) or 0)
                                split = min(max(split,0),len(scores))
                                for v in scores[split:]:
                                    try:
                                        res_scores.append(float(v))
                                    except Exception:
                                        pass
                            clash_flag = any(v <= -5.0 for v in res_scores)
                    display_name = base_name
                    if self.config.edda:
                        f.write(f"{display_name};{chain_id};{res_num};{proton_state};"
                                f"{display_name}{res_num};{is_flipped};{str(clash_flag).upper()};{result.binding_energy:.6f};"
                                f"{result.cap_correction:.6f}\n")
                    else:
                        f.write(f"{display_name};{chain_id};{res_num};{proton_state};"
                                f"{display_name}{res_num};{is_flipped};{result.binding_energy:.6f};"
                                f"{result.cap_correction:.6f}\n")
                #compute totals with a single, shared implementation.
                #(important for --no-cap where we must pick exactly one standard variant per residue.)
                total_interaction,_ = self._compute_total_interaction_energy_and_list_res(results)
                ligand_int_energies[lig_idx] += total_interaction
                f.write("\n")
        self.logger.info(f"Results written to {output_dir}")
        with open(output_dir / "stats.txt", "w") as f:
            for lig_idx,(ligand,results) in enumerate(zip(ligand_states,binding_results)):
                if (self.config.protonation_mode == ProtonationMode.RDKIT and ligand.smiles):
                    f.write(f"{ligand.smiles}:\n")
                else:
                    f.write(f"Ligand state {lig_idx}/{len(ligand_states)}: \n")
                IPA_binding_energy = 0.0
                #deformation contribution
                f.write(f"    Deformation Energy:          {ligand_states[lig_idx].deformation_energy:.2f} kcal/mol\n")
                IPA_binding_energy += ligand_states[lig_idx].deformation_energy
                #distortion (reported only)
                f.write(f"    Distortion:                  {ligand_states[lig_idx].distortion_energy:.2f} kcal/mol\n")
                #energy Chemical State (relative IP + acid–base correction)
                f.write(f"    Energy Chemical State:       {((ligand_states[lig_idx].energy_ip - ligand_min_energy + ligand_states[lig_idx].acid_base_correction)*self.Eh2kcalmol):.2f} kcal/mol\n")
                #IPA binding energy includes: deformation + (IP – min_IP + acid–base correction) + IPA interaction
                IPA_binding_energy += (ligand_states[lig_idx].energy_ip - ligand_min_energy + ligand_states[lig_idx].acid_base_correction)*self.Eh2kcalmol
                #IPA interaction energy
                f.write(f"    IPA Interaction Energy:      {ligand_int_energies[lig_idx]:.2f} kcal/mol\n")
                IPA_binding_energy += ligand_int_energies[lig_idx]
                #final IPA binding energy
                f.write(f"    IPA Binding Energy:          {IPA_binding_energy:.2f} kcal/mol\n")
                #final IPA model score
                Edist = ligand_states[lig_idx].distortion_energy
                if Edist < 1.0e-5: Edist = 1.0
                f.write(f"    IPA Model Score:             {(IPA_binding_energy/Edist):.2f}\n")
                #final IPA binding score
                Edef = ligand_states[lig_idx].deformation_energy
                if Edef < 1.0e-5: Edef = 1.0
                f.write(f"    IPA Binding Score:           {(IPA_binding_energy/Edef):.2f}\n")
                #raw QM chemical-state energy
                f.write(f"    Chemical State Total Energy: {ligand_states[lig_idx].energy_ip:.7f} Hartree\n")
                #RMSDs
                f.write(f"    RMSD IPA:                    {ligand_states[lig_idx].rmsd['ip']:.2f} Å\n")
                f.write(f"    RMSD full-relax:             {ligand_states[lig_idx].rmsd['opt']:.2f} Å\n")
                f.write("\n")
        #EDDA-only outputs
        if self.config.edda:
            try:
                self._write_clash_reports(ligand_states,binding_results)
            except Exception as e:
                self.logger.warning("Failed to write clash reports: %s",e)

    def _write_clash_reports(self, ligand_states: List[LigandState], binding_results: List[List[BindingResult]]) -> None:
        """Write EDDA-based clash reports under <output>/clash/.
           Creates:
             - clash/<ligand_name>.txt : human-readable per-ligand summary
             - clash/full_output.json  : machine-readable dump of all ClashAnalysis objects"""
        output_dir = self.config.output_folder
        clash_dir = output_dir / "clash"
        clash_dir.mkdir(parents = True,exist_ok = True)
        full_dump = {}
        for lig_idx,ligand in enumerate(ligand_states):
            lig_name = Path(ligand.file_path).stem
            #keep filenames tidy: strip trailing _ip if present (we add it in some places)
            lig_tag = lig_name.replace("_ip","")
            txt_path = clash_dir / f"{lig_tag}.txt"
            #collect per-residue entries
            per_residue = {}
            per_residue_dump = {}
            for result in binding_results[lig_idx]:
                cl = getattr(result,"clash",None)
                if cl is None:
                    continue
                split = int(getattr(cl,"split_index",0) or 0)
                labels = list(getattr(cl,"atom_labels",[]) or [])
                scores = list(getattr(cl,"short_contact_indices",[]) or [])
                #basic sanity, fall back to element list if labels missing
                if not labels or not scores or len(labels) != len(scores):
                    labels = list(getattr(cl,"heavy_atom_list",[]) or [])
                n = min(len(labels),len(scores))
                labels = labels[:n]
                scores = scores[:n]
                split = min(max(split,0),n)
                #split ligand/residue portions
                lig_scores = list(scores[:split])
                res_scores = list(scores[split:])
                lig_labels = list(labels[:split])
                res_labels = list(labels[split:])
                #element list (best effort): heavy_atom_list is aligned with short_contact_indices
                elems = list(getattr(cl,"heavy_atom_list",[]) or [])
                if len(elems) != len(scores):
                    elems = ["?"]*len(scores)
                lig_elems = elems[:split]

                def _to_float(x):
                    try:
                        return float(x)
                    except (TypeError,ValueError):
                        return None

                #build ligand entries: (element, index(1-based), value)
                lig_weird = []
                lig_subopt = []
                lig_clash = []
                for i,(el,lab,rawv) in enumerate(zip(lig_elems,lig_labels,lig_scores),start = 1):
                    v = _to_float(rawv)
                    if v is None: continue
                    if v >= 2.0:
                        lig_weird.append((el,i,v))
                    if (-5.0 <= v <= -2.0):
                        lig_subopt.append((el,i,v))
                    if v <= -5.0:
                        lig_clash.append((el,i,v))
                #build residue clash entries: (atomtype, value)
                #res_labels already look like "TYR_B253_1:OH" etc from pair code
                res_clash = []
                for lab,rawv in zip(res_labels,res_scores):
                    v = _to_float(rawv)
                    if v is None: continue
                    if v <= -5.0:
                        #atomtype is the part after "residue_id:"
                        atomtype = lab.split(":",1)[1] if ":" in lab else lab
                        res_clash.append((atomtype,v))
                if not (lig_weird or lig_subopt or lig_clash or res_clash):
                    continue
                #residue label for grouping (BindingResult name is fine)
                res_key = result.residue_name
                per_residue.setdefault(res_key,{"lig_weird": [],"lig_subopt": [],"lig_clash": [],"res_clash": []})
                per_residue[res_key]["lig_weird"].extend(lig_weird)
                per_residue[res_key]["lig_subopt"].extend(lig_subopt)
                per_residue[res_key]["lig_clash"].extend(lig_clash)
                per_residue[res_key]["res_clash"].extend(res_clash)
                #full dump (raw + labels)
                per_residue_dump[res_key] = {"residue_state": getattr(result,"residue_state",""),
                                             "binding_energy": getattr(result,"binding_energy",None),
                                             "cap_correction": getattr(result,"cap_correction",None),
                                             "split_index": split,
                                             "repulsion_energy": getattr(cl,"repulsion_energy",[]),
                                             "EDDA_binding_energy": getattr(cl,"EDDA_binding_energy",[]),
                                             "repulsion_contacts": getattr(cl,"repulsion_contacts",[]),
                                             "interaction_contacts": getattr(cl,"interaction_contacts",[]),
                                             "residue_id": getattr(cl,"residue_id",""),
                                             "heavy_atom_list": getattr(cl,"heavy_atom_list",[]),
                                             "short_contact_indices": getattr(cl,"short_contact_indices",[]),
                                             "atom_labels": getattr(cl,"atom_labels",[]),
                                             "split_index": getattr(cl,"split_index",0)}
            #write human-readable summary
            with open(txt_path,"w") as fh:
                fh.write(f"EDDA clash report for ligand: {lig_tag}\n")
                if ligand.smiles:
                    fh.write(f"SMILES: {ligand.smiles}\n")
                fh.write(f"Charge: {ligand.charge}\n\n")
                if not per_residue:
                    fh.write("No sub-optimal interactions, clashes, or weird behaviour detected.\n")
                else:
                    for res_key in sorted(per_residue.keys()):
                        fh.write(f"Residue: {res_key}\n")
                        #1) weird behaviour (ligand only)
                        items = per_residue[res_key]["lig_weird"]
                        fh.write(f"  Weird behaviour (score >= +2) (n={len(items)}):\n")
                        for el,idx,v in sorted(items,key = lambda x: x[2],reverse = True):
                            fh.write(f"    Ligand: {el:>2s}, {idx:>3d}, {v: .6f}\n")
                        #2) sub-optimal interactions (ligand only)
                        items = per_residue[res_key]["lig_subopt"]
                        fh.write(f"  Sub-optimal interactions (-5 <= score <= -2) (n={len(items)}):\n")
                        for el,idx,v in sorted(items,key = lambda x: x[2]):
                            fh.write(f"    Ligand: {el:>2s}, {idx:>3d}, {v: .6f}\n")
                        #3) clashes (ligand + protein)
                        lig_items = per_residue[res_key]["lig_clash"]
                        res_items = per_residue[res_key]["res_clash"]
                        fh.write(f"  Clashes (score <= -5) (n={len(lig_items) + len(res_items)}):\n")
                        for el,idx,v in sorted(lig_items,key = lambda x: x[2]):
                            fh.write(f"    Ligand:  {el:>2s}, {idx:>3d}, {v: .6f}\n")
                        for atomtype, v in sorted(res_items, key=lambda x: x[1]):
                            fh.write(f"    Protein: {atomtype:>4s}, {v: .6f}\n")
                        fh.write("\n")
            #store BOTH: (1) categorized human-friendly lists for CSV, (2) full raw dump for reference
            lig_charge = getattr(ligand_states[lig_idx],"charge","")
            full_dump[lig_tag] = {"ligand_index": lig_idx,
                                  "charge": lig_charge,
                                  "residues": per_residue,      # <- this is what the CSV flattener expects
                                  "raw": per_residue_dump}      # keep your detailed/raw content (not used for CSV)
        #global machine-readable dump (CSV, Excel-friendly)
        #NOTE: this is the *complete* per-atom EDDA score list (ligand atoms + residue atoms), not only the filtered "clash/suboptimal/weird" subsets.
        csv_path = clash_dir / "full_output.csv"
        import csv
        with open(csv_path,"w",newline = "",encoding = "utf-8") as fh:
            writer = csv.writer(fh)
            #header (no "category" on purpose)
            writer.writerow(["ligand_name",
                             "ligand_index",
                             "ligand_charge",
                             "residue_id",
                             "entity",          #ligand / protein
                             "atom_label",      #e.g. LIG:C1 or TYR_B253_1:OH
                             "atom_type",       #ligand element OR residue atomtype
                             "atom_index",      #1-based index within the entity (heavy atoms)
                             "value"])
            for lig_name,lig_data in full_dump.items():
                lig_idx = lig_data.get("ligand_index","")
                lig_charge = lig_data.get("charge","")
                #we use the "raw" section because it retains the full label/value list
                for res_id,dump in (lig_data.get("raw",{}) or {}).items():
                    labels = list(dump.get("atom_labels",[]) or [])
                    values = list(dump.get("short_contact_indices",[]) or [])
                    split = int(dump.get("split_index",0) or 0)
                    #best-effort element list aligned with values (often matches the ligand part)
                    elems = list(dump.get("heavy_atom_list",[]) or [])
                    if len(elems) != len(values):
                        elems = [""]*len(values)
                    n = min(len(labels),len(values))
                    split = min(max(split,0),n)
                    for i in range(n):
                        lab = labels[i]
                        try:
                            v = float(values[i])
                        except Exception:
                            continue
                        if i < split:
                            entity = "ligand"
                            atom_label = lab
                            atom_type = elems[i] or (lab.split(":",1)[1] if ":" in lab else "")
                            atom_index = i + 1
                        else:
                            entity = "protein"
                            atom_label = lab
                            atom_type = lab.split(":",1)[1] if ":" in lab else lab
                            atom_index = (i - split) + 1
                        writer.writerow([lig_name,lig_idx,lig_charge,res_id,entity,atom_label,atom_type,atom_index,v])

    def run_analysis(self) -> None:
        """Run the complete in-pocket analysis."""
        analysis_start = time.time()
        self.logger.info("Starting In-Pocket Analysis (%s)",IPA_VERSION)
        self.logger.info("="*60)
        native_ignored_msg = getattr(self.config,"native_ignored_flags_message","")
        if native_ignored_msg:
            self.logger.info(native_ignored_msg)
        try:
            #validate inputs
            step_start = time.time()
            self.validate_inputs()
            self._log_timing("Input validation",step_start)
            #read structure
            step_start = time.time()
            active_residue,protein_residues = self.read_structure()
            self._log_timing("Structure reading",step_start)
            #process ligand states
            step_start = time.time()
            ligand_states = self.process_ligand_states(active_residue)
            self._log_timing("Ligand state processing",step_start)
            #find neighbor residues
            step_start = time.time()
            neighbor_residues = self.find_neighbour_residues(active_residue,protein_residues,ligand_states)
            self._log_timing("Neighbor residue identification",step_start)
            if not neighbor_residues:
                self.logger.warning("No neighbor residues found within critical distance!")
                return
            #calculate deformation energies
            step_start = time.time()
            self.calculate_deformation_energies(ligand_states)
            self._log_timing("Deformation energy calculation",step_start)
            if not neighbor_residues:
                self.logger.warning("No neighbor residues found within critical distance!")
                return
            #process neighbor residues
            step_start = time.time()
            processed_residues = self.process_neighbor_residues(neighbor_residues)
            self._log_timing("Neighbor residue processing",step_start)
            #map processed residue variants by residue_id for later disambiguation
            #(used when --no-cap is active to select exactly one "standard" variant per residue)
            try:
                self._processed_residues_by_id = {getattr(r,"residue_id"): r
                                                  for r in (processed_residues or [])
                                                  if getattr(r,"residue_id",None)}
            except Exception:
                self._processed_residues_by_id = {}
            #calculate binding energies
            step_start = time.time()
            binding_results = self.calculate_binding_energies(ligand_states,processed_residues)
            self._log_timing("Binding energy calculation",step_start)
            #log interaction summary
            self._log_interaction_summary(ligand_states,binding_results)
            #generate per-ligand residue lists when requested
            # - no --best       -> True
            # - --best-debug    -> True
            # - --best only     -> False
            create_all_lists = (not self.config.best_mode) or self.config.best_debug
            if create_all_lists:
                try:
                    list_dir = self.config.output_folder / "list_residues"
                    list_dir.mkdir(parents = True,exist_ok = True)
                    for state in ligand_states:
                        ligand_name = Path(state.file_path).stem + "_ip"
                        ligand_charge = state.charge
                        csv_path = list_dir / f"{ligand_name}.csv"
                        self._write_residue_csv_for_ligand(ligand_name,ligand_charge,processed_residues,csv_path)
                except Exception as e:
                    self.logger.warning("Failed to generate per-ligand residue lists: %s",e)
            #generate per-ligand IP-optimised SDFs when requested
            # - no --best       -> True
            # - --best-debug    -> True
            # - --best only     -> False
            try:
                self._write_ligand_states_ip_sdfs(ligand_states,binding_results)
            except Exception as e:
                self.logger.warning("Failed to generate ligand_states_ip SDFs: %s",e)
            #write results
            step_start = time.time()
            self.write_results(ligand_states,binding_results)
            self._log_timing("Results writing",step_start)
            #final geometry cleanup to mirror serial layout
            self._final_geometry_cleanup()
            #if densities were calculated, convert all .cub -> .mrc
            if self.config.qm_params.calculate_density:
                step_start = time.time()
                self.logger.info("Cube to MRC: conversion started")
                try:
                    geometries_dir = self.config.output_folder / self.config.geometry_folder
                    cub_dir = self.config.output_folder / "cub_files"
                    mrc_dir = self.config.output_folder / "mrc_files"
                    cub_dir.mkdir(parents = True,exist_ok = True)
                    mrc_dir.mkdir(parents = True,exist_ok = True)
                    #Strategy: if any *_dens.cub are still in geometries/, move them into cub_files/ so they also get converted and geometries/ becomes clean.
                    if geometries_dir.exists():
                        for p in geometries_dir.rglob("*_dens.cub"):
                            try:
                                dest = cub_dir / p.name
                                if dest.exists():
                                    dest.unlink()
                                shutil.move(str(p),str(dest))
                            except Exception as e:
                                self.logger.warning("Failed to move stray cub from geometries/: %s",e)
                    #now onvert cubes found in cub_files/ to MRC files
                    n_converted = convert_cub_tree_to_mrc(cub_dir = cub_dir,mrc_dir = mrc_dir,move_cub_to_archive = False)  # cubes already live in cub_files/
                    #move any .mrc accidentally written under geometries/ into mrc_files/
                    if geometries_dir.exists():
                        for p in geometries_dir.rglob("*.mrc"):
                            try:
                                dest = mrc_dir / p.name
                                if dest.exists():
                                    dest.unlink()
                                shutil.move(str(p),str(dest))
                            except Exception as e:
                                self.logger.warning("Failed to move stray mrc from geometries/: %s",e)
                    #if we keep cubes, leave cub_files/ there; otherwise delete it
                    if not getattr(self.config,"keep_cub",False):
                        shutil.rmtree(cub_dir,ignore_errors = True)
                    self.logger.info("Cube to MRC: converted %d .cub files; mrc stored in %s",n_converted,mrc_dir)
                except Exception as e:
                    self.logger.warning("Cube to MRC conversion failed: %s",e)
                self._log_timing("Cube to MRC conversion",step_start)
            #optional best-ligand post-processing
            if self.config.best_mode:
                #generate compact best_* outputs
                #best_ligand uses the same Eint + residue-picking rules
                #implemented in this class (including --no-cap behavior).  
                try:
                    summary = generate_best_outputs(self.config.output_folder,
                                                    self.config.geometry_folder,
                                                    ligand_states,
                                                    binding_results,
                                                    logger = self.logger,
                                                    compute_total_interaction_and_list_res = self._compute_total_interaction_energy_and_list_res)
                except TypeError:
                    #backwards compatible with older best_ligand.generate_best_outputs signature
                    summary = generate_best_outputs(self.config.output_folder,
                                                    self.config.geometry_folder,
                                                    ligand_states,
                                                    binding_results,
                                                    logger = self.logger)
                #create best_list_residues.csv for the selected ligand
                try:
                    best_idx = summary.index
                    if 0 <= best_idx < len(ligand_states):
                        best_state = ligand_states[best_idx]
                        best_name = Path(best_state.file_path).stem + "_ip"
                        best_charge = best_state.charge
                        best_csv = self.config.output_folder / "best_list_residues.csv"
                        self._write_residue_csv_for_ligand(best_name,best_charge,processed_residues,best_csv)
                        #optional: also create best_list_residues_unique.csv (one entry per residue)
                        if getattr(self.config,"best_unique",False):
                            best_unique_csv = (self.config.output_folder / "best_list_residues_unique.csv")
                            self._write_best_unique_residue_csv(ligand_name = best_name,
                                                                ligand_charge = best_charge,
                                                                residues = processed_residues,
                                                                best_results = binding_results[best_idx],
                                                                target_file = best_unique_csv)
                except Exception as e:
                    self.logger.warning("Failed to write best residue list CSVs: %s",e)
                #log a short summary before (possibly) deleting logs
                try:
                    self.logger.info("[best_ligand] Best ligand index: %d "
                                     "(IPA Binding Energy: %.2f kcal/mol, Total interaction: %.2f kcal/mol)",
                                     summary.index + 1,
                                     summary.ipa_binding_energy,
                                     summary.total_interaction_energy)
                except Exception:
                    #fail-safe: don't let logging issues kill the run
                    pass
                #clean up outputs depending on debug flag
                #but keep geometries/ligand_states/ as a stable cache.
                cache_dir = self.config.output_folder / self.config.geometry_folder / "ligand_states"
                cache_tmp = self.config.output_folder / "_ligand_states_cache_tmp"
                moved_cache = False
                try:
                    #if we're going to delete stuff (best_mode + not best_debug),
                    #move the cache out of geometries/ before cleanup.
                    if self.config.best_mode and (not self.config.best_debug) and cache_dir.exists():
                        if cache_tmp.exists():
                            #remove any stale tmp cache from previous runs
                            shutil.rmtree(cache_tmp,ignore_errors = True)
                        cache_dir.rename(cache_tmp)
                        moved_cache = True
                except Exception as e:
                    self.logger.warning("Failed to protect ligand_states cache before cleanup: %s",e)
                #protect required geometries/ip/*.xyz referenced by best_list_residues*.csv
                #(only does something in compact --best, not in --best-debug)
                self._protect_and_restore_best_ip_xyz()
                cleanup_outputs(self.config.output_folder,
                                geometry_folder = self.config.geometry_folder,
                                output_folder_qm = self.config.output_folder_qm,
                                keep_debug = self.config.best_debug,
                                logger = self.logger)
                #restore the protected IP xyz files after cleanup
                try:
                    tmp_keep = getattr(self,"_ip_keep_tmp",None)
                    if tmp_keep and Path(tmp_keep).exists():
                        ip_dir = (self.config.output_folder / self.config.geometry_folder / "ip")
                        ip_dir.mkdir(parents = True,exist_ok = True)
                        for xyz in Path(tmp_keep).glob("*.xyz"):
                            shutil.copy2(xyz,ip_dir / xyz.name)
                        shutil.rmtree(tmp_keep,ignore_errors = True)
                        delattr(self,"_ip_keep_tmp")
                except Exception as e:
                    self.logger.warning("Failed to restore protected IP xyz: %s",e)
                #restore the cache after cleanup
                try:
                    if moved_cache:
                        cache_dir.parent.mkdir(parents = True,exist_ok = True)
                        if cache_dir.exists():
                            shutil.rmtree(cache_dir,ignore_errors = True)
                        cache_tmp.rename(cache_dir)
                except Exception as e:
                    self.logger.warning("Failed to restore ligand_states cache after cleanup: %s",e)
            #log total analysis time
            total_time = time.time() - analysis_start
            #strict failure policy: if --edda is present and any EDDA pair fails, the run FAILED.
            bc = getattr(self,"binding_calculator",None)
            failed = bool(getattr(bc,"analysis_failed",False))
            pair_fail = int(getattr(bc,"pair_failures",0) or 0)
            edda_fail = int(getattr(bc,"edda_failures",0) or 0)
            if failed:
                #drop a marker file that pipelines can detect
                try:
                    fail_path = self.config.output_folder / "FAILED.flag"
                    with open(fail_path,"w") as fh:
                        fh.write("In-Pocket Analysis FAILED\n")
                        fh.write(f"pair_failures = {pair_fail}\n")
                        fh.write(f"edda_failures = {edda_fail}\n")
                        fh.write("See output/debug/ for EDDA failure bundles.\n")
                except Exception:
                    pass
                self.logger.error("="*60)
                self.logger.error("Analysis finished with failures in %.2f seconds (pair_failures = %d, edda_failures = %d). "
                                  "See output/debug/ for details.",total_time,pair_fail,edda_fail)
                raise InPocketAnalyzerError(f"Analysis finished with failures (pair_failures = {pair_fail}, edda_failures = {edda_fail}). "
                    f"See {self.config.output_folder}/debug/ and FAILED.flag.")
            self.logger.info("="*60)
            self.logger.info(f"Analysis completed successfully in {total_time:.2f} seconds")
        except InPocketAnalyzerError as e:
            self.logger.error(f"Analysis failed: {e}")
            raise
        except Exception as e:
            self.logger.error(f"Unexpected error: {e}",exc_info = True)
            raise InPocketAnalyzerError(f"Unexpected error during analysis: {e}")

def create_argument_parser() -> argparse.ArgumentParser:
    """Create and configure the argument parser."""
    parser = argparse.ArgumentParser(
        description = "In-Pocket Analysis: Protein-ligand interaction analysis tool",
        formatter_class = argparse.RawDescriptionHelpFormatter,
        epilog = """Examples:
                    %(prog)s --pdb protein.pdb --active-residue 123 --chain A
                    %(prog)s --pdb complex.pdb --active-residue 456 --distance 8.0 --ph 6.5
                    %(prog)s --pdb structure.cif --active-residue LIG --protonation rdkit --tautomers""")
    #required arguments
    parser.add_argument("--pdb",
                        "--pdb-file",
                        "--pdbfile",
                        type = str,
                        required = True,
                        help = "Input PDB or mmCIF file")
    parser.add_argument("--active-residue",
                        "--a-residue",
                        "--aresidue",
                        "--ares",
                        type = int,
                        default = None,
                        help = "Active residue number")
    parser.add_argument("--active-residue-name",
                        "--a-residue-name",
                        "--aresidue-name",
                        "--aresn",
                        type = str,
                        default = None,
                        help = "Active residue name")
    #optional file arguments
    parser.add_argument("--output",
                        "--output-folder",
                        type = str,
                        default = "analysis_output",
                        help = "Output folder (default: analysis_output)")
    parser.add_argument("--png",
                        action = "store_true",
                        help = "Keep RDKit 2D PNG previews (low resolution) in <output_folder>/png/ (default: delete)")
    parser.add_argument("--png-hq",
                        action = "store_true",
                        help = "Keep RDKit 2D PNG previews (high resolution) in <output_folder>/png/ (implies --png)")
    parser.add_argument("--use-ip-structure",
                        "--useipstructure",
                        "--use-ip-struct",
                        "--useipstruct",
                        action = "store_true",
                        help = "Use in-pocket optimized structures for all optimisations (default: False)")
    parser.add_argument("--use-ip-residues",
                        "--useipresidues",
                        "--use-ip-res",
                        "--useipres",
                        action = "store_true",
                        help = "Use in-pocket optimized residues for all optimisations (default: False)")
    parser.add_argument("--chain",
                        "--active-chain",
                        "--a-chain",
                        "--achain",
                        type = str,
                        default = "A",
                        help = "Chain ID of active residue (default: A)")
    #analysis parameters
    parser.add_argument("--distance",
                        "--critical-distance",
                        "--within",
                        "--maxdist",
                        "--rmax",
                        type = float,
                        default = 5.0,
                        help = "Critical distance for neighbor search in Å (default: 5.0)")
    parser.add_argument("--charge",
                        "--chrg",
                        "--chr",
                        type = int,
                        default = None,
                        help = "Active-moiety charge (default: 0)")
    parser.add_argument("--force-user-charge",
                        action = "store_true",
                        help = ("If a backend-defined active-moiety charge disagrees with --charge, "
                                "use the user-provided charge directly and skip the interactive prompt."))
    parser.add_argument("--ph",
                        type = float,
                        default = 7.0,
                        help = "pH for protonation (default: 7.0)")
    #QM parameters
    parser.add_argument("--alpha",
                        type = float,
                        default = 1.0,
                        help = "Alpha parameter for in-pocket (default: 1.0)")
    parser.add_argument("--kappa",
                        type = float,
                        default = 0.025,
                        help = "Kappa parameter for in-pocket (default: 0.025)")
    parser.add_argument("--etemperature",
                        "--electronic-temperature",
                        "--etemp",
                        "--ETemp",
                        type = float,
                        default = 300.0,
                        help = "Electronic temperature in K (default: 300.0); you can use this to increase the multireference character of the system, thus simulate better transition metals")
    #protonation options
    parser.add_argument("--protonation",
                        choices=["rdkit","obabel","vsepr","ccd","original","native"],
                        default = "rdkit",
                        help = ("Protonation algorithm (default: rdkit). Use 'original'/'native' to preserve "
                                "the input chemical state, skip protonation/tautomer generation, and keep "
                                "existing hydrogens; capping is still applied, and residue termini may be "
                                "completed with VSEPR only to add missing cap hydrogens."))
    parser.add_argument("--ccd-file",
                        "--ccdfile",
                        "--ccd",
                        type = str,
                        default = None,
                        help = "CCD file for CCD protonation")
    #processing options
    parser.add_argument("--no-cap",
                        "--nocap",
                        action = "store_true",
                        help = "Don't cap residue termini")
    parser.add_argument("--include-water",
                        "--include-h2o",
                        "--include-hoh",
                        action = "store_true",
                        help = "Include water molecules in analysis")
    parser.add_argument("--no-flip",
                        "--noflip",
                        action = "store_true",
                        help = "Don't flip residues")
    parser.add_argument("--no-flip-HIS",
                        "--noflipHIS",
                        action = "store_true",
                        help = "Don't flip HIS residues")
    parser.add_argument("--no-flip-ASN",
                        "--noflipASN",
                        action = "store_true",
                        help = "Don't flip ASN residues")
    parser.add_argument("--no-flip-GLN",
                        "--noflipGLN",
                        action = "store_true",
                        help = "Don't flip GLN residues")
    parser.add_argument("--no-flip-SER",
                        "--noflipSER",
                        action = "store_true",
                        help = "Don't flip SER residues")
    parser.add_argument("--no-flip-THR",
                        "--noflipTHR",
                        action = "store_true",
                        help = "Don't flip THR residues")
    parser.add_argument("--no-flip-TYR",
                        "--noflipTYR",
                        action = "store_true",
                        help = "Don't flip TYR residues")
    parser.add_argument("--no-flip-NA",
                        "--noflipNA",
                        action = "store_true",
                        help = "Don't flip nucleic acid residues (phosphates)")
    parser.add_argument("--tautomers",
                        "--tautomerise",
                        "--tautomerize",
                        "--get-tautomers",
                        "--gettautomers",
                        action = "store_true",
                        help = "Generate and analyze tautomers")
    parser.add_argument("--max-tautomers",
                        type = int,
                        default = 10,
                        help = "Maximum number of tautomers (default: 10)")
    parser.add_argument("--check-residue-structures",
                        "--check-residues",
                        "--checkresidues",
                        action = "store_true",
                        help = "Pause the run to check the structures of protonated neighbouring residues")
    #advanced options
    parser.add_argument("--no-solvation",
                        action = "store_true",
                        help = "Disable implicit solvation")
    parser.add_argument("--solvent",
                        type = str,
                        default = "water",
                        help = "Solvent for implicit solvation (default: water)")
    parser.add_argument("--edda",
                        "--qm-clash-analysis",
                        "--qm-clash-score",
                        "--qm-clash",
                        action = "store_true",
                        help = "perform EDDA based QM analysis")
    parser.add_argument("--energy-threshold",
                        type = float,
                        default = 5.0e-6,
                        help = "Energy convergence threshold (default: 5.0e-6 Eh)")
    parser.add_argument("--gradient-threshold",
                        type = float,
                        default = 5.0e-3,
                        help = "Gradient convergence threshold (default: 5.0e-3 Eh/A)")
    parser.add_argument("--altloc",
                        "--alternate-location",
                        "--alt-loc",
                        type = str,
                        default = None,
                        help = "Alternative location identifier (default: A)")
    parser.add_argument("--no-parallel",
                        action = "store_true",
                        help = "Disable parallel processing")
    parser.add_argument("--max-workers",
                        type = int,
                        help = "Maximum number of parallel workers")
    parser.add_argument("--inclh",
                        "--inclH",
                        "--include-hydrogens",
                        "--includehydrogens",
                        action = "store_true",
                        help = "Whether to include hydrogen atoms in the bias potential (default: false)")
    parser.add_argument("--keep-input-h",
                        action = "store_true",
                        help = "Do NOT strip hydrogens already present in the input PDB/mmCIF (default: strip input H for RDKit robustness).")
    parser.add_argument("--keep-external-caps",
                        action = "store_true",
                        help = ("Do NOT strip external/artificial caps already present in the input structure "
                                "(e.g., ACE/NMA/NME). Default: strip only if they collide with a real residue id."))
    parser.add_argument("--noip",
                        "--no-ip",
                        "--no-in-pocket",
                        "--no-inpocket",
                        action = "store_true",
                        help = "Instead of in-pocket, do unconstrained optimisation")
    parser.add_argument("--eval",
                        "--evaluate",
                        "--evaluate-only",
                        action = "store_true",
                        help = "Evaluate ligand-residue interactions only, without any sort of optimisation")
    parser.add_argument("--skipdef",
                        "--skip-def",
                        "--skip-deformation",
                        action = "store_true",
                        help = "Whether to completely skip the calculation of deformation energy")
    parser.add_argument("--badcs",
                        "--badCS",
                        "--bad_chemical_states",
                        "--bad_chem_states",
                        action = "store_true",
                        help = "Whether to avoid skipping the bad chemical states defined by RDKit")
    parser.add_argument("--calculate-density",
                        "--calc-density",
                        "--calc-dens",
                        "--calcdens",
                        action = "store_true",
                        help = "Whether to calculate and store electronic densities")
    parser.add_argument("--keep-cap-uncorrected",
                        "--keep-cap-uncorr",
                        "--keepcapuncorr",
                        action = "store_true",
                        help = "Whether to remove the capping corrections")
    parser.add_argument("--ligand-topology",
                        "--ligand-topol",
                        "--lig_topol",
                        "--ligtopol",
                        type = str,
                        default = None,
                        help = "The user can specify a certain base topology, either from an SDF file or a SMILES. Accepts link to .sdf file, link to .smi file or a SMILES string")
    #best-ligand compact output modes
    parser.add_argument("--best",
                        action = "store_true",
                        help = ("Generate only best_interaction_energies.csv and best_lig.sdf as "
                                "compact outputs (other files/folders are removed at the end)."))
    parser.add_argument("--best-debug",
                        action = "store_true",
                        help = ("Generate best_interaction_energies.csv and best_lig.sdf but keep "
                                "all standard outputs (logs, geometries, QM files) for debugging."))
    parser.add_argument("--best-unique",
                        action = "store_true",
                        help = ("When used together with --best/--best-debug, also write "
                                "best_list_residues_unique.csv (one entry per residue, selecting the "
                                "best protonation/flip variant by most favorable interaction energy)."))
    #density options
    parser.add_argument("--cub",
                        "--keep-cub",
                        action = "store_true",
                        help = ("Keep original .cub density files. When used with --calculate-density, "
                                "cubes are stored in <output>/cub_files and are NOT deleted at the end."))
    #output output options
    parser.add_argument("--debug",
                        action = "store_true",
                        help = "Enable debug logging (very verbose per-residue/per-cap details).")
    parser.add_argument("--verbose",
                        "-v",
                        action = "store_true",
                        help = "Enable verbose logging")
    return parser

def parse_arguments_to_config(args: argparse.Namespace) -> AnalysisConfig:
    """Convert parsed arguments to AnalysisConfig object."""
    config = AnalysisConfig(pdb_file = Path(args.pdb))
    #checking that alternative arguments are given
    if (args.active_residue is None) and (args.active_residue_name is None):
        print("\nEither --active-residue or --active-residue-name must be provided",file = sys.stderr)
        sys.exit(1)
    #file paths
    config.pdb_file = Path(args.pdb)
    config.output_folder = Path(args.output)
    config.active_residue_chain = args.chain
    config.keep_png = bool(getattr(args,"png",False) or getattr(args,"png_hq",False))
    config.keep_png_hq = bool(getattr(args,"png_hq",False))
    #ensure the output path exists
    if (not config.output_folder.exists()) or (not config.output_folder.is_dir()):
        config.output_folder.mkdir(parents = True,exist_ok = True)
    config.remove_cap_effect = not args.keep_cap_uncorrected
    #active residue identifiers
    config.active_residue = args.active_residue
    config.active_residue_name = args.active_residue_name
    #analysis parameters
    config.critical_distance = args.distance
    config.ph = args.ph
    #QM parameters
    qm = config.qm_params
    object.__setattr__(qm,"alpha",args.alpha)
    object.__setattr__(qm,"kappa",args.kappa)
    object.__setattr__(qm,"electronic_temperature",args.etemperature)
    object.__setattr__(qm,"ignore_hydrogens",not args.inclh)
    #input sanitization flags
    config.keep_input_h = bool(getattr(args,"keep_input_h",False))
    config.debug = bool(getattr(args,"debug",False))
    config.keep_external_caps = bool(getattr(args,"keep_external_caps",False))
    #processing options
    config.cap_residues = not args.no_cap
    config.skip_water = not args.include_water
    config.flip_residues = not args.no_flip
    config.flip_ASN = not args.no_flip_ASN
    config.flip_GLN = not args.no_flip_GLN
    config.flip_HIS = not args.no_flip_HIS
    config.flip_SER = not args.no_flip_SER
    config.flip_THR = not args.no_flip_THR
    config.flip_TYR = not args.no_flip_TYR
    config.flip_NA = not args.no_flip_NA
    config.tautomerise = args.tautomers
    config.check_residue_structures = args.check_residue_structures
    config.max_tautomers = args.max_tautomers
    config.ligand_topology = args.ligand_topology
    config.edda = args.edda
    #parallelisation options
    config.parallel_processing = not args.no_parallel
    config.max_workers = args.max_workers
    #if parallel is enabled and user did NOT specify max_workers, use half of the available cores (minimum 1) by default
    if config.parallel_processing and config.max_workers is None:
        try:
            n_cores = os.cpu_count() or 1
        except Exception:
            n_cores = 1
        config.max_workers = max(1,n_cores//2)
    #QM options
    solvent_name = SOLV.SolventList(args.solvent)
    use_solvation = not args.no_solvation
    if solvent_name == "gas":
        use_solvation = False
    #update QMParameters
    object.__setattr__(qm,"use_solvation",use_solvation)
    object.__setattr__(qm,"solvent_name",solvent_name)
    object.__setattr__(qm,"energy_threshold",args.energy_threshold)
    object.__setattr__(qm,"gradient_threshold",args.gradient_threshold)
    object.__setattr__(qm,"calculate_density",args.calculate_density)
    #density output options
    config.keep_cub = bool(getattr(args,"cub",False))
    #charge-related options stay on config (they’re not part of QMParameters)
    config.charge = args.charge
    config.user_defined_charge = args.charge is not None
    config.force_user_charge = bool(getattr(args,"force_user_charge",False))
    if config.charge is None:
        config.charge = 0
    config.use_ip_structure = args.use_ip_structure
    config.use_ip_residues = args.use_ip_residues
    config.do_in_pocket = not args.noip
    config.evaluate_only = args.eval
    config.exclude_rdkit_bad_structures = not args.badcs
    config.skip_deformation_energy = args.skipdef
    #protonation
    prot_arg = args.protonation
    if prot_arg == "native":
        prot_arg = "original"
    config.protonation_mode = ProtonationMode(prot_arg)
    if config.protonation_mode == ProtonationMode.ORIGINAL:
        native_ignored_flags = []
        if args.tautomers:
            native_ignored_flags.append("--tautomers")
        if args.ligand_topology not in (None,""):
            native_ignored_flags.append("--ligand-topology")
        if native_ignored_flags:
            config.native_ignored_flags_message = (f"Ignoring {', '.join(native_ignored_flags)} because --protonation {args.protonation} "
                                                   f"preserves the original/native input chemistry.")
        config.keep_input_h = True
        config.flip_residues = False
        config.tautomerise = False
    config.ccd_file = args.ccd_file
    if args.protonation == "ccd" and (args.ccd_file is None):
        print("\nIf you give --protonation ccd then you must provide a valid --ccd-file",file = sys.stderr)
        sys.exit(1)
    #file format
    config.altloc = args.altloc
    if config.altloc is None:
        config.altloc = "A"
    config.check_altloc = True
    #best-ligand post-processing options
    #any of --best or --best-debug turns on best_mode;
    #best_debug controls whether we keep the full output tree.
    config.best_mode = bool(getattr(args,"best",False) or getattr(args,"best_debug",False))
    config.best_debug = bool(getattr(args,"best_debug",False))
    #only meaningful if best_mode is enabled, but harmless otherwise
    config.best_unique = bool(getattr(args,"best_unique",False))
    return config

def main():
    """Main entry point for the script."""
    parser = create_argument_parser()
    args = parser.parse_args()
    #set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    else:
        logging.getLogger().setLevel(logging.INFO)
    try:
        #convert arguments to configuration
        config = parse_arguments_to_config(args)
        #create and run analyzer
        analyzer = InPocketAnalyzer(config)
        #run analysis
        analyzer.run_analysis()
    except InPocketAnalyzerError as e:
        print(f"Error: {e}",file = sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("\nAnalysis interrupted by user",file = sys.stderr)
        sys.exit(130)
    except Exception as e:
        print(f"Unexpected error: {e}",file = sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
