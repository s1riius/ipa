"""In-Pocket extension to include complete chemical environments"""

import argparse
import logging
import copy
import sys
import os
import shutil
import math
import subprocess
from dataclasses import dataclass,field
from typing import List,Tuple
from enum import Enum
from pathlib import Path
import numpy as np
import util.BaseAminoAcids as AA
import util.SolventList as SOLVENT
import util.Formatting as FORMAT
import util.Processing as PROCESS
import util.Geometry as GEOM
import util.QMwrapper as QM
import util.Protoro as PROTORO
import util.Tautoro as TAUTORO
import util.VSEPR as VSEPR
import util.Scoring as SCORE
import util.Monomers as MONOMER
import util.UnitConversion as UC
import util.Residues as RES

class ProtonationMode(Enum):
    """Enumeration for protonation algorithms."""
    RDKIT = "rdkit"
    OBABEL = "obabel"
    VSEPR = "vsepr"
    CCD = "ccd"

class CalculationMode(Enum):
    """Enumeration for calculation modes."""
    NOOPT = "no-opt"
    IPA = "ipa"
    OPT = "opt"

Eh_to_kcalmol = UC.Hartree2Kcalmol

def parse_arguments(argv = None) -> argparse.Namespace:
    """Parse command-line arguments for fragment analysis."""
    parser = argparse.ArgumentParser(
        description = "In-Pocket fragment mode for analysis of a protein-ligand system.",
        formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    #core inputs
    parser.add_argument(
        "--pdb",
        dest = "pdb_file",
        required = True,
        help = "Input PDB/mmCIF file containing the protein-ligand complex.")
    parser.add_argument(
        "--active-residue-a",
        dest = "active_residues_a",
        required = True,
        help = "Ligand residue(s), e.g., 'A501' or 'A501,B502'.")
    parser.add_argument(
        "--active-residue-b",
        dest = "active_residues_b",
        required = True,
        help = "Pocket residue(s), e.g., 'A70,A71,A75'.")
    #file/storage control
    parser.add_argument(
        "--save-file",
        dest = "save_file",
        default = "protein_cut",
        help = "Template name to store the output of the calculations.")
    parser.add_argument(
        "--summary-file",
        dest = "summary_file",
        default = "summary_file.csv",
        help = "Summary file name.")
    parser.add_argument(
        "--output-folder",
        dest = "output_folder",
        default = "",
        help = "Optional prefix folder in which to write geometries and output.")
    #optimisation thresholds
    parser.add_argument(
        "--ethresh",
        type = float,
        default = 5.0e-6,
        help = "Energy threshold for convergence, in Eh (forwarded as ethresh/energythresh).")
    parser.add_argument(
        "--gthresh",
        type = float,
        default = 5.0e-3,
        help = "Gradient threshold for convergence in Eh/A (forwarded as gthresh/gradthresh).")
    #connectivity/capping/intermediates
    parser.add_argument(
        "--no-cap",
        dest = "cap",
        action = "store_false",
        help = "Disable capping of residue termini (default: capped).")
    parser.set_defaults(cap = True)
    parser.add_argument(
        "--no-connect-residues",
        dest = "connect_residues",
        action = "store_false",
        help = "Do NOT connect residues when fragmenting (default: connect).")
    parser.set_defaults(connect_residues = True)
    parser.add_argument(
        "--no-intermediate-pdb",
        dest = "write_interm_pdb",
        action = "store_false",
        help = "Disable writing of intermediate PDB structures.")
    parser.set_defaults(write_interm_pdb = True)
    #in-pocket parameters
    parser.add_argument(
        "--alpha",
        type = float,
        default = 1.0,
        help = "Alpha parameter for in-pocket optimisation.")
    parser.add_argument(
        "--kappa",
        type = float,
        default = 0.025,
        help = "Kappa parameter for in-pocket optimisation.")
    parser.add_argument(
        "--ignore-protons",
        dest = "ignore_protons",
        action = "store_true",
        help = "Ignore protons (legacy ignh = 1).")
    parser.add_argument(
        "--include-protons",
        dest = "ignore_protons",
        action = "store_false",
        help = "Include protons explicitly (legacy ignh = 0).")
    parser.set_defaults(ignore_protons = True)
    parser.add_argument(
        "--calculation-mode",
        choices = ["no-opt","ipa","opt"],
        default = "ipa",
        help = "Type of optimisation ran: no optimisation of the ligand-pocket complex (no-opt); restrained optimisation of the complex (ipa); free, unconstrained optimisation of the complex (opt).")
    #protonation and pH
    parser.add_argument(
        "--ph",
        type = float,
        default = 7.0,
        help = "pH for protonation.")
    parser.add_argument(
        "--ph-spread",
        type = float,
        default = 4.0,
        help = "pH spread for protonation (pH +/- pH-spread).")
    parser.add_argument(
        "--protonation-mode",
        choices = ["obabel","rdkit","vsepr","ccd"],
        default = "obabel",
        help = "Algorithm used for protonation of the structures.")
    parser.add_argument(
        "--tautomerise",
        dest = "tautomerise",
        action = "store_true",
        help = "Whether to search for tautomers.")
    parser.add_argument(
        "--max-tautomers",
        type = int,
        default = 10,
        help = "Maximum number of conformers considered.")
    parser.add_argument(
        "--max-natom-tautomers",
        type = int,
        default = 3,
        help = "Maximum number of ATOM residue types allowed for tautomer search.")
    #flipping of ASN/GLN/HIS residues
    parser.add_argument(
        "--flip-residues",
        dest = "flip_residues",
        action = "store_true",
        help = "Enable flipping of residues (ASN/GLN/HIS) where appropriate.")
    parser.add_argument(
        "--no-flip-residues",
        dest = "flip_residues",
        action = "store_false",
        help = "Disable flipping of residues.")
    parser.set_defaults(flip_residues = True)
    parser.add_argument(
        "--flip-asn",
        dest = "flip_asn",
        action = "store_true",
        help = "Allow flipping of ASN residues.")
    parser.add_argument(
        "--no-flip-asn",
        dest = "flip_asn",
        action = "store_false",
        help = "Disallow flipping of ASN residues.")
    parser.set_defaults(flip_asn = False)
    parser.add_argument(
        "--flip-gln",
        dest = "flip_gln",
        action = "store_true",
        help = "Allow flipping of GLN residues.")
    parser.add_argument(
        "--no-flip-gln",
        dest = "flip_gln",
        action = "store_false",
        help = "Disallow flipping of GLN residues.")
    parser.set_defaults(flip_gln = False)
    parser.add_argument(
        "--flip-his",
        dest = "flip_his",
        action = "store_true",
        help = "Allow flipping of HIS residues.")
    parser.add_argument(
        "--no-flip-his",
        dest = "flip_his",
        action = "store_false",
        help = "Disallow flipping of HIS residues.")
    parser.set_defaults(flip_his = False)
    #solvent/solvation
    parser.add_argument(
        "--use-solvent",
        dest = "use_solvent",
        action = "store_true",
        help = "Enable ALPB solvation.")
    parser.add_argument(
        "--no-solvent",
        dest = "use_solvent",
        action = "store_false",
        help = "Disable ALPB solvation.")
    parser.set_defaults(use_solvent = True)
    parser.add_argument(
        "--solvent",
        dest = "solvent_name",
        default = "water",
        help = "Solvent to use (water, chloroform, methanol, etc.).")
    #charge
    parser.add_argument(
        "--charge-a",
        type = int,
        default = 0,
        help = "Charge for subsystem A.")
    parser.add_argument(
        "--charge-b",
        type = int,
        default = 0,
        help = "Charge for subsystem B.")
    #density
    parser.add_argument(
        "--calculate-density",
        dest = "calculate_density",
        action = "store_true",
        help = "Request density calculation.")
    parser.add_argument(
        "--no-calculate-density",
        dest = "calculate_density",
        action = "store_false",
        help = "Do not calculate density.")
    parser.set_defaults(calculate_density = False)
    #binding/scoring
    parser.add_argument(
        "--no-ebind",
        dest = "calculate_ebind",
        action = "store_false",
        help = "Do not calculate binding energy between fragments.")
    parser.set_defaults(calculate_ebind = True)
    parser.add_argument(
        "--no-edda",
        dest = "calculate_edda",
        action = "store_false",
        help = "Do not calculate EDDA and Short-Contact ccore.")
    parser.set_defaults(calculate_edda = True)
    #electronic temperature
    parser.add_argument(
        "--telec",
        type = float,
        default = 300.0,
        help = "Electronic temperature in Kelvin (Telec).")
    #CCD/altloc options
    parser.add_argument(
        "--ccd-file",
        dest = "ccd_file",
        default = "",
        help = "Path to the CCD CIF file used for building structures.")
    parser.add_argument(
        "--use-ccd",
        dest = "use_ccd",
        action = "store_true",
        help = "Enable use of CCD information when building structures.")
    parser.add_argument(
        "--no-use-ccd",
        dest = "use_ccd",
        action = "store_false",
        help = "Disable use of CCD information when building structures.")
    parser.set_defaults(use_ccd = False)
    #altLoc selection
    parser.add_argument(
        "--select-altloc",
        dest = "altloc",
        choices = ["A", "B", "C", "none"],
        default = "A",
        help = "AltLoc identifier to select from PDB/mmCIF (A, B, C, or none).")
    #verbosity
    parser.add_argument(
        "--verbose",
        action = "store_true",
        help = "Enable verbose logging.")
    return parser.parse_args(argv)

@dataclass
class AnalysisConfig:
    """Configuration class for In-Pocket Analysis fragment mode parameters."""
    #arguments and default-like values,
    pdb_file:            Path = field(default_factory = Path)
    geom_folder:         str = "geometries"
    out_folder:          str = "output"
    output_folder:       Path = field(default_factory = lambda: Path.cwd() / "fragment_analysis")
    geometry_folder:     str = "geometries"
    output_file_folder:  str = "output"
    savefile:            str = "protein_cut"
    summaryfile:         str = "summary_file.csv"
    mmcif_file:          bool = False
    #basic processing and parameters
    cap_residues:        bool = True
    connectresidues:     bool = True
    writeIntermPDB:      bool = True
    pH:                  float = 7.0
    pH_spread:           float = 4.0
    protonation_mode:    ProtonationMode = ProtonationMode.OBABEL
    tautomerise:         bool = False
    max_tautomers:       int = 10
    max_natom_tautomers: int = 3
    #bias potential options
    alpha: float = 1.0
    kappa: float = 0.025
    ignh:  bool = True
    calculation_mode:    CalculationMode = CalculationMode.IPA
    #QM options and thresholds
    Telec:              float = 300.0
    energy_threshold:   float = 5.0e-6
    gradient_threshold: float = 5.0e-3
    calculate_density:  bool = False
    #flipping
    flipresidues: bool = True
    flipASN:      bool = True
    flipGLN:      bool = True
    flipHIS:      bool = True
    #residues of interest
    activeresA: str = ""
    activeresB: str = ""
    #charge
    chargeA: int = 0
    chargeB: int = 0
    #solvation
    solvation:  bool = True
    solventname: str = "water"
    #energetics
    calculate_Ebind: bool = True
    calculate_EDDA:  bool = True
    #CCD-related
    ccdfile: str = ""
    useccd:  bool = False
    #altlocs
    altloc:      str = "A"
    checkaltloc: bool = False

    def __post_init__(self):
        """Post-initialization processing."""
        #convert paths to Path objects
        if isinstance(self.pdb_file,str):
            self.pdb_file = Path(self.pdb_file)
        if isinstance(self.output_folder,str):
            self.output_folder = Path(self.output_folder)

@dataclass
class FragmentState:
    """Represents a fragment chemical state."""
    file_path_xyz: List[Path] = field(default_factory = list)
    file_path_sdf: List[Path] = field(default_factory = list)
    file_path_xyz_H: List[Path] = field(default_factory = list)
    file_path_pdb: List[Path] = field(default_factory = list)
    file_path_pdb_H: List[Path] = field(default_factory = list)
    file_path_ip: List[Path] = field(default_factory = list)
    file_path_ebind: List[Path] = field(default_factory = list)
    sdf_files_are_available: bool = False
    number_HETATM: List[int] = field(default_factory = list)
    number_ATOM: List[int] = field(default_factory = list)
    charge: List[int] = field(default_factory = lambda: [0,0,0])
    energy_ip: List[float] = field(default_factory = lambda: [0.0,0.0,0.0])
    RMSD: List[float] = field(default_factory = lambda: [0.0,0.0,0.0])
    converged_ip: List[bool] = field(default_factory = lambda: [False,False,False])
    total_polarisability: List[float] = field(default_factory = lambda: [0.0,0.0,0.0])
    atom_list: List[int] = field(default_factory = lambda: [0,0,0])
    atomic_polarisabilities: List[float] = field(default_factory = lambda: [0.0,0.0,0.0])
    atomic_charges: List[float] = field(default_factory = lambda: [0.0,0.0,0.0])
    repulsion_energy: float = 0.0
    EDDA_binding_energy: float = 0.0
    heavy_atom_list: List[int] = field(default_factory = lambda: [0,0,0])
    repulsion_contacts: List[float] = field(default_factory = lambda: [0,0,0])
    interaction_contacts: List[float] = field(default_factory = lambda: [0,0,0])
    short_contact_indices: List[float] = field(default_factory = lambda: [0,0,0])

def parse_arguments_to_config(args: argparse.Namespace) -> AnalysisConfig:
    """Convert parsed arguments to AnalysisConfig object."""
    config = AnalysisConfig()
    #arguments and default-like values,
    config.pdb_file      = Path(args.pdb_file)
    config.geom_folder   = Path("geometries/")
    config.out_folder    = Path("output/")
    config.output_folder = Path(args.output_folder or "fragment_analysis/")
    config.output_folder.mkdir(parents = True,exist_ok = True)
    config.savefile      = Path(args.save_file)
    config.summaryfile   = Path(args.summary_file)
    config.mmcif_file    = False  # will be set based on file extension later
    config.geometry_folder = config.output_folder/config.geom_folder
    config.output_file_folder = config.output_folder/config.out_folder
    #basic processing and parameters
    config.cap_residues        = args.cap
    config.connectresidues     = args.connect_residues
    config.writeIntermPDB      = args.write_interm_pdb
    config.pH                  = args.ph
    config.pH_spread           = args.ph_spread
    config.protonation_mode    = ProtonationMode(args.protonation_mode)
    config.tautomerise         = args.tautomerise
    config.max_tautomers       = args.max_tautomers
    config.max_natom_tautomers = args.max_natom_tautomers
    #bias potential options
    config.alpha              = args.alpha
    config.kappa              = args.kappa
    config.ignh               = 1 if args.ignore_protons else 0
    config.calculation_mode   = CalculationMode(args.calculation_mode)
    #system charge
    config.chargeA = args.charge_a
    config.chargeB = args.charge_b
    #QM options and thresholds
    config.Telec              = args.telec
    config.energy_threshold   = args.ethresh
    config.gradient_threshold = args.gthresh
    config.calculate_density  = 1 if args.calculate_density else 0
    #flipping
    config.flipresidues       = args.flip_residues
    config.flipASN            = args.flip_asn
    config.flipGLN            = args.flip_gln
    config.flipHIS            = args.flip_his
    #residues of interest
    config.activeresA     = args.active_residues_a
    config.activeresB     = args.active_residues_b
    #checking that alternative arguments are given
    if (config.activeresA is None) or (config.activeresB is None) or (config.activeresA == "") or (config.activeresB == ""):
        print("\nAll active residues must be defined",file = sys.stderr)
        sys.exit(1)
    #solvation
    config.solvation      = 1 if args.use_solvent else 0
    config.solventname    = SOLVENT.SolventList(args.solvent_name)
    #energetics
    config.calculate_Ebind  = args.calculate_ebind
    config.calculate_EDDA   = args.calculate_edda
    #CCD-related
    config.ccdfile       = args.ccd_file
    config.useccd        = 1 if args.use_ccd else 0
    if args.protonation_mode == "ccd" and (args.ccd_file is None):
        print("\nIf you give --protonation ccd then you must provide a valid --ccd-file",file = sys.stderr)
        sys.exit(1)
    #altlocs
    config.altloc        = args.altloc if args.altloc != "none" else "A"
    config.checkaltloc   = (args.altloc != "none")
    return config

class FragmentAnalyzerError(Exception):
    """Custom exception for InPocketAnalyzer errors."""
    pass

class FragmentAnalyzer:
    """Main class for performing fragment analysis."""
    run_from_previous_IPA: bool
    run_from_geometry: bool
    path_to_IPA_geometries: str
    set_from_IPA: str
    ligand_name_from_IPA: str
    residue_capping: str
    ligand_charge_from_IPA: int
    residue_states_from_IPA: List[List[str]] = field(default_factory = list)
    extension_from_IPA: str

    def __init__(self, config: AnalysisConfig):
        self.config = config
        self.logger = self._setup_logging()
        self._setup_directories()
        self.run_from_previous_IPA = False
        self.run_from_geometry = False
        self.path_to_IPA_geometries = ""
        self.set_from_IPA = ""
        self.ligand_name_from_IPA = ""
        self.ligand_charge_from_IPA = 0
        self.Eh2kcalmol = Eh_to_kcalmol
        self.minimum_energy = float("inf")  #track minimum energy
    
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
            log_file = self.config.output_folder / "fragment_analysis.log"
            file_handler = logging.FileHandler(log_file)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
            logger.setLevel(logging.INFO)
        return logger

    def _setup_directories(self) -> None:
        """Create necessary directories for the analysis."""
        directories = [self.config.output_folder,
                       self.config.geometry_folder,
                       self.config.output_file_folder,
                       self.config.geometry_folder / "start_xyz",
                       self.config.geometry_folder / "start_sdf",
                       self.config.geometry_folder / "start_pdb",
                       self.config.geometry_folder / "ip",
                       self.config.geometry_folder / "opt"]
        for directory in directories:
            directory.mkdir(parents = True,exist_ok = True)
            self.logger.debug(f"Created directory: {directory}")

    def validate_inputs(self) -> None:
        """Validate input parameters and files."""
        self.logger.info("Validating inputs...")
        #check PDB file
        if not self.config.pdb_file or not self.config.pdb_file.exists():
            #try with extensions
            if (self.config.pdb_file.suffix != ".pdb") and (self.config.pdb_file.suffix != ".cif") and (self.config.pdb_file.suffix != ".mmcif"):
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
                    raise FragmentAnalyzerError(f"Structure file not found: {self.config.pdb_file}")
        #determine file format
        if self.config.pdb_file.suffix.lower() in [".cif",".mmcif"]:
            self.config.mmcif_file = True
            self.logger.info("Detected mmCIF format")
        #validate parameters
        if self.config.alpha <= 0:
            raise FragmentAnalyzerError("Alpha parameter must be positive")
        if not 0 <= self.config.pH <= 14:
            self.logger.warning(f"pH value {self.config.pH} is outside typical range [0-14]")
        self.logger.info("Input validation completed successfully")

    def _restructure_residues(self, geometry, atoms, atmtp, residues, resnumbr, chainid, restp, seqidcodes, pdblines, terminal) -> Tuple[List[RES.ResidueData],List[Tuple[int,int,str]]]:
        """Transform a full list of residues into an organised dictionary."""
        return RES.restructure_residues(geometry,atoms,atmtp,residues,resnumbr,chainid,restp,seqidcodes,pdblines,terminal)

    def read_structure(self) -> Tuple[List[RES.ResidueData],List[Tuple[int,int,str]]]:
        """Read and parse the input structure file."""
        self.logger.info(f"Reading structure from {self.config.pdb_file}")
        try:
            if self.config.mmcif_file:
                geometry,atoms,atmtp,residues,resnumbr,chainid,restp,seqidcodes,pdblines = (GEOM.ReadCIFStructure(str(self.config.pdb_file)))
                try:
                    geometry,atoms,atmtp,residues,resnumbr,chainid,restp,seqidcodes,pdblines = (
                        GEOM.SelectOneAltLocOnly(geometry,atoms,atmtp,residues,resnumbr,chainid,restp,seqidcodes,pdblines,self.config.altloc))
                    altloc_warnings = GEOM.GetAltLocSelectionWarnings(pdblines)
                    for msg in altloc_warnings:
                        self.logger.warning(msg)
                except Exception as e:
                    self.logger.warning(f"Failed to process altloc: {e}")
            else:
                geometry,atoms,atmtp,residues,resnumbr,chainid,restp,seqidcodes,pdblines = (GEOM.ReadPDB_base(str(self.config.pdb_file),self.config.altloc))
        except Exception as e:
            raise FragmentAnalyzerError(f"Failed to read structure file: {e}")
        #process altloc if needed
        if self.config.checkaltloc and self.config.mmcif_file:
            try:
                geometry,atoms,atmtp,residues,resnumbr,chainid,restp,pdblines = (GEOM.SelectOneAltLocOnly(geometry,atoms,atmtp,residues,resnumbr,chainid,restp,seqidcodes,pdblines,self.config.altloc))
            except Exception as e:
                self.logger.warning(f"Failed to process altloc: {e}")
        #calculate termini information
        terminal = ["N"]
        for idx in range(1,len(resnumbr)):
            if (resnumbr[idx] != resnumbr[idx - 1]) and (resnumbr[idx - 1] + 1 != resnumbr[idx]):
                terminal.append("N")
            else:
                terminal.append("0")
        for idx in range(1,len(resnumbr)):
            if terminal[idx] == "N":
                terminal[idx - 1] = "C"
        #restructure residue data
        residues, non_prot_data = self._restructure_residues(geometry,atoms,atmtp,residues,resnumbr,chainid,restp,seqidcodes,pdblines,terminal)
        return residues,non_prot_data

    def get_active_res_type(self, listactiveres, ires_other, totalnumberofresidues) -> List[str]:
        self.logger.info(f"Getting the types of active residues")
        return RES.GetActiveResType(listactiveres,ires_other,totalnumberofresidues)

    def get_connectors(self, listactiveres, activerestype) -> List[Tuple[int,str]]:
        self.logger.info(f"Extracting connecting residues")
        return RES.GetConnectors(listactiveres,activerestype)

    def ProcessActiveResidues(self, activeres) -> List[Tuple[int,str]]:
        """function that processes and parses the input residues into a list that can be used to extract what matters."""
        self.logger.info(f"Getting the active residue list")
        return PROCESS.process_active_residues_base(activeres)
    
    def FetchSystemFromIPA(self, activeresA, activeresB) -> Tuple[List[Tuple[int,str]],List[Tuple[int,str]]]:
        """Getting the pocket information from previous lists/calculations (IPA residue csv)."""
        self.logger.info(f"Fetching pocket information from IPA run.")
        datafile = activeresA
        if (not activeresA.endswith(".txt")) and (not activeresA.endswith(".csv")):
            datafile = activeresB
        #defaults for new-format CSVs
        self.run_from_previous_IPA = True
        self.residue_states_from_IPA = []
        self.path_to_IPA_geometries = ""
        self.extension_from_IPA = ".xyz"
        self.set_from_IPA = ""              #may be absent in new csv
        self.residue_capping = False        #cap: may be absent in new csv
        self.ligand_name_from_IPA = ""
        self.ligand_charge_from_IPA = 0
        with open(datafile, "r") as rfile:
            for raw in rfile:
                line = raw.strip()
                if (not line) or line.startswith("#"):
                    continue
                #some examples contain "..." placeholders -> ignore safely
                if line == "..." or line.startswith("..."):
                    continue
                if line.startswith("path-to-geometries:"):
                    self.path_to_IPA_geometries = line.replace("path-to-geometries:","").strip()
                    continue
                if line.startswith("extension:"):
                    self.extension_from_IPA = line.replace("extension:","").strip()
                    continue
                if line.startswith("set:"):
                    self.set_from_IPA = line.replace("set:","").strip()
                    continue
                if line.startswith("cap:"):
                    capval = line.replace("cap:","").strip().lower()
                    self.residue_capping = (capval in ("true","1","yes","y"))
                    continue
                if line.startswith("ligand:"):
                    #new format: ligand: Z9999_prot_2_tautomer_2_ip -2
                    data = line.replace("ligand:","").strip().split()
                    if len(data) >= 1:
                        self.ligand_name_from_IPA = data[0]
                    if len(data) >= 2:
                        try:
                            self.ligand_charge_from_IPA = int(data[1])
                        except Exception:
                            self.ligand_charge_from_IPA = 0
                    continue
                if line.startswith("r:"):
                    #new format: r:MET;A;69;0;MET_A69_ip
                    parts = line.replace("r:","").strip().split(";")
                    #normalize to 5 fields: [RES, CHAIN, RESNUM, STATE, BASENAME]
                    # - tolerate extra columns (keep first 4 + last as basename)
                    # - tolerate missing basename (best effort)
                    if len(parts) >= 5:
                        resname = parts[0].strip()
                        chain = parts[1].strip()
                        resnum = parts[2].strip()
                        state = parts[3].strip()
                        basename = parts[-1].strip()
                        self.residue_states_from_IPA.append([resname,chain,resnum,state,basename])
                    elif len(parts) == 4:
                        resname = parts[0].strip()
                        chain = parts[1].strip()
                        resnum = parts[2].strip()
                        state = parts[3].strip()
                        basename = f"{resname}_{chain}{resnum}_ip"
                        self.residue_states_from_IPA.append([resname,chain,resnum,state,basename])
                    else:
                        #unexpected line -> skip safely
                        continue
        #ensure geometry path ends with "/"
        if self.path_to_IPA_geometries and (not self.path_to_IPA_geometries.endswith("/")):
            self.path_to_IPA_geometries += "/"
        #build the list of active residues
        #ligand "active residue" is mostly dummy for this mode, but keep previous behaviour.
        list_active_res_A = [[9999,"A"]]
        if self.ligand_name_from_IPA:
            try:
                if (self.ligand_name_from_IPA[0].isalpha()) and (self.ligand_name_from_IPA[1:].isdigit()):
                    list_active_res_A = [[int(self.ligand_name_from_IPA[1:]),self.ligand_name_from_IPA[0]]]
                elif (self.ligand_name_from_IPA[-1].isalpha()) and (self.ligand_name_from_IPA[:-1].isdigit()):
                    list_active_res_A = [[int(self.ligand_name_from_IPA[:-1]),self.ligand_name_from_IPA[-1]]]
            except Exception:
                pass
        list_active_res_B = []
        for residue in self.residue_states_from_IPA:
            try:
                list_active_res_B.append([int(residue[2]),residue[1]])
            except Exception:
                continue
        return list_active_res_A,list_active_res_B

    def get_fragments_from_residues(self, residues, non_protein_residue_data) -> Tuple[List[Tuple[int,str]],List[Tuple[int,str]],List[Tuple[int,str]],List[Tuple[int,str]]]:
        """Parse the input structure into fragments."""
        self.logger.info(f"Splitting structure into fragments")
        #get the active residues as a list
        list_of_activeres_A = []
        list_of_activeres_B = []
        connectorsA = []
        connectorsB = []
        skip_logging = False
        if (self.config.activeresA.endswith(".txt") or self.config.activeresA.endswith(".csv")) or (self.config.activeresB.endswith(".txt") or self.config.activeresB.endswith(".csv")):
            list_of_activeres_A,list_of_activeres_B = self.FetchSystemFromIPA(self.config.activeresA,self.config.activeresB)
        elif (self.config.activeresA.endswith(".xyz") or self.config.activeresA.endswith(".pdb") or self.config.activeresA.endswith(".sdf")) and (self.config.activeresB.endswith(".xyz") or self.config.activeresB.endswith(".pdb") or self.config.activeresB.endswith(".sdf")):
            list_of_activeres_A = self.config.activeresA
            list_of_activeres_B = self.config.activeresB
            self.run_from_geometry = True
            skip_logging = True
        else:
            list_of_activeres_A = self.ProcessActiveResidues(self.config.activeresA)
            list_of_activeres_B = self.ProcessActiveResidues(self.config.activeresB)
            totalnumberofresidues = len(residues)
            #get the type of the active residues
            activerestypeA = self.get_active_res_type(list_of_activeres_A,non_protein_residue_data,totalnumberofresidues)
            activerestypeB = self.get_active_res_type(list_of_activeres_B,non_protein_residue_data,totalnumberofresidues)
            #connect the residues?
            connectorsA = self.get_connectors(list_of_activeres_A,activerestypeA) if self.config.connectresidues else []
            connectorsB = self.get_connectors(list_of_activeres_B,activerestypeB) if self.config.connectresidues else []
        #logging
        if not skip_logging:
            self.logger.info(f"list of selected active residues for block A:")
            for iactA in list_of_activeres_A:
                self.logger.info(f"    {iactA[1]}{iactA[0]}")
            self.logger.info(f"list of connectors for block A:")
            for iactA in connectorsA:
                self.logger.info(f"    {iactA[1]}{iactA[0]}")
            self.logger.info(f"list of selected active residues for block B:")
            for iactB in list_of_activeres_B:
                self.logger.info(f"    {iactB[1]}{iactB[0]}")
            self.logger.info(f"list of connectors for block B:")
            for iactB in connectorsB:
                self.logger.info(f"    {iactB[1]}{iactB[0]}")
        return list_of_activeres_A,connectorsA,list_of_activeres_B,connectorsB

    def ExtractSubSystem(self, residues, list_of_activeres, connectors, cap_residues, flip_HIS, flip_ASN, flip_GLN, connecting_residues, identifier) -> Tuple[List[RES.ResidueData],List[int],List[int],List[int]]:
        """Function that extracts part of a system of interest."""
        self.logger.info(f"Extracting subsystem " + identifier)
        #allocate the data
        new_residues: List[RES.ResidueData] = []
        #loop over residues
        icounter = 0
        active_site_size = len(list_of_activeres) + len(connectors)
        for idx in range(len(residues)):
            res_x = residues[idx]
            identifier_x = [res_x.number,res_x.chain_id]
            #add whole residue
            if (identifier_x in list_of_activeres):
                if cap_residues and (idx > 0) and (res_x.terminal != "N") and (res_x.residue_type != "HETATM"):
                    idw = idx - 1
                    res_w = residues[idw]
                    identifier_w = [res_w.number,res_w.chain_id]
                    if (not (identifier_w in list_of_activeres)) and (not (identifier_w in connectors)):
                        for idw_atm in range(len(res_w.atoms)):
                            atom_tp_w = res_w.atom_types[idw_atm]
                            if (atom_tp_w == "C") or (atom_tp_w == "O"):
                                res_x.cap_atoms.append(res_w.atoms[idw_atm])
                                res_x.cap_geometry.append(res_w.geometry[idw_atm])
                new_residues.append(res_x)
                icounter += 1
            #add only main chain atoms
            elif connecting_residues and (identifier_x in connectors):
                new_geom = []
                new_atoms = []
                new_atom_types = []
                new_pdb_lines = []
                for idx_atm in range(len(res_x.atoms)):
                    atom_tp_x = res_x.atom_types[idx_atm]
                    if (atom_tp_x == "C") or (atom_tp_x == "O") or (atom_tp_x == "N") or (atom_tp_x == "CA"):
                        new_geom.append(res_x.geometry[idx_atm])
                        new_atoms.append(res_x.atoms[idx_atm])
                        new_atom_types.append(res_x.atom_types[idx_atm])
                        new_pdb_lines.append(res_x.pdb_lines[idx_atm])
                new_residue = RES.ResidueData(geometry = new_geom,
                                              atoms = new_atoms,
                                              atom_types = new_atom_types,
                                              name = res_x.name,
                                              number = res_x.number,
                                              chain_id = res_x.chain_id,
                                              residue_type = res_x.residue_type,
                                              pdb_lines = new_pdb_lines,
                                              charge = res_x.charge,
                                              terminal = res_x.terminal,
                                              residue_id = res_x.residue_id,
                                              protonation_state = res_x.protonation_state,
                                              processed_atoms = res_x.processed_atoms,
                                              processed_geometry = res_x.processed_geometry)
                new_residues.append(new_residue)
                icounter += 1
            if icounter == active_site_size:
                break
        #check for flipping residues
        list_ASN_flip = []
        list_GLN_flip = []
        list_HIS_flip = []
        if flip_HIS or flip_ASN or flip_GLN:
            for idx, res in enumerate(new_residues):
                if res.name == "ASN" and flip_ASN:
                    list_ASN_flip.append(idx)
                elif res.name == "GLN" and flip_GLN:
                    list_GLN_flip.append(idx)
                elif res.name == "HIS" and flip_HIS:
                    list_HIS_flip.append(idx)
        return new_residues,list_ASN_flip,list_GLN_flip,list_HIS_flip

    def assemble_system(self, residues, identifier) -> Tuple[List[Tuple[float,float,float]],List[str],List[str],List[str],List[int],List[str],List[str],float]:
        self.logger.info(f"Assembling geometry {identifier} from the collection of residues")
        return RES.AssembleSystem(residues,False,str(self.config.geometry_folder) + "/start_pdb/")

    def _flip_residue(self, residue: RES.ResidueData) -> RES.ResidueData:
        """Flip ASN, GLN, or HIS residue according to original logic."""
        import copy
        flipped = copy.deepcopy(residue)
        if residue.name == "ASN":
            #swap OD1 with ND2 or OD2 with ND1
            for i, atp in enumerate(residue.atom_types):
                if atp == "OD1":
                    flipped.atom_types[i] = "ND2"
                    flipped.atoms[i] = "N"
                    coordinate_to_find = str("{:.3f}".format(flipped.geometry[i][2]))
                    string_to_replace = flipped.pdb_lines[i][flipped.pdb_lines[i].rfind(coordinate_to_find):]
                    new_string = string_to_replace.replace("O","N")
                    flipped.pdb_lines[i] = flipped.pdb_lines[i].replace("OD1","ND2").replace(string_to_replace,new_string)
                elif atp == "OD2":
                    flipped.atom_types[i] = "ND1"
                    flipped.atoms[i] = "N"
                    coordinate_to_find = str("{:.3f}".format(flipped.geometry[i][2]))
                    string_to_replace = flipped.pdb_lines[i][flipped.pdb_lines[i].rfind(coordinate_to_find):]
                    new_string = string_to_replace.replace("O","N")
                    flipped.pdb_lines[i] = flipped.pdb_lines[i].replace("OD2","ND1").replace(string_to_replace,new_string)
                elif atp == "ND1":
                    flipped.atom_types[i] = "OD2"
                    flipped.atoms[i] = "O"
                    coordinate_to_find = str("{:.3f}".format(flipped.geometry[i][2]))
                    string_to_replace = flipped.pdb_lines[i][flipped.pdb_lines[i].rfind(coordinate_to_find):]
                    new_string = string_to_replace.replace("N","O")
                    flipped.pdb_lines[i] = flipped.pdb_lines[i].replace("ND1","OD2").replace(string_to_replace,new_string)
                elif atp == "ND2":
                    flipped.atom_types[i] = "OD1"
                    flipped.atoms[i] = "O"
                    coordinate_to_find = str("{:.3f}".format(flipped.geometry[i][2]))
                    string_to_replace = flipped.pdb_lines[i][flipped.pdb_lines[i].rfind(coordinate_to_find):]
                    new_string = string_to_replace.replace("N","O")
                    flipped.pdb_lines[i] = flipped.pdb_lines[i].replace("ND2","OD1").replace(string_to_replace,new_string)
        elif residue.name == "GLN":
            #swap OE1 with NE2 or OE2 with NE1
            for i, atp in enumerate(residue.atom_types):
                if atp == "OE1":
                    flipped.atom_types[i] = "NE2"
                    flipped.atoms[i] = "N"
                    coordinate_to_find = str("{:.3f}".format(flipped.geometry[i][2]))
                    string_to_replace = flipped.pdb_lines[i][flipped.pdb_lines[i].rfind(coordinate_to_find):]
                    new_string = string_to_replace.replace("O","N")
                    flipped.pdb_lines[i] = flipped.pdb_lines[i].replace("OE1","NE2").replace(string_to_replace,new_string)
                elif atp == "OE2":
                    flipped.atom_types[i] = "NE1"
                    flipped.atoms[i] = "N"
                    coordinate_to_find = str("{:.3f}".format(flipped.geometry[i][2]))
                    string_to_replace = flipped.pdb_lines[i][flipped.pdb_lines[i].rfind(coordinate_to_find):]
                    new_string = string_to_replace.replace("O","N")
                    flipped.pdb_lines[i] = flipped.pdb_lines[i].replace("OE2","NE1").replace(string_to_replace,new_string)
                elif atp == "NE1":
                    flipped.atom_types[i] = "OE2"
                    flipped.atoms[i] = "O"
                    coordinate_to_find = str("{:.3f}".format(flipped.geometry[i][2]))
                    string_to_replace = flipped.pdb_lines[i][flipped.pdb_lines[i].rfind(coordinate_to_find):]
                    new_string = string_to_replace.replace("N","O")
                    flipped.pdb_lines[i] = flipped.pdb_lines[i].replace("NE1","OE2").replace(string_to_replace,new_string)
                elif atp == "NE2":
                    flipped.atom_types[i] = "OE1"
                    flipped.atoms[i] = "O"
                    coordinate_to_find = str("{:.3f}".format(flipped.geometry[i][2]))
                    string_to_replace = flipped.pdb_lines[i][flipped.pdb_lines[i].rfind(coordinate_to_find):]
                    new_string = string_to_replace.replace("N","O")
                    flipped.pdb_lines[i] = flipped.pdb_lines[i].replace("NE2","OE1").replace(string_to_replace,new_string)
        elif residue.name == "HIS":
            #complex HIS flipping, 4 atoms swapped
            for i, atp in enumerate(residue.atom_types):
                if atp == "ND1":
                    flipped.atom_types[i] = "CD2"
                    flipped.atoms[i] = "C"
                    coordinate_to_find = str("{:.3f}".format(flipped.geometry[i][2]))
                    string_to_replace = flipped.pdb_lines[i][flipped.pdb_lines[i].rfind(coordinate_to_find):]
                    new_string = string_to_replace.replace("N","C")
                    flipped.pdb_lines[i] = flipped.pdb_lines[i].replace("ND1","CD2").replace(string_to_replace,new_string)
                elif atp == "ND2":
                    flipped.atom_types[i] = "CD1"
                    flipped.atoms[i] = "C"
                    coordinate_to_find = str("{:.3f}".format(flipped.geometry[i][2]))
                    string_to_replace = flipped.pdb_lines[i][flipped.pdb_lines[i].rfind(coordinate_to_find):]
                    new_string = string_to_replace.replace("N","C")
                    flipped.pdb_lines[i] = flipped.pdb_lines[i].replace("ND2","CD1").replace(string_to_replace,new_string)
                elif atp == "CD1":
                    flipped.atom_types[i] = "ND2"
                    flipped.atoms[i] = "N"
                    coordinate_to_find = str("{:.3f}".format(flipped.geometry[i][2]))
                    string_to_replace = flipped.pdb_lines[i][flipped.pdb_lines[i].rfind(coordinate_to_find):]
                    new_string = string_to_replace.replace("C","N")
                    flipped.pdb_lines[i] = flipped.pdb_lines[i].replace("CD1","ND2").replace(string_to_replace,new_string)
                elif atp == "CD2":
                    flipped.atom_types[i] = "ND1"
                    flipped.atoms[i] = "N"
                    coordinate_to_find = str("{:.3f}".format(flipped.geometry[i][2]))
                    string_to_replace = flipped.pdb_lines[i][flipped.pdb_lines[i].rfind(coordinate_to_find):]
                    new_string = string_to_replace.replace("C","N")
                    flipped.pdb_lines[i] = flipped.pdb_lines[i].replace("CD2","ND1").replace(string_to_replace,new_string)
                elif atp == "NE1":
                    flipped.atom_types[i] = "CE2"
                    flipped.atoms[i] = "C"
                    coordinate_to_find = str("{:.3f}".format(flipped.geometry[i][2]))
                    string_to_replace = flipped.pdb_lines[i][flipped.pdb_lines[i].rfind(coordinate_to_find):]
                    new_string = string_to_replace.replace("N","C")
                    flipped.pdb_lines[i] = flipped.pdb_lines[i].replace("NE1","CE2").replace(string_to_replace,new_string)
                elif atp == "NE2":
                    flipped.atom_types[i] = "CE1"
                    flipped.atoms[i] = "C"
                    coordinate_to_find = str("{:.3f}".format(flipped.geometry[i][2]))
                    string_to_replace = flipped.pdb_lines[i][flipped.pdb_lines[i].rfind(coordinate_to_find):]
                    new_string = string_to_replace.replace("N","C")
                    flipped.pdb_lines[i] = flipped.pdb_lines[i].replace("NE2","CE1").replace(string_to_replace,new_string)
                elif atp == "CE1":
                    flipped.atom_types[i] = "NE2"
                    flipped.atoms[i] = "N"
                    coordinate_to_find = str("{:.3f}".format(flipped.geometry[i][2]))
                    string_to_replace = flipped.pdb_lines[i][flipped.pdb_lines[i].rfind(coordinate_to_find):]
                    new_string = string_to_replace.replace("N","C")
                    flipped.pdb_lines[i] = flipped.pdb_lines[i].replace("CE1","NE2").replace(string_to_replace,new_string)
                elif atp == "CE2":
                    flipped.atom_types[i] = "NE1"
                    flipped.atoms[i] = "N"
                    coordinate_to_find = str("{:.3f}".format(flipped.geometry[i][2]))
                    string_to_replace = flipped.pdb_lines[i][flipped.pdb_lines[i].rfind(coordinate_to_find):]
                    new_string = string_to_replace.replace("N","C")
                    flipped.pdb_lines[i] = flipped.pdb_lines[i].replace("CE2","NE1").replace(string_to_replace,new_string)
        return flipped
    
    def extract_whole_connectivity(self, residues) -> Tuple[List[Tuple[float,float,float]],List[str],List[str],List[Tuple[str,str,str,bool]]]:
        """Function extracting the connectivity for a whole block"""
        self.logger.info(f"Generating fragment connectivity")
        connectivity = []
        geometry_fragment = []
        atoms_fragment = []
        atomtypes_fragment = []
        for ires in range(len(residues)):
            residue = residues[ires]
            skip_side_chain = (residue.name != "GLY") and (len(residue.atoms) == 4)
            for iatm in range(len(residue.cap_atoms)):
                if residue.cap_atoms[iatm] == "C": residue.cap_atoms[iatm] = "C*"
                elif residue.cap_atoms[iatm] == "O": residue.cap_atoms[iatm] = "O*"
            terminal_res = ["0"]
            if ires + 1 == 1:
                if len(residue.cap_atoms) != 0: terminal_res.append("N-cap")
                else: terminal_res.append("N-term")
            if ires + 1 > 1:
                residue_prev = residues[ires - 1]
                if (residue_prev.number + 1 == residue.number) and (residue_prev.residue_type == "ATOM") and (residue.residue_type == "ATOM"):
                    terminal_res.append("connect prev")
                    #this can lead to issues with non-standard aa, like SEP and so on
                if (residue_prev.number + 1 != residue.number):
                    if len(residue.cap_atoms) != 0: terminal_res.append("N-cap")
                    else: terminal_res.append("N-term")
            if ires + 1 < len(residues):
                residue_next = residues[ires + 1]
                if (residue_next.number == residue.number + 1) and (residue_next.residue_type == "ATOM") and (residue.residue_type == "ATOM"):
                    terminal_res.append("connect next")
                if (residue_next.number - 1 != residue.number):
                    terminal_res.append("C-cap")
            if ires + 1 == len(residues):
                if "OXT" in residue.atom_types: terminal_res.append("C-term")
                else: terminal_res.append("C-cap")
            for iatm in range(len(residue.atoms)):
                if residue.atoms[iatm] == "H": continue
                geometry_fragment.append(residue.geometry[iatm])
                atoms_fragment.append(residue.atoms[iatm])
                atomtypes_fragment.append(residue.atom_types[iatm] + "_" + residue.chain_id + str(residue.number))
            for iatm in range(len(residue.cap_atoms)):
                if residue.cap_atoms[iatm] == "H": continue
                geometry_fragment.append(residue.cap_geometry[iatm])
                atomtypes_fragment.append(residue.cap_atoms[iatm] + "_" + residue.chain_id + str(residue.number))
                atoms_fragment.append(residue.cap_atoms[iatm].replace("*",""))
            local_connectivity = []
            if not AA.IsBaseAminoAcid(residue.name):
                local_connectivity = MONOMER.GeneralisedExtractCCDconnectivityGEMMI(str(self.config.pdb_file),residue.name)
            else: 
                local_connectivity = AA.BaseAminoAcidConnectivity(residue.name,terminal_res)
            #basic connectivity
            for pair in local_connectivity:
                atomA = pair[0]
                atomB = pair[1]
                side_chain_A = not AA.main_chain_atom(atomA)
                side_chain_B = not AA.main_chain_atom(atomB)
                if skip_side_chain and (side_chain_A or side_chain_B): continue
                bond_order = pair[2]
                aromatic = pair[3]
                connectivity.append([atomA + "_" + residue.chain_id + str(residue.number),atomB + "_" + residue.chain_id + str(residue.number),bond_order,aromatic])
            #terminal information: caps and other termini are done, so now we need to connect residues
            if "connect prev" in terminal_res:
                #connect N to C
                connectivity.append(["N_" + residue.chain_id + str(residue.number),"C_" + residues[ires - 1].chain_id + str(residues[ires - 1].number),"sing",False])
        return geometry_fragment, atoms_fragment, atomtypes_fragment, connectivity

    def count_hetatm_and_atom(self, residues, identifier) -> List[int]:
        """Function that counts the number of HETATM and ATOM residues in a residue list"""
        self.logger.info(f"Counting the number of ATOM and HETATM in residue list {identifier}")
        n_atom = 0
        n_hetatm = 0
        for residue in residues:
            if residue.residue_type.lower() == "atom": n_atom += 1
            elif residue.residue_type.lower() == "hetatm": n_hetatm += 1
        return [n_atom,n_hetatm]

    def generate_fragment_geometries(self, residues, list_of_activeres_A, connectorsA, list_of_activeres_B, connectorsB) -> List[FragmentState]:
        """Generating fragment geometry from the available data"""
        self.logger.info(f"Generating fragment structures")
        #base files to store data
        xyz_files = [[self.config.geometry_folder / f"{self.config.savefile}_AB",
                      self.config.geometry_folder / f"{self.config.savefile}_A",
                      self.config.geometry_folder / f"{self.config.savefile}_B"]]
        #to speed up, check whether there is anything to flip
        flip_HIS_A = self.config.flipHIS
        flip_ASN_A = self.config.flipASN
        flip_GLN_A = self.config.flipGLN
        connect_res_A = self.config.connectresidues
        if len(list_of_activeres_A) + len(connectorsA) == 1:
            #assume this is a ligand
            flip_HIS_A = False
            flip_ASN_A = False
            flip_GLN_A = False
            connect_res_A = False
        flip_HIS_B = self.config.flipHIS
        flip_ASN_B = self.config.flipASN
        flip_GLN_B = self.config.flipGLN
        connect_res_B = self.config.connectresidues
        if len(list_of_activeres_B) + len(connectorsB) == 1:
            #assume this is a ligand
            flip_HIS_B = False
            flip_ASN_B = False
            flip_GLN_B = False
            connect_res_B = False
        #collect all the residues of interest for the system
        residues_A,list_ASN_flip_A,list_GLN_flip_A,list_HIS_flip_A = self.ExtractSubSystem(residues,list_of_activeres_A,connectorsA,self.config.cap_residues,flip_HIS_A,flip_ASN_A,flip_GLN_A,connect_res_A,"A")
        residues_B,list_ASN_flip_B,list_GLN_flip_B,list_HIS_flip_B = self.ExtractSubSystem(residues,list_of_activeres_B,connectorsB,self.config.cap_residues,flip_HIS_B,flip_ASN_B,flip_GLN_B,connect_res_B,"B")
        types_list_A = self.count_hetatm_and_atom(residues_A,"A")
        types_list_B = self.count_hetatm_and_atom(residues_B,"B")
        #merge list of residues to flip
        list_ASN_flip = [*list_ASN_flip_A,*list_ASN_flip_B]
        list_GLN_flip = [*list_GLN_flip_A,*list_GLN_flip_B]
        list_HIS_flip = [*list_HIS_flip_A,*list_HIS_flip_B]
        list_flips_A = [*list_ASN_flip_A,*list_GLN_flip_A,*list_HIS_flip_A]
        list_flips = [*list_ASN_flip,*list_GLN_flip,*list_HIS_flip]
        #assemble the system into a single geometry
        cutgeom_A,cutatmlist_A,cutatomtypelist_A,cutresiduelist_A,cutresiduenumbers_A,cutchainidentifiers_A,cutresiduetypes_A,_ = self.assemble_system(residues_A,"A")
        cutgeom_B,cutatmlist_B,cutatomtypelist_B,cutresiduelist_B,cutresiduenumbers_B,cutchainidentifiers_B,cutresiduetypes_B,_ = self.assemble_system(residues_B,"B")
        #save the xyz
        GEOM.Write2XYZ(cutatmlist_A,cutgeom_A,str(xyz_files[0][1]))
        GEOM.Write2XYZ(cutatmlist_B,cutgeom_B,str(xyz_files[0][2]))
        #save the PDB
        GEOM.WriteGemmiPDB(cutgeom_A,cutatmlist_A,cutatomtypelist_A,cutresiduelist_A,cutresiduenumbers_A,cutchainidentifiers_A,cutresiduetypes_A,str(xyz_files[0][1]))
        GEOM.WriteGemmiPDB(cutgeom_B,cutatmlist_B,cutatomtypelist_B,cutresiduelist_B,cutresiduenumbers_B,cutchainidentifiers_B,cutresiduetypes_B,str(xyz_files[0][2]))
        sdf_files_available = False
        if self.config.mmcif_file:
            geometry_fragment_A, atoms_fragment_A, atomtypes_fragment_A, connectivity_A = self.extract_whole_connectivity(residues_A)
            geometry_fragment_B, atoms_fragment_B, atomtypes_fragment_B, connectivity_B = self.extract_whole_connectivity(residues_B)
            GEOM.WriteSDF(atomtypes_fragment_A,atoms_fragment_A,geometry_fragment_A,connectivity_A,[],str(str(xyz_files[0][1])))
            GEOM.WriteSDF(atomtypes_fragment_B,atoms_fragment_B,geometry_fragment_B,connectivity_B,[],str(str(xyz_files[0][2])))
            sdf_files_available = True
        if len(list_flips) > 0:
            #number of possible flips
            nHIS = len(list_HIS_flip)
            nGLN = len(list_GLN_flip)
            nASN = len(list_ASN_flip)
            nflips = (np.power(2,nHIS)*np.power(2,nGLN)*np.power(2,nASN))
            nres_permut = nHIS + nGLN + nASN
            self.logger.info(f"Performing residue flipping.")
            self.logger.info(f"    nASN: {nASN}")
            for asn in list_ASN_flip:
                self.logger.info(f"      {asn}")
            self.logger.info(f"    nGLN: {nGLN}")
            for gln in list_GLN_flip:
                self.logger.info(f"      {gln}")
            self.logger.info(f"    nHIS: {nHIS}")
            for his in list_HIS_flip:
                self.logger.info(f"      {his}")
            self.logger.info(f"Total number of configurations to test: {nflips}")
            import itertools
            combinations = [list(i) for i in itertools.product([0,1],repeat = nres_permut)]
            self.logger.info(f"Combinations:")
            icomb = 0
            for comb in combinations:
                self.logger.info(f"    {comb}")
                #check if this is the zero combination, since that was already done
                sum = 0
                for iel in comb: sum += iel
                if sum != 0:
                    #generate deep copies of the residue data
                    flip_residues_A = copy.deepcopy(residues_A)
                    flip_residues_B = copy.deepcopy(residues_B)
                    #flip the residues
                    for iel in range(len(comb)):
                        if comb[iel] != 0:
                            #modify this residue
                            if list_flips[iel] in list_flips_A:
                                flipped_residue = self._flip_residue(flip_residues_A[list_flips[iel]])
                                flipped_residue.name = flip_residues_A[list_flips[iel]].name + "_flip"
                                flipped_residue.residue_id += "f"
                                flip_residues_A[list_flips[iel]] = flipped_residue
                            else:
                                flipped_residue = self._flip_residue(flip_residues_B[list_flips[iel]])
                                flipped_residue.name = flip_residues_B[list_flips[iel]].name + "_flip"
                                flipped_residue.residue_id += "f"
                                flip_residues_B[list_flips[iel]] = flipped_residue
                    #now assemble the real "structures"
                    cutgeom_A,cutatmlist_A,cutatomtypelist_A,cutresiduelist_A,cutresiduenumbers_A,cutchainidentifiers_A,cutresiduetypes_A,_ = self.assemble_system(flip_residues_A,"A")
                    cutgeom_B,cutatmlist_B,cutatomtypelist_B,cutresiduelist_B,cutresiduenumbers_B,cutchainidentifiers_B,cutresiduetypes_B,_ = self.assemble_system(flip_residues_B,"B")
                    #save the xyz
                    icomb += 1
                    flip_xyz_files = [self.config.geometry_folder / f"{self.config.savefile}_flip_{icomb}_AB",
                                      self.config.geometry_folder / f"{self.config.savefile}_flip_{icomb}_A",
                                      self.config.geometry_folder / f"{self.config.savefile}_flip_{icomb}_B"]
                    #save the xyz
                    GEOM.Write2XYZ(cutatmlist_A,cutgeom_A,str(flip_xyz_files[1]))
                    GEOM.Write2XYZ(cutatmlist_B,cutgeom_B,str(flip_xyz_files[2]))
                    #save the PDB
                    GEOM.WriteGemmiPDB(cutgeom_A,cutatmlist_A,cutatomtypelist_A,cutresiduelist_A,cutresiduenumbers_A,cutchainidentifiers_A,cutresiduetypes_A,str(flip_xyz_files[1]))
                    GEOM.WriteGemmiPDB(cutgeom_B,cutatmlist_B,cutatomtypelist_B,cutresiduelist_B,cutresiduenumbers_B,cutchainidentifiers_B,cutresiduetypes_B,str(flip_xyz_files[2]))
                    xyz_files.append(flip_xyz_files)
        list_of_fragments = []
        for files in xyz_files:
            newfrag = FragmentState()
            newfrag.file_path_xyz = files
            newfrag.file_path_pdb = files
            if sdf_files_available: 
                newfrag.file_path_sdf = files
                newfrag.sdf_files_are_available = True
            newfrag.number_ATOM = [types_list_A[0] + types_list_B[0],types_list_A[0],types_list_B[0]]
            newfrag.number_HETATM = [types_list_A[1] + types_list_B[1],types_list_A[1],types_list_B[1]]
            list_of_fragments.append(newfrag)
        return list_of_fragments
    
    def fragment_analysis_core(self,list_of_fragments) -> None:
        """Function running the core of the fragment analysis"""
        self.logger.info("Fragment Analysis core calculation")
        optimisation_type = "1"
        if self.config.calculation_mode is CalculationMode.OPT:
            optimisation_type = "0"
        elif self.config.calculation_mode is CalculationMode.NOOPT:
            optimisation_type = "-1"
        for idx,fragment in enumerate(list_of_fragments,start = 1):
            filename = str(fragment.file_path_xyz_H[0]).replace(".sdf", "") + ".xyz"
            densfile = (filename.replace(".xyz","_ipdens").replace(str(self.config.geometry_folder) + "/start_xyz",str(self.config.geometry_folder) + "/ip"))
            optgfile_base = (filename.replace(".xyz","_ip").replace(str(self.config.geometry_folder) + "/start_xyz",str(self.config.geometry_folder) + "/ip"))
            ip_dir = Path(optgfile_base).parent
            ip_dir.mkdir(parents = True,exist_ok = True)
            #the file we *want* to exist after RunIP
            expected_xyz = Path(optgfile_base + ".xyz")
            #if NOOPT: create the expected xyz directly
            if optimisation_type == "-1":
                subprocess.run(["cp",filename,str(expected_xyz)],check = False)
            charge = fragment.charge[0]
            #run RunIP inside ip_dir so any default xtb outputs land here.
            #IMPORTANT: make paths absolute before chdir(), otherwise relative paths break.
            filename_abs = os.path.abspath(filename)
            densfile_abs = os.path.abspath(densfile)
            optgfile_base_abs = os.path.abspath(optgfile_base)
            #define expected output geometry path (deterministic)
            expected_xyz = Path(optgfile_base_abs + ".xyz")
            old_cwd = os.getcwd()
            os.chdir(str(ip_dir))
            try:
                out,err = QM.RunIP(filename_abs,
                                str(charge),
                                str(self.config.alpha),
                                str(self.config.kappa),
                                str(self.config.ignh),
                                optimisation_type,
                                str(self.config.Telec),
                                str(self.config.solvation),
                                self.config.solventname,
                                #str(self.config.calculate_density),
                                #densfile_abs,
                                optgfile_base_abs,
                                str(self.config.energy_threshold),
                                str(self.config.gradient_threshold))
                #recenter the whole block to the complex frame
                GEOM.ReCenterAndWrite(str(filename_abs),str(optgfile_base_abs) + ".xyz",str(optgfile_base_abs))
                if self.config.calculate_density:
                    QM.GetDensity(str(optgfile_base_abs) + ".xyz",str(charge),str(self.config.Telec),str(self.config.solvation),self.config.solventname,densfile_abs)
            finally:
                os.chdir(old_cwd)
            #parse energies etc. (keep original behaviour)
            etotal_ip,rmsd_ip,totalpol_ip,atmlist,charges_ip,polarisabilities_ip,converged_ip = QM.ParseOutput(out)
            fragment.energy_ip[0] = etotal_ip
            fragment.RMSD[0] = rmsd_ip
            fragment.total_polarisability[0] = totalpol_ip
            fragment.atom_list = atmlist
            fragment.atomic_charges = charges_ip
            fragment.atomic_polarisabilities = polarisabilities_ip
            fragment.converged_ip[0] = converged_ip
            #ensure the expected geometry exists; recover common xtb default names
            if not expected_xyz.exists():
                #look for typical xtb outputs in ip_dir (since we forced cwd there)
                base_name = Path(optgfile_base_abs).name
                candidates = [ip_dir / "xtbopt.xyz",
                              ip_dir / "xtbopt.coord",
                              ip_dir / f"{base_name}.xyz",
                              ip_dir / base_name,
                              Path(optgfile_base_abs),            # absolute base (no .xyz)
                              Path(optgfile_base_abs + ".xyz")]   # absolute expected (redundant but safe)
                found = next((c for c in candidates if c.exists()),None)
                if found is not None:
                    #prefer renaming to expected_xyz so downstream is deterministic
                    try:
                        shutil.move(str(found),str(expected_xyz))
                    except Exception:
                        shutil.copy2(str(found),str(expected_xyz))
                else:
                    #clean actionable error: include stderr excerpt + folder listing
                    err_txt = ""
                    try:
                        if isinstance(err,(bytes,bytearray)):
                            err_txt = err.decode("utf-8",errors = "replace")
                        else:
                            err_txt = str(err)
                    except Exception:
                        err_txt = "<unable to decode stderr>"
                    listing = []
                    try:
                        listing = sorted([p.name for p in ip_dir.iterdir()])[:200]
                    except Exception:
                        listing = ["<unable to list ip_dir>"]
                    raise FragmentAnalyzerError("RunIP finished but ip geometry not found.\n"
                                                f"Expected: {expected_xyz}\n"
                                                f"ip_dir: {ip_dir}\n"
                                                f"ip_dir listing (first 200): {listing}\n"
                                                f"RunIP stderr (first 400 chars): {err_txt[:400]}")
            #store REAL path for downstream ReadXYZ()
            fragment.file_path_ip = [expected_xyz]
            self.logger.info(f"    Fragment Analysis on {filename} done")
            if not converged_ip:
                self.logger.warning("in-pocket not converged")
        
    def select_fragment_states(self, list_smiles, fragment_id) -> List[int]:
        self.logger.info(f"WARNING: several chemical states identified for fragment {fragment_id}:")
        for smile in list_smiles:
            self.logger.info(smile)
        codes_2_do = input("please select 'first', 'all', or comma-separated input list: ")
        codes_to_do = []
        if codes_2_do.lower() == "all": 
            for idx in range(len(list_smiles)): codes_to_do.append(idx)
        elif codes_2_do.lower() == "first": 
            codes_to_do = [0]
        else: 
            completelist = codes_2_do.split(",")
            for idl in completelist:
                if "-" in idl:
                    newsplit = idl.split("-")
                    for idx in range(int(newsplit[0]),int(newsplit[1]) + 1):
                        index = idx - 1
                        codes_to_do.append(index)
                else:
                    index = int(idl) - 1
                    codes_to_do.append(index)
        return codes_to_do
    
    def tautomer_search(self, files_old, charges_old, smiles_old, geom_base, identifier) -> Tuple[List[str],List[int],List[str]]:
        """Performing tautomer search and updating the structure list"""
        self.logger.info(f"Searching for tautomers in fragment {identifier}")
        new_files = []
        new_charges = []
        new_smiles = []
        for iprot in range(len(files_old)):
            print(str(geom_base),str(files_old[iprot]))
            tautomers, smiles = TAUTORO.process_tautomers(str(files_old[iprot]),str(geom_base),self.config.max_tautomers)
            if not tautomers: 
                new_files.append(files_old[iprot])
                new_charges.append(charges_old[iprot])
                new_smiles.append(smiles_old[iprot])
            else: 
                for itaut in range(len(tautomers)):
                    find_files = files_old[iprot].rfind("/")
                    find_ext = files_old[iprot].rfind(".")
                    str_to_remove = files_old[iprot][find_files + 1:find_ext] + "_tautomers/"
                    tautomer_moved = tautomers[itaut].replace(str_to_remove,"")
                    subprocess.run(["mv",tautomers[itaut],tautomer_moved])
                    new_files.append(tautomer_moved)
                    new_charges.append(charges_old[iprot])
                    new_smiles.append(smiles[itaut])
        return new_files,new_charges,new_smiles
    
    def protonate_structures(self,list_of_fragments) -> List[FragmentState]:
        """Protonate all blocks of structures and get their charges"""
        self.logger.info("Protonating all substructures")
        import copy
        new_list_of_fragments = []
        if self.config.protonation_mode is ProtonationMode.RDKIT: 
            PROTORO.print_header_protoro()
            for fragment in list_of_fragments:
                fileAB = str(fragment.file_path_xyz[0])
                fileA = str(fragment.file_path_xyz[1])
                fileB = str(fragment.file_path_xyz[2])
                chargeA = 0
                chargeB = 0
                #protonate each species isolatedly, to make extraction easier later
                ext = ".pdb"
                if fragment.sdf_files_are_available: ext = ".sdf"
                A_files,A_charges,A_smiles = PROTORO.RDKitProtonation(fileA + ext,self.config.pH,"lowres",pH_spread = self.config.pH_spread)
                B_files,B_charges,B_smiles = PROTORO.RDKitProtonation(fileB + ext,self.config.pH,"lowres",pH_spread = self.config.pH_spread)
                if self.config.tautomerise:
                    self.logger.info("Searching for tautomers")
                    geom_base = self.config.geometry_folder
                    if (fragment.number_HETATM[1] > 0) and (fragment.number_ATOM[1] < self.config.max_natom_tautomers):
                        A_files,A_charges,A_smiles = self.tautomer_search(A_files,A_charges,A_smiles,geom_base,"A")
                    if (fragment.number_HETATM[2] > 0) and (fragment.number_ATOM[2] < self.config.max_natom_tautomers):
                        B_files,B_charges,B_smiles = self.tautomer_search(B_files,B_charges,B_smiles,geom_base,"B")
                #if many protonation states are selected, let the user select some
                codes_to_do_A = [0]
                if len(A_smiles) > 1:
                    codes_to_do_A = self.select_fragment_states(A_smiles,"A")
                codes_to_do_B = [0]
                if len(B_smiles) > 1:
                    codes_to_do_B = self.select_fragment_states(B_smiles,"B")
                for istate_A in range(len(codes_to_do_A)):
                    idx_A = codes_to_do_A[istate_A]
                    for jstate_B in range(len(codes_to_do_B)):
                        idx_B = codes_to_do_B[jstate_B]
                        temp_fragment = copy.deepcopy(fragment)
                        #assign sdf files
                        filename_A = A_files[idx_A].replace(".sdf","")
                        filename_B = B_files[idx_B].replace(".sdf","")
                        temp_fragment.file_path_sdf[1] = filename_A
                        temp_fragment.file_path_sdf[2] = filename_B
                        #get the xyz
                        atoms_A,geometry_A = GEOM.ReadSDF(filename_A)
                        GEOM.Write2XYZ(atoms_A,geometry_A,filename_A)
                        atoms_B,geometry_B = GEOM.ReadSDF(filename_B)
                        GEOM.Write2XYZ(atoms_B,geometry_B,filename_B)
                        B_file_simplified = B_files[idx_B].replace(str(self.config.geometry_folder),"").replace("/","_")
                        AB_file = A_files[idx_A].replace(".sdf",B_file_simplified).replace(str(self.config.geometry_folder),str(self.config.geometry_folder) + "/start_xyz").replace(".sdf","")
                        A_file = filename_A.replace(str(self.config.geometry_folder),str(self.config.geometry_folder) + "/start_xyz")
                        B_file = filename_B.replace(str(self.config.geometry_folder),str(self.config.geometry_folder) + "/start_xyz")
                        temp_fragment.file_path_xyz_H = [AB_file,A_file,B_file]
                        #get charges
                        temp_fragment.charge = [A_charges[idx_A] + B_charges[idx_B],A_charges[idx_A],B_charges[idx_B]]
                        #build the pair
                        atoms_AB = [*atoms_A,*atoms_B]
                        temp_fragment.atom_list = atoms_AB
                        geometry_AB = [*geometry_A,*geometry_B]
                        #write the geometry
                        GEOM.Write2XYZ(atoms_AB,geometry_AB,AB_file)
                        #store fragment
                        new_list_of_fragments.append(temp_fragment)
        elif self.config.protonation_mode is ProtonationMode.OBABEL: 
            for fragment in list_of_fragments:
                fileAB = str(fragment.file_path_xyz[0])
                fileA = str(fragment.file_path_xyz[1])
                fileB = str(fragment.file_path_xyz[2])
                #protonate each species isolatedly, to make extraction easier later
                chargeA = GEOM.ObabelProtonation(fileA + ".xyz",self.config.pH)
                GEOM.MOL2toXYZ(fileA + "H.mol2","dummy")
                chargeB = GEOM.ObabelProtonation(fileB + ".xyz",self.config.pH)
                GEOM.MOL2toXYZ(fileB + "H.mol2","dummy")
                charge = (input("Please check the consistency of the protonation in the files " + fileA + "H.xyz and " + fileB + "H.xyz enter the total charges separated by comma.\n (Predicted charges " + str(chargeA) + "," + str(chargeB) + "): ")).split(",")
                #store charges
                chargeA = int(charge[0])
                chargeB = int(charge[1])
                fragment.charge = [chargeA + chargeB,chargeA,chargeB]
                #new file names
                fileAB += "H"
                fileA += "H"
                fileB += "H"
                fileAB_new = fileAB.replace(str(self.config.geometry_folder),str(self.config.geometry_folder) + "/start_xyz")
                fileA_new = fileA.replace(str(self.config.geometry_folder),str(self.config.geometry_folder) + "/start_xyz")
                fileB_new = fileB.replace(str(self.config.geometry_folder),str(self.config.geometry_folder) + "/start_xyz")
                fragment.file_path_xyz_H = [Path(fileAB_new),Path(fileA_new),Path(fileB_new)]
                fragment.file_path_pdb_H = [Path(fileAB_new),Path(fileA_new),Path(fileB_new)]
                fragment.file_path_xyz = [Path(fileAB_new),Path(fileA_new),Path(fileB_new)]
                fragment.file_path_pdb = [Path(fileAB_new),Path(fileA_new),Path(fileB_new)]
                #assemble the pair
                atomsA,geometryA = GEOM.ReadXYZ(fileA)
                atomsB,geometryB = GEOM.ReadXYZ(fileB)
                atomsAB = [*atomsA,*atomsB]
                fragment.atom_list = atomsAB
                geometryAB = [*geometryA,*geometryB]
                #write geometries
                GEOM.Write2XYZ(atomsAB,geometryAB,fileAB)
                temp_fragment = copy.deepcopy(fragment)
                new_list_of_fragments.append(temp_fragment)
        #elif self.config.protonation_mode is ProtonationMode.VSEPR: 
                ####FMCM: this is for VSEPR
                ####elif (self.config.protonation_mode is ProtonationMode.VSEPR) or (self.config.protonation_mode is ProtonationMode.CCD):
                ####    print("in1")
                ####    if (self.config.protonation_mode is ProtonationMode.VSEPR):
                ####        print("in2")
                ####        chargeA = VSEPRProtonation(fileA + ".xyz",fileA + "H")
                ####        chargeB = VSEPRProtonation(fileB + ".xyz",fileB + "H")
                ####    else:
                ####        if useccd == "A":
                ####            if self.config.protonation_mode is ProtonationMode.CCD:
                ####                chargeA = CCDProtonation(fileA + ".pdb",fileA + ".pdb",ccdfile,fileA + "H")
                ####                subprocess.run(["rm",fileA + "H.sdf"])
                ####            elif self.config.protonation_mode is ProtonationMode.VSEPR:
                ####                chargeA = VSEPRProtonation(fileA + ".pdb",fileA + "H")
                ####            chargeB = ProteinProtonation(fileB + ".pdb",fileB + "H")
                ####        elif useccd == "B":
                ####            if self.config.protonation_mode is ProtonationMode.CCD:
                ####                chargeB = CCDProtonation(fileB + ".pdb",fileB + ".pdb",ccdfile,fileB + "H")
                ####                subprocess.run(["rm",fileB + "H.sdf"])
                ####            elif self.config.protonation_mode is ProtonationMode.VSEPR:
                ####                chargeB = VSEPRProtonation(fileB + ".pdb",fileB + "H")
                ####            chargeA = ProteinProtonation(fileA + ".pdb",fileA + "H")
        #move files
        list_of_files = os.listdir(str(self.config.geometry_folder))
        for file in list_of_files:
            if file.endswith(".xyz"):
                subprocess.run(["mv",str(self.config.geometry_folder) + "/" + file,str(self.config.geometry_folder) + "/start_xyz/" + file])
            elif file.endswith(".pdb"):
                subprocess.run(["mv",str(self.config.geometry_folder) + "/" + file,str(self.config.geometry_folder) + "/start_pdb/" + file])
            elif file.endswith(".sdf"):
                subprocess.run(["mv",str(self.config.geometry_folder) + "/" + file,str(self.config.geometry_folder) + "/start_sdf/" + file])
        return new_list_of_fragments
    
    def reassemble_system_from_previous_run(self, file_A, file_B) -> List[FragmentState]:
        """Function to fetch fragment system from a previous run (geometries)"""
        self.logger.info("Assembling geometry systems into fragments")
        atoms_A,geometry_A = GEOM.ReadGeometry(file_A)
        atoms_B,geometry_B = GEOM.ReadGeometry(file_B)
        xyz_files = [str(self.config.geometry_folder) + "/start_xyz/" + f"{self.config.savefile}_AB",
                     str(self.config.geometry_folder) + "/start_xyz/" + f"{self.config.savefile}_A",
                     str(self.config.geometry_folder) + "/start_xyz/" + f"{self.config.savefile}_B"]
        complex_atoms = [*atoms_A,*atoms_B]
        complex_geometry = [*geometry_A,*geometry_B]
        GEOM.Write2XYZ(complex_atoms,complex_geometry,str(xyz_files[0]))
        GEOM.Write2XYZ(atoms_A,geometry_A,str(xyz_files[1]))
        GEOM.Write2XYZ(atoms_B,geometry_B,str(xyz_files[2]))
        newfrag = FragmentState()
        newfrag.file_path_xyz = xyz_files
        newfrag.file_path_pdb = xyz_files
        newfrag.file_path_xyz_H = xyz_files
        newfrag.file_path_pdb_H = xyz_files
        newfrag.sdf_files_are_available = False
        chargeA = self.config.chargeA
        chargeB = self.config.chargeB
        newfrag.charge = [chargeA + chargeB,chargeA,chargeB]
        newfrag.atom_list = complex_atoms
        return [newfrag]
    
    def reassemble_system_from_IPA_run(self) -> List[FragmentState]:
        """Function to reassemble fragment system from a previous IPA run"""
        self.logger.info("Assembling IPA system into fragments")
        #reassembling ligand file
        if not self.path_to_IPA_geometries.endswith("/"): self.path_to_IPA_geometries += "/"
        ligand_file = self.path_to_IPA_geometries + self.ligand_name_from_IPA
        #new CSV already stores basenames like "..._ip".
        #old CSV used "set: ip" to request adding "_ip".
        if getattr(self, "set_from_IPA", "").lower() == "ip" and (not ligand_file.endswith("_ip")):
            ligand_file += "_ip"
        ligand_atoms,ligand_geometry = GEOM.ReadXYZ(ligand_file)
        #reassembling pocket system
        status,pocket_atoms,pocket_geometry,pocket_charge = GEOM.RebuildPocketFromResidues(self.residue_states_from_IPA,self.path_to_IPA_geometries)
        if status == 0:
            self.logger.error(f"Residue repeated in the supplied residue list. This means ambiguous residue protonation state declared.")
            sys.exit(1)
        xyz_files = [str(self.config.geometry_folder) + "/start_xyz/" + f"{self.config.savefile}_AB",
                     str(self.config.geometry_folder) + "/start_xyz/" + f"{self.config.savefile}_A",
                     str(self.config.geometry_folder) + "/start_xyz/" + f"{self.config.savefile}_B"]
        GEOM.Write2XYZ(ligand_atoms,ligand_geometry,str(xyz_files[1]))
        GEOM.Write2XYZ(pocket_atoms,pocket_geometry,str(xyz_files[2]))
        complex_atoms = [*ligand_atoms,*pocket_atoms]
        complex_geometry = [*ligand_geometry,*pocket_geometry]
        GEOM.Write2XYZ(complex_atoms,complex_geometry,str(xyz_files[0]))
        newfrag = FragmentState()
        newfrag.file_path_xyz = xyz_files
        newfrag.file_path_pdb = xyz_files
        newfrag.file_path_xyz_H = xyz_files
        newfrag.file_path_pdb_H = xyz_files
        newfrag.sdf_files_are_available = False
        chargeA = self.ligand_charge_from_IPA
        chargeB = pocket_charge
        newfrag.charge = [chargeA + chargeB,chargeA,chargeB]
        newfrag.atom_list = complex_atoms
        return [newfrag]

    def process_results(self, list_of_fragments) -> None:
        """Function to process all the results"""
        self.logger.info("Processing results")
        if self.config.calculate_Ebind:
            index = 1
            for fragment in list_of_fragments:
                self.logger.info(f"Processing fragment {index}")
                #get the ip-optimised pair
                AB_ip_opt = str(fragment.file_path_ip[0]).replace(".xyz","")
                atomsAB,geometryAB = GEOM.ReadXYZ(AB_ip_opt)
                #get the pre-optimised monomers
                fileA = str(fragment.file_path_xyz_H[1])
                fileB = str(fragment.file_path_xyz_H[2])
                #create files
                A_ip = str(fileA).replace(str(self.config.geometry_folder) + "/start_xyz",str(self.config.geometry_folder) + "/ip") + "_ip"
                B_ip = str(fileB).replace(str(self.config.geometry_folder) + "/start_xyz",str(self.config.geometry_folder) + "/ip") + "_ip"
                #read them
                atomsA_bad,_ = GEOM.ReadXYZ(fileA)
                atomsB_bad,_ = GEOM.ReadXYZ(fileB)
                #prepare to store new geometries
                atomsA = []
                geometryA = []
                atomsB = []
                geometryB = []
                for idatm in range(len(atomsA_bad)):
                    atomsA.append(atomsAB[idatm])
                    geometryA.append(geometryAB[idatm])
                for idbtm in range(len(atomsB_bad)):
                    atomsB.append(atomsAB[len(atomsA_bad) + idbtm])
                    geometryB.append(geometryAB[len(atomsA_bad) + idbtm])
                GEOM.Write2XYZ(atomsA,geometryA,A_ip)
                GEOM.Write2XYZ(atomsB,geometryB,B_ip)
                fragment.file_path_ebind = [AB_ip_opt,A_ip,B_ip]
                optgfile = "dummy"
                A_dens = A_ip + "_dens.xyz"
                B_dens = B_ip + "_dens.xyz"
                chargeA = fragment.charge[1]
                chargeB = fragment.charge[2]
                chargeAB = fragment.charge[0]
                out,_ = QM.RunIP(A_ip + ".xyz",str(chargeA),
                              str(self.config.alpha),
                              str(self.config.kappa),
                              str(self.config.ignh),
                              "-1",
                              str(self.config.Telec),
                              str(self.config.solvation),
                              self.config.solventname,
                              #str(self.config.calculate_density),
                              #A_dens,
                              optgfile,
                              str(self.config.energy_threshold),
                              str(self.config.gradient_threshold))
                etotal_A,_,totalpol_A,_,_,_,_ = QM.ParseOutput(out)
                out,_ = QM.RunIP(B_ip + ".xyz",
                              str(chargeB),
                              str(self.config.alpha),
                              str(self.config.kappa),
                              str(self.config.ignh),
                              "-1",
                              str(self.config.Telec),
                              str(self.config.solvation),
                              self.config.solventname,
                              #str(self.config.calculate_density),
                              #B_dens,
                              optgfile,
                              str(self.config.energy_threshold),
                              str(self.config.gradient_threshold))
                if self.config.calculate_density:
                    QM.GetDensity(A_ip + ".xyz",str(chargeA),str(self.config.Telec),str(self.config.solvation),self.config.solventname,A_dens)
                    QM.GetDensity(B_ip + ".xyz",str(chargeB),str(self.config.Telec),str(self.config.solvation),self.config.solventname,B_dens)
                etotal_B,_,totalpol_B,_,_,_,_ = QM.ParseOutput(out)
                fragment.energy_ip[1] = etotal_A
                fragment.energy_ip[2] = etotal_B
                fragment.total_polarisability[1] = totalpol_A
                fragment.total_polarisability[2] = totalpol_B
                self.logger.info(f"Monomers for fragment {index} processed")
                if self.config.calculate_EDDA:
                    self.logger.info(f"Running EDDA for fragment {index}")
                    out,err = QM.RunEDDA(AB_ip_opt + ".xyz",str(chargeAB),A_ip + ".xyz",str(chargeA),B_ip + ".xyz",str(chargeB))
                    RepulsionE,InteractionE,LocalRepulsionA,LocalRepulsionB,LocalInteractionA,LocalInteractionB = QM.ParseEDDAOutput(out,atomsA,atomsB)
                    fragment.repulsion_energy = RepulsionE
                    fragment.EDDA_binding_energy = InteractionE
                    heavyatm_list_A,heavyatm_rep_A,heavyatm_int_A,shortcontact_A = SCORE.ShortContactScore(atomsA,LocalRepulsionA,LocalInteractionA,geometryA)
                    heavyatm_list_B,heavyatm_rep_B,heavyatm_int_B,shortcontact_B = SCORE.ShortContactScore(atomsB,LocalRepulsionB,LocalInteractionB,geometryB)
                    fragment.heavy_atom_list = [*heavyatm_list_A,*heavyatm_list_B]
                    fragment.repulsion_contacts = [*heavyatm_rep_A,*heavyatm_rep_B]
                    fragment.interaction_contacts = [*heavyatm_int_A,*heavyatm_int_B]
                    fragment.short_contact_indices = [*shortcontact_A,*shortcontact_B]
                    self.logger.info(f"Running EDDA for fragment {index} finished")
                index += 1
        return None

    def write_results(self, list_of_fragments) -> None:
        """Function to write the results"""
        self.logger.info("Writing results")
        writing_file = str(self.config.summaryfile)
        if not writing_file.endswith(".csv"): writing_file += ".csv"
        path_to_summary_file = str(self.config.output_file_folder)
        if not path_to_summary_file.endswith("/"): path_to_summary_file += "/"
        wfile = open(path_to_summary_file + writing_file,"w")
        energy_minimum = float("inf")
        for fragment in list_of_fragments:
            if fragment.energy_ip[0] < energy_minimum: energy_minimum = fragment.energy_ip[0]
        wfile.write("EnergyAB(Eh);ErelativeAB(kcal/mol);Ebind(kcal/mol);Erep(kcal/mol);RMSD(A);TotalPolarisabilityAB(a0^3);TotalPolarisabilityA(a0^3);TotalPolarisabilityB(a0^3);Converged\n")
        for fragment in list_of_fragments:
            energy_relative = (fragment.energy_ip[0] - energy_minimum)*Eh_to_kcalmol
            energy_binding = 0.0
            polarisability_A = 0.0
            polarisability_B = 0.0
            energy_repulsion = 0.0
            if self.config.calculate_Ebind: 
                energy_binding = (fragment.energy_ip[0] - fragment.energy_ip[1] - fragment.energy_ip[2])*Eh_to_kcalmol
                polarisability_AB = fragment.total_polarisability[0]
                polarisability_A = fragment.total_polarisability[1]
                polarisability_B = fragment.total_polarisability[2]
            if self.config.calculate_EDDA: energy_repulsion = fragment.repulsion_energy
            energy_relative_string = "               " + FORMAT.ReformatingString(f'{energy_relative:.5}',6) + ";"
            energy_binding_string = ""
            if energy_binding >= 0.0: energy_binding_string = " "
            energy_binding_string += "        " + FORMAT.ReformatingString(f'{energy_binding:.5}',6) + ";"
            energy_repulsion_string = "        " + FORMAT.ReformatingString(f'{energy_repulsion:.5}',6) + ";"
            TotalPolAB_string = "                    " + FORMAT.ReformatingString(f'{polarisability_AB:.6}',7) + ";"
            TotalPolA_string = "                   " + FORMAT.ReformatingString(f'{polarisability_A:.6}',7) + ";"
            TotalPolB_string = "                   " + FORMAT.ReformatingString(f'{polarisability_B:.6}',7) + ";"
            wfile.write(str(fragment.energy_ip[0]) + ";" + energy_relative_string + energy_binding_string + energy_repulsion_string + f"{fragment.RMSD[0]:.3};" + TotalPolAB_string + TotalPolA_string + TotalPolB_string + str(fragment.converged_ip[0]) + "\n")
            wfile2 = open(str(fragment.file_path_ip[0]).replace(".xyz","").replace(str(self.config.geometry_folder) + "/ip",str(self.config.output_file_folder)) + ".csv","w")
            wfile2.write("AtomicNumber;Charge;Polarisability(a0^3);lrep(kcal/mol);lint(kcal/mol);sc\n")
            localrepulsion_AB = fragment.repulsion_contacts
            localinteraction_AB = fragment.interaction_contacts
            shortcontact_AB = fragment.short_contact_indices
            iheavy = 0
            for idatm in range(len(fragment.atom_list)):
                lrep = 0.0
                lint = 0.0
                sc = 0.0
                if self.config.calculate_EDDA:
                    if (fragment.atom_list[idatm] != 1) and (fragment.atom_list[idatm] != "H"):
                        lrep = localrepulsion_AB[iheavy]*Eh_to_kcalmol
                        lint = localinteraction_AB[iheavy]*Eh_to_kcalmol
                        sc = math.fabs(shortcontact_AB[iheavy])
                        iheavy += 1
                lrep_string = "        " + FORMAT.ReformatingString(f'{lrep:.5}',6) + ";"
                lint_string = "        " + FORMAT.ReformatingString(f'{lint:.5}',6) + ";"
                sc_string = "        " + FORMAT.ReformatingString(f'{sc:.5}',6)
                atomnumber_string = ""
                if fragment.atom_list[idatm] < 10: atomnumber_string = " "
                atomnumber_string += str(fragment.atom_list[idatm]) + ";"
                atomiccharge_string = "        " + FORMAT.ReformatingString(f'{fragment.atomic_charges[idatm]:.3}',4) + ";"
                atomicpol_string = "        " + FORMAT.ReformatingString(f'{fragment.atomic_polarisabilities[idatm]:.3}',4) + ";"
                wfile2.write(atomnumber_string + atomiccharge_string + atomicpol_string + lrep_string + lint_string + sc_string + "\n")
            wfile2.close()
        wfile.close()
        return None

    def fragment_analysis(self) -> None:
        """Run the in-pocket fragment analysis."""
        self.logger.info("Starting Fragment Analysis")
        try:
            #validate inputs
            self.validate_inputs()
            #read structure
            residues,non_prot_data = self.read_structure()
            #get fragments and connectors
            active_residues_A,connectorsA,active_residues_B,connectorsB = self.get_fragments_from_residues(residues,non_prot_data)
            #get the fragment geometries and protonate
            if self.run_from_previous_IPA:
                list_of_fragments = self.reassemble_system_from_IPA_run()
            elif self.run_from_geometry:
                list_of_fragments = self.reassemble_system_from_previous_run(active_residues_A,active_residues_B)
            else: 
                list_of_fragments = self.generate_fragment_geometries(residues,active_residues_A,connectorsA,active_residues_B,connectorsB)
                list_of_fragments = self.protonate_structures(list_of_fragments)
            #run Fragment Analysis
            #input_ = input("continue")
            self.fragment_analysis_core(list_of_fragments)
            #process data
            self.process_results(list_of_fragments)
            #write results
            self.write_results(list_of_fragments)
        except FragmentAnalyzerError as e:
            self.logger.error(f"Analysis failed: {e}")
            raise
        except Exception as e:
            self.logger.error(f"Unexpected error: {e}")
            raise FragmentAnalyzerError(f"Unexpected error: {e}")
        return None

def main():
    """Main entry point for the script."""
    #getting arguments
    args = parse_arguments()
    #set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    else:
        logging.getLogger().setLevel(logging.INFO)
    try:
        #convert arguments to configuration
        config = parse_arguments_to_config(args)
        #create analyzer
        analyzer = FragmentAnalyzer(config)
        #run analysis
        analyzer.fragment_analysis()
    except FragmentAnalyzerError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("\nAnalysis interrupted by user", file=sys.stderr)
        sys.exit(130)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
