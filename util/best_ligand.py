# ------------------------------------------------------------
# Utilities to post-process In-Pocket results and generate
# compact "best ligand" outputs:
# - best_interaction_energies.csv
# - best_lig.sdf  (SDF with IP-optimized coordinates + metadata)
# ------------------------------------------------------------

from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import List,Optional,Sequence,Tuple,Callable
import logging
import shutil
import re
import util.UnitConversion as UC

try:
    import numpy as _np         #optional
    _INT_TYPES = (int,_np.integer)
except Exception:
    _INT_TYPES = (int,)

#same conversion factor as in InPocketAnalyzer
EH2KCAL = UC.Hartree2Kcalmol

@dataclass
class BestLigandSummary:
    """Small container for the best ligand information."""
    index: int                         #0-based index in ligand_states
    ipa_binding_energy: float          #kcal/mol
    total_interaction_energy: float    #kcal/mol
    residues: List[str]                #e.g. ["A51", "B253", ...]

#charge-state helpers (mirror in_pocket.py)
def _format_charge_state(charge: int) -> str:
    """Format an integer charge as a stable residue-state string ('0', '+1', '-2', ...)."""
    try:
        q = int(charge)
    except Exception:
        q = 0
    return "0" if q == 0 else f"{q:+d}"

def _parse_charge_state(state) -> int:
    """Parse residue_state strings from ('+','--') or ('+1','-2','0','1') formats."""
    if state is None:
        return 0
    if isinstance(state,_INT_TYPES):
        return int(state)
    s = str(state).strip()
    if not s:
        return 0
    if re.fullmatch(r"\++",s):
        return len(s)
    if re.fullmatch(r"-+",s):
        return -len(s)
    try:
        return int(s)
    except Exception:
        u = s.upper()
        if u in ("POS","PLUS","P"):
            return 1
        if u in ("NEG","MINUS","M"):
            return -1
        return 0

def _parse_residue_id(residue_name: str, residue_state: str) -> Tuple[str, str, str, bool]:
    """Parse the residue identifier used in BindingResult.residue_name.
       Examples:
           "TYR_B253"          -> ("TYR", "B", "253", False)
           "TYR_B253_1"        -> ("TYR", "B", "253", False)
           "ASN_B368_flip"     -> ("ASN", "B", "368", True)
           "ASN_B368_1_flip"   -> ("ASN", "B", "368", True)"""
    #remove explicit "_flip" annotation for parsing
    name = residue_name.replace("_flip","")
    parts = name.split("_")
    base_name = parts[0]                              #e.g., "TYR"
    chain_res = parts[1] if len(parts) > 1 else "A1"  #e.g., "B253" or "A1"
    #ignore any trailing protonation index in chain_res part, if present
    #(e.g. "B253_1" -> "B253")
    chain_res = chain_res.split("_")[0]
    #extract chain and residue number
    chain_id = "".join(c for c in chain_res if not c.isdigit()).replace("f","")
    res_num = "".join(c for c in chain_res if c.isdigit())
    #flipped?
    is_flipped = residue_name.endswith("f") or ("_flip" in residue_name)
    #proton_state = residue_state
    return base_name,chain_id,res_num,is_flipped

def _picked_to_list_res_sorted(picked_res: Sequence[str]) -> List[str]:
    """Convert picked residue entries like 'MET_A69_0' to LIST_RES entries like 'A69',
       sorted lexicographically (string sort), matching in_pocket conventions."""
    list_res_simple: List[str] = []
    for r in picked_res or []:
        parts = str(r).split("_")
        if len(parts) >= 2:
            chain_res = parts[1].replace("f","")
            list_res_simple.append(chain_res)
    return sorted(list_res_simple)

def _compute_clash_flag(result) -> bool:
    """Best-effort reproduction of in_pocket's EDDA clash boolean for one residue result.
       Returns True if any residue-side short_contact_index <= -5.0."""
    cl = getattr(result,"clash",None)
    if cl is None:
        return False
    scores = list(getattr(cl,"short_contact_indices",[]) or [])
    if not scores:
        return False
    labels = list(getattr(cl,"atom_labels",[]) or [])
    rid = str(getattr(cl,"residue_id","") or "")
    res_scores: List[float] = []
    used_label_split = False
    #prefer label-based selection if labels look residue-prefixed
    if rid and labels and len(labels) == len(scores):
        for lab,v in zip(labels,scores):
            if isinstance(lab,str) and lab.startswith(rid + ":"):
                try:
                    res_scores.append(float(v))
                    used_label_split = True
                except Exception:
                    pass
    #fallback: split_index separates ligand atoms (first) from residue atoms (second)
    if not used_label_split:
        try:
            split = int(getattr(cl,"split_index",0) or 0)
        except Exception:
            split = 0
        split = min(max(split,0),len(scores))
        for v in scores[split:]:
            try:
                res_scores.append(float(v))
            except Exception:
                pass
    return any(v <= -5.0 for v in res_scores)

def _compute_ligand_interactions_and_ipa(
        ligand_states: Sequence,
        binding_results: Sequence[Sequence],
        logger: Optional[logging.Logger] = None,
        compute_total_interaction_and_list_res: Optional[Callable[[Sequence], Tuple[float, List[str]]]] = None
    ) -> Tuple[List[float], List[float], List[List[str]]]:
    """Reproduce the IPA binding energy logic from InPocketAnalyzer.write_results,
       but entirely in memory.
       Returns:
           ligand_int_energies:   total IPA interaction energy per ligand (kcal/mol)
           ipa_binding_energies:  IPA binding energy per ligand (kcal/mol)
           residues_per_ligand:   list of residue lists (e.g. ["A51","B253",...]) per ligand"""
    if not ligand_states:
        raise ValueError("No ligand states provided to best_ligand module")
    #global minimum chemical state energy (Hartree), as in write_results
    valid_energies = [s.energy_ip for s in ligand_states if s.energy_ip is not None]
    if not valid_energies:
        raise ValueError("No valid ligand energy_ip values found")
    ligand_min_energy = min(valid_energies)
    ligand_int_energies: List[float] = []
    ipa_binding_energies: List[float] = []
    residues_per_ligand: List[List[str]] = []
    for lig_idx,(ligand,results) in enumerate(zip(ligand_states,binding_results)):
        #prefer using in_pocket's own implementation (passed as callback) to guarantee consistency
        if compute_total_interaction_and_list_res is not None:
            total_interaction,picked_res = compute_total_interaction_and_list_res(results)
            list_res_sorted = _picked_to_list_res_sorted(picked_res)
        else:
            #fallback implementation (best-effort match to current in_pocket rules; no --no-cap disambiguation)
            total_interaction = 0.0
            residues_done: List[str] = []
            picked_res: List[str] = []
            for result in results:
                base_name,chain_id,res_num,is_flipped = _parse_residue_id(result.residue_name,result.residue_state)
                if is_flipped:
                    continue
                q = _parse_charge_state(getattr(result,"residue_state",None))
                take = False
                b = str(base_name or "").upper()
                if b in ("ARG","LYS"):
                    take = (q == 1)
                elif b == "TYR":
                    take = (q == 0)
                elif b in ("ASP","GLU"):
                    take = (q == -1)
                elif b == "HIS":
                    key = f"{b}{res_num}"
                    if (q == 0) and (key not in residues_done):
                        take = True
                        residues_done.append(key)
                else:
                    take = True
                if take:
                    total_interaction += float(getattr(result,"binding_energy",0.0) or 0.0)
                    picked_res.append(f"{b}_{chain_id}{res_num}_{_format_charge_state(q)}")
            list_res_sorted = _picked_to_list_res_sorted(picked_res)
        ligand_int_energies.append(total_interaction)
        residues_per_ligand.append(list_res_sorted)
        #IPA binding energy = Deformation + (E_ip - E_min + acid_base_correction)*Eh2kcal + IPA_interaction
        deformation = ligand.deformation_energy if ligand.deformation_energy is not None else 0.0
        acid_base = getattr(ligand,"acid_base_correction",None)
        acid_base = float(acid_base) if acid_base is not None else 0.0
        if ligand.energy_ip is None:
            if logger:
                logger.warning(f"Ligand {lig_idx+1} has no energy_ip; setting chemical state term to 0")
            chem_state = 0.0
        else:
            #chem_state = (ligand.energy_ip - ligand_min_energy) * EH2KCAL 260223
            chem_state = (ligand.energy_ip - ligand_min_energy + acid_base) * EH2KCAL
        ipa_binding = deformation + chem_state + total_interaction
        ipa_binding_energies.append(ipa_binding)
        if logger:
            logger.info(f"[best_ligand] Ligand {lig_idx+1}: "
                        f"Deformation = {deformation:.2f} kcal/mol, "
                        f"ChemState = {chem_state:.2f} kcal/mol, "
                        f"IPA_interaction = {total_interaction:.2f} kcal/mol, "
                        f"IPA_binding = {ipa_binding:.2f} kcal/mol")
    return ligand_int_energies,ipa_binding_energies,residues_per_ligand

def _write_best_interaction_csv(output_dir: Path,ligand,results: Sequence,filename: str = "best_interaction_energies.csv") -> None:
    """Write best_interaction_energies.csv containing only the best ligand block,
       formatted like interaction_energies.csv."""
    out_path = output_dir / filename
    with out_path.open("w") as f:
        #first line: use SMILES if available, otherwise a generic label
        if getattr(ligand,"smiles",""):
            f.write(f"{ligand.smiles}:\n")
        else:
            f.write("best ligand state:\n")
        #match interaction_energies.csv header conventions (charge + optional clash column)
        edda = any(getattr(r,"clash",None) is not None for r in (results or []))
        if edda:
            f.write("res;chain;resn;charge;res+resn;flip;clash;Eint(kcal/mol);CapCorr(kcal/mol)\n")
        else:
            f.write("res;chain;resn;charge;res+resn;flip;Eint(kcal/mol);CapCorr(kcal/mol)\n")
        for result in results:
            base_name,chain_id,res_num,is_flipped = _parse_residue_id(result.residue_name,result.residue_state)
            #normalize charge/state formatting to match in_pocket outputs
            q = _parse_charge_state(getattr(result,"residue_state",None))
            proton_state = _format_charge_state(q)
            display_name = base_name
            if edda:
                clash_flag = _compute_clash_flag(result)
                f.write(f"{display_name};{chain_id};{res_num};{proton_state};"
                        f"{display_name}{res_num};{is_flipped};{str(clash_flag).upper()};"
                        f"{result.binding_energy:.6f};{result.cap_correction:.6f}\n")
            else:
                f.write(f"{display_name};{chain_id};{res_num};{proton_state};"
                        f"{display_name}{res_num};{is_flipped};"
                        f"{result.binding_energy:.6f};{result.cap_correction:.6f}\n")

def _write_best_sdf(output_dir: Path,geometry_folder: str,ligand,ipa_binding_energy: float,residues: Sequence[str],filename: str = "best_lig.sdf",logger: Optional[logging.Logger] = None) -> None:
    """Build best_lig.sdf by taking the original SDF file for the best ligand state,
       replacing coordinates with those from the IP-optimized XYZ, and appending
       SMILES / IPA_EBIND / LIST_RES blocks.
       The SDF formatting (columns, number of decimals, etc.) is preserved."""
    geometry_root = Path(output_dir) / geometry_folder
    sdf_hint = Path(ligand.file_path)
    stem = sdf_hint.stem
    #try a few candidate locations for the SDF, in order:
    #   1) geometries/start_sdf/{stem}.sdf
    #   2) ligand.file_path itself if it's an existing .sdf
    #   3) geometries/{stem}.sdf  (last resort / backward compatibility)
    candidates = [geometry_root / "start_sdf" / f"{stem}.sdf"]
    if sdf_hint.suffix.lower() == ".sdf":
        candidates.append(sdf_hint)
    candidates.append(geometry_root / f"{stem}.sdf")
    sdf_path = None
    for cand in candidates:
        if cand.is_file():
            sdf_path = cand
            break
    if sdf_path is None:
        msg = (f"Could not locate SDF for best ligand. Tried:\n" + "\n".join(str(c) for c in candidates))
        if logger:
            logger.error(msg)
        raise FileNotFoundError(msg)
    #use geometries/ip as IP directory regardless of where the SDF lives.
    #(this matches both serial and parallel layouts.)
    ip_dir = geometry_root / "ip"
    ip_xyz = ip_dir / f"{sdf_path.stem}_ip.xyz"
    use_ip_coords = True
    if not ip_xyz.is_file():
        if logger:
            logger.warning(f"IP-optimized structure {ip_xyz} not found; "
                           f"best_lig.sdf will keep original SDF coordinates.")
        use_ip_coords = False
    #read original SDF
    sdf_lines = sdf_path.read_text().splitlines(keepends = True)
    #locate atom and bond blocks
    if len(sdf_lines) < 4:
        raise ValueError(f"SDF file {sdf_path} is too short / malformed")
    counts_line = sdf_lines[3]
    try:
        num_atoms = int(counts_line[0:3])
    except ValueError as exc:
        raise ValueError(f"Cannot parse counts line in {sdf_path}: {counts_line!r}") from exc
    atom_start = 4
    atom_end = atom_start + num_atoms
    if len(sdf_lines) < atom_end:
        raise ValueError(f"Not enough atom lines in SDF file {sdf_path}")
    atom_lines = sdf_lines[atom_start:atom_end]
    #find "M  END" and "$$$$" indices
    m_end_idx: Optional[int] = None
    mol_end_idx: Optional[int] = None
    for i in range(atom_end,len(sdf_lines)):
        stripped = sdf_lines[i].strip()
        if stripped.startswith("M  END") and m_end_idx is None:
            m_end_idx = i
        if stripped == "$$$$":
            mol_end_idx = i
            break
    if m_end_idx is None:
        raise ValueError(f"'M  END' not found in SDF file {sdf_path}")
    if mol_end_idx is None:
        mol_end_idx = len(sdf_lines)    #we will add "$$$$" ourselves
    bond_and_other = sdf_lines[atom_end:m_end_idx]
    m_end_line = sdf_lines[m_end_idx]
    #replace coordinates if IP structure is available
    if use_ip_coords:
        xyz_lines = ip_xyz.read_text().splitlines()
        try:
            num_xyz_atoms = int(xyz_lines[0].strip())
        except ValueError as exc:
            raise ValueError(f"Cannot parse atom count in XYZ file {ip_xyz}") from exc
        if num_xyz_atoms != num_atoms:
            msg = (f"Atom count mismatch between SDF ({num_atoms}) and IP XYZ "
                   f"({num_xyz_atoms}) for {sdf_path}")
            if logger:
                logger.error(msg)
            raise ValueError(msg)
        coords = []
        for i in range(num_xyz_atoms):
            parts = xyz_lines[2 + i].split()
            if len(parts) < 4:
                raise ValueError(f"Malformed XYZ line {i+1} in {ip_xyz}: {xyz_lines[2+i]!r}")
            x,y,z = map(float,parts[1:4])
            coords.append((x,y,z))
        new_atom_lines = []
        for (x,y,z), old in zip(coords,atom_lines):
            #coordinates occupy columns 1–10, 11–20, 21–30 in standard V2000
            prefix = f"{x:10.4f}{y:10.4f}{z:10.4f}"
            rest = old[30:] if len(old) > 30 else "\n"
            new_atom_lines.append(prefix + rest)
    else:
        new_atom_lines = atom_lines
    #assemble new SDF molecule:
    #   - keep header and counts line (first 4)
    #   - replaced atom block
    #   - original bond block (and any lines between atom block and M  END)
    #   - M  END
    new_lines: List[str] = []
    new_lines.extend(sdf_lines[:atom_start])     #header + counts
    new_lines.extend(new_atom_lines)             #atoms
    new_lines.extend(bond_and_other)             #bonds and anything up to M  END
    new_lines.append(m_end_line if m_end_line.endswith("\n") else m_end_line + "\n")
    #drop any existing property block between M  END and $$$$ and replace it with our three fields
    #SMILES might be empty if non-RDKit protonation was used; handle gracefully
    smiles = getattr(ligand,"smiles","") or ""
    new_lines.append(">  <SMILES>\n")
    new_lines.append(smiles + "\n")
    new_lines.append(">  <IPA_EBIND>\n")
    new_lines.append(f"{ipa_binding_energy:.3f}\n")
    new_lines.append(">  <LIST_RES>\n")
    new_lines.append(",".join(residues) + "\n")
    new_lines.append("$$$$\n")
    out_path = Path(output_dir) / filename
    out_path.write_text("".join(new_lines))

def generate_best_outputs(
        output_dir: Path,
        geometry_folder: str,
        ligand_states: Sequence,
        binding_results: Sequence[Sequence],
        logger: Optional[logging.Logger] = None,
        compute_total_interaction_and_list_res: Optional[Callable[[Sequence], Tuple[float, List[str]]]] = None
    ) -> BestLigandSummary:
    """High-level helper to generate:
        - best_interaction_energies.csv
        - best_lig.sdf
       Returns a BestLigandSummary with the chosen ligand index and energies."""
    output_dir = Path(output_dir)
    #1) Compute total interaction and IPA binding energies for all states
    (ligand_int_energies,ipa_binding_energies,residues_per_ligand) = _compute_ligand_interactions_and_ipa(ligand_states,binding_results,logger,compute_total_interaction_and_list_res = compute_total_interaction_and_list_res)
    #2) Select best ligand (most negative IPA binding energy)
    if not ipa_binding_energies:
        raise ValueError("No IPA binding energies computed")
    best_idx = min(range(len(ipa_binding_energies)),key = lambda i: ipa_binding_energies[i])
    best_lig = ligand_states[best_idx]
    best_results = binding_results[best_idx]
    best_ipa = ipa_binding_energies[best_idx]
    best_int = ligand_int_energies[best_idx]
    best_residues = residues_per_ligand[best_idx]
    if logger:
        logger.info(f"[best_ligand] Best ligand index = {best_idx} "
                    f"(1-based: {best_idx+1}), IPA_binding = {best_ipa:.3f} kcal/mol")
    #3) Write CSV with interactions only for best ligand
    _write_best_interaction_csv(output_dir,best_lig,best_results)
    #4) Write SDF with IP coordinates + metadata
    _write_best_sdf(output_dir = output_dir,
                    geometry_folder = geometry_folder,
                    ligand = best_lig,
                    ipa_binding_energy = best_ipa,
                    residues = best_residues,
                    filename = "best_lig.sdf",
                    logger = logger)
    return BestLigandSummary(index = best_idx,ipa_binding_energy = best_ipa,total_interaction_energy = best_int,residues = best_residues)

def cleanup_outputs(output_dir: Path, geometry_folder: str = "geometries", output_folder_qm: str = "output", keep_debug: bool = False, logger: Optional[logging.Logger] = None) -> None:
    """Optional helper to clean the output folder so that, in the strict
       `--best` mode, *only* the two summary files remain.
       If keep_debug is True (e.g. for --best-debug), nothing is removed."""
    if keep_debug:
        if logger:
            logger.info("[best_ligand] keep_debug=True -> no cleanup performed")
        return
    output_dir = Path(output_dir)
    #files to remove if present
    extra_files = ["interaction_energies.csv","stats.txt","in_pocket_analysis.log"]
    for fname in extra_files:
        path = output_dir / fname
        if path.exists():
            try:
                path.unlink()
                if logger:
                    logger.info(f"[best_ligand] Removed file {path}")
            except Exception as e:
                if logger:
                    logger.warning(f"[best_ligand] Could not remove file {path}: {e}")
    #directories to remove recursively
    extra_dirs = [output_dir / geometry_folder,output_dir / output_folder_qm]
    for d in extra_dirs:
        if d.is_dir():
            try:
                shutil.rmtree(d)
                if logger:
                    logger.info(f"[best_ligand] Removed directory {d}")
            except Exception as e:
                if logger:
                    logger.warning(f"[best_ligand] Could not remove directory {d}: {e}")
