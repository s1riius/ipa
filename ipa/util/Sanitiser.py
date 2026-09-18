# ------------------------------------------------------------
# Sanitiser
# ------------------------------------------------------------

import logging
import re
from pathlib import Path
from typing import Dict,List,Optional,Tuple
import util.AtomicData as AD

def _get_logger(logger: Optional[logging.Logger] = None) -> logging.Logger:
    return logger if logger is not None else logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Generic text helpers
# ---------------------------------------------------------------------------
def _split_sdf_records(lines: List[str]) -> List[Tuple[List[str], bool]]:
    """Split an SDF into records.
       Returns a list of (record_lines_without_terminator, had_terminator)."""
    records: List[Tuple[List[str],bool]] = []
    current: List[str] = []
    for ln in lines:
        if ln.strip() == "$$$$":
            records.append((current,True))
            current = []
        else:
            current.append(ln)
    if current:
        records.append((current,False))
    return records

def _first_index(lines: List[str], predicate, start: int = 0) -> int:
    for i in range(start,len(lines)):
        if predicate(lines[i]):
            return i
    return -1

# ---------------------------------------------------------------------------
# V2000 helpers
# ---------------------------------------------------------------------------
def _parse_v2000_counts_line(line: str) -> Tuple[int, int]:
    try:
        return int(line[0:3]),int(line[3:6])
    except Exception as e:
        raise ValueError(f"Failed to parse V2000 counts line: {line!r}") from e

def _rebuild_v2000_counts_line(original: str, n_atoms: int, n_bonds: int) -> str:
    tail = original[6:] if len(original) > 6 else "  0  0  0  0  0  0  0  0999 V2000\n"
    return f"{n_atoms:>3}{n_bonds:>3}{tail}"

def _extract_v2000_element(atom_line: str) -> str:
    return (atom_line[31:34].strip() or atom_line[30:34].strip() or "?").upper()

def _rewrite_v2000_bond_line(original: str, a1_new: int, a2_new: int) -> str:
    tail = original[6:] if len(original) > 6 else "  1  0  0  0  0\n"
    return f"{a1_new:>3}{a2_new:>3}{tail}"

def _remap_v2000_m_record(line: str,atom_map: Dict[int,int]) -> Optional[str]:
    """Remap common atom-indexed V2000 M-records.
       Supports:
         M  CHG
         M  ISO
         M  RAD
       Unknown lines are returned unchanged.
       If a remapped record becomes empty, return None to drop it."""
    stripped = line.strip()
    if not stripped.startswith("M  "):
        return line
    parts = stripped.split()
    if len(parts) < 3:
        return line
    if parts[1] not in {"CHG","ISO","RAD"}:
        return line
    tag = parts[1]
    try:
        n_pairs = int(parts[2])
    except Exception:
        return line
    ints: List[int] = []
    for tok in parts[3:]:
        try:
            ints.append(int(tok))
        except Exception:
            return line
    if len(ints) < 2*n_pairs:
        return line
    kept_pairs: List[Tuple[int,int]] = []
    for i in range(n_pairs):
        old_idx = ints[2*i]
        val = ints[2*i + 1]
        new_idx = atom_map.get(old_idx)
        if new_idx is None:
            continue
        kept_pairs.append((new_idx,val))
    if not kept_pairs:
        return None
    rebuilt = f"M  {tag}{len(kept_pairs):>3}"
    for idx,val in kept_pairs:
        rebuilt += f"{idx:>4}{val:>4}"
    return rebuilt + "\n"

def _process_v2000_record(record_lines: List[str]) -> Tuple[List[str], int]:
    if len(record_lines) < 5:
        raise ValueError("V2000 record too short")
    counts = record_lines[3]
    n_atoms,n_bonds = _parse_v2000_counts_line(counts)
    atom_start = 4
    atom_end = atom_start + n_atoms
    bond_start = atom_end
    bond_end = bond_start + n_bonds
    if len(record_lines) < bond_end:
        raise ValueError("V2000 record truncated before end of bond block")
    atom_lines = record_lines[atom_start:atom_end]
    bond_lines = record_lines[bond_start:bond_end]
    keep_old_indices: List[int] = []
    for i,ln in enumerate(atom_lines,start = 1):
        elem = _extract_v2000_element(ln)
        if not AD._is_hydrogen_symbol(elem):
            keep_old_indices.append(i)
    h_removed = n_atoms - len(keep_old_indices)
    if h_removed == 0:
        return record_lines,0
    if not keep_old_indices:
        raise ValueError("All atoms in SDF record are hydrogens after sanitization")
    atom_map: Dict[int,int] = {old_idx: new_idx for new_idx,old_idx in enumerate(keep_old_indices,start = 1)}
    new_atom_lines: List[str] = [atom_lines[old_idx - 1] for old_idx in keep_old_indices]
    new_bond_lines: List[str] = []
    for ln in bond_lines:
        try:
            a1_old = int(ln[0:3])
            a2_old = int(ln[3:6])
        except Exception as e:
            raise ValueError(f"Failed to parse V2000 bond line: {ln!r}") from e
        a1_new = atom_map.get(a1_old)
        a2_new = atom_map.get(a2_old)
        if a1_new is None or a2_new is None:
            continue
        new_bond_lines.append(_rewrite_v2000_bond_line(ln,a1_new,a2_new))
    m_end_idx = _first_index(record_lines,lambda s: s.strip() == "M  END",start = bond_end)
    if m_end_idx < 0:
        #no M END found; keep everything after bond block as "tail"
        pre_m_end = record_lines[bond_end:]
        tail_after_m_end: List[str] = []
        m_end_line = None
    else:
        pre_m_end = record_lines[bond_end:m_end_idx]
        m_end_line = record_lines[m_end_idx]
        tail_after_m_end = record_lines[m_end_idx + 1:]
    remapped_pre_m_end: List[str] = []
    for ln in pre_m_end:
        remapped = _remap_v2000_m_record(ln,atom_map)
        if remapped is None:
            continue
        remapped_pre_m_end.append(remapped)
    out: List[str] = []
    out.extend(record_lines[:4])
    out[3] = _rebuild_v2000_counts_line(record_lines[3],n_atoms = len(new_atom_lines),n_bonds = len(new_bond_lines))
    out.extend(new_atom_lines)
    out.extend(new_bond_lines)
    out.extend(remapped_pre_m_end)
    if m_end_line is not None:
        out.append(m_end_line)
    out.extend(tail_after_m_end)
    return out,h_removed

# ---------------------------------------------------------------------------
# V3000 helpers
# ---------------------------------------------------------------------------
_V30_PREFIX = "M  V30 "

def _v30_body(line: str) -> str:
    s = line.rstrip("\n")
    if s.startswith(_V30_PREFIX):
        return s[len(_V30_PREFIX):]
    return s

def _rewrite_v30_counts_line(line: str, n_atoms: int, n_bonds: int) -> str:
    body = _v30_body(line)
    body = re.sub(r"\bCOUNTS\s+\d+\s+\d+",f"COUNTS {n_atoms} {n_bonds}",body,count = 1)
    return _V30_PREFIX + body + "\n"

def _remap_v30_collection_line(line: str, atom_map: Dict[int, int], bond_map: Dict[int, int]) -> Optional[str]:
    """Remap common V3000 collection lines containing:
         ATOMS=(n ...)
         BONDS=(n ...)
       If a collection becomes empty after remapping, return None."""
    body = _v30_body(line)
    changed = False
    became_empty = False

    def _replace_block(match, which_map):
        nonlocal changed,became_empty
        items_raw = match.group(2).strip()
        if not items_raw:
            kept: List[int] = []
        else:
            kept = []
            for tok in items_raw.split():
                try:
                    old_idx = int(tok)
                except Exception:
                    #leave original block unchanged if parsing is unexpected
                    return match.group(0)
                new_idx = which_map.get(old_idx)
                if new_idx is not None:
                    kept.append(new_idx)
        changed = True
        if not kept:
            became_empty = True
            return f"{match.group(1)}(0)"
        return f"{match.group(1)}({len(kept)} {' '.join(str(x) for x in kept)})"

    body = re.sub(r"(ATOMS=)\((\d+)\s*([^\)]*)\)",lambda m: _replace_block(re.match(r"(ATOMS=)\((.*)\)", f"ATOMS=({m.group(2)} {m.group(3)})"),atom_map),body)
    body = re.sub(r"(BONDS=)\((\d+)\s*([^\)]*)\)",lambda m: _replace_block(re.match(r"(BONDS=)\((.*)\)", f"BONDS=({m.group(2)} {m.group(3)})"),bond_map),body)
    if changed and became_empty:
        return None
    return _V30_PREFIX + body + "\n"

def _process_v3000_record(record_lines: List[str]) -> Tuple[List[str], int]:
    if len(record_lines) < 5:
        raise ValueError("V3000 record too short")
    counts_idx = _first_index(record_lines,lambda s: "M  V30 COUNTS" in s)
    atom_begin_idx = _first_index(record_lines,lambda s: s.strip() == "M  V30 BEGIN ATOM")
    atom_end_idx = _first_index(record_lines,lambda s: s.strip() == "M  V30 END ATOM")
    bond_begin_idx = _first_index(record_lines,lambda s: s.strip() == "M  V30 BEGIN BOND")
    bond_end_idx = _first_index(record_lines,lambda s: s.strip() == "M  V30 END BOND")
    if min(counts_idx,atom_begin_idx,atom_end_idx,bond_begin_idx,bond_end_idx) < 0:
        raise ValueError("Failed to locate mandatory V3000 sections")
    atom_lines = record_lines[atom_begin_idx + 1:atom_end_idx]
    bond_lines = record_lines[bond_begin_idx + 1:bond_end_idx]
    atom_map: Dict[int,int] = {}
    new_atom_lines: List[str] = []
    h_removed = 0
    next_atom_idx = 1
    for ln in atom_lines:
        body = _v30_body(ln).strip()
        toks = body.split()
        if len(toks) < 6:
            raise ValueError(f"Unexpected V3000 atom line: {ln!r}")
        old_idx = int(toks[0])
        elem = toks[1].upper()
        if AD._is_hydrogen_symbol(elem):
            h_removed += 1
            continue
        atom_map[old_idx] = next_atom_idx
        new_atom_lines.append(f"{_V30_PREFIX}{next_atom_idx} {' '.join(toks[1:])}\n")
        next_atom_idx += 1
    if h_removed == 0:
        return record_lines,0
    if not atom_map:
        raise ValueError("All atoms in V3000 record are hydrogens after sanitization")
    bond_map: Dict[int,int] = {}
    new_bond_lines: List[str] = []
    next_bond_idx = 1
    for ln in bond_lines:
        body = _v30_body(ln).strip()
        toks = body.split()
        if len(toks) < 4:
            raise ValueError(f"Unexpected V3000 bond line: {ln!r}")
        old_bond_idx = int(toks[0])
        order = toks[1]
        a1_old = int(toks[2])
        a2_old = int(toks[3])
        rest = toks[4:]
        a1_new = atom_map.get(a1_old)
        a2_new = atom_map.get(a2_old)
        if a1_new is None or a2_new is None:
            continue
        bond_map[old_bond_idx] = next_bond_idx
        rebuilt = f"{_V30_PREFIX}{next_bond_idx} {order} {a1_new} {a2_new}"
        if rest:
            rebuilt += " " + " ".join(rest)
        new_bond_lines.append(rebuilt + "\n")
        next_bond_idx += 1
    collection_begin_idx = _first_index(record_lines,lambda s: s.strip() == "M  V30 BEGIN COLLECTION",start = bond_end_idx + 1)
    collection_end_idx = _first_index(record_lines,lambda s: s.strip() == "M  V30 END COLLECTION",start = bond_end_idx + 1)
    out: List[str] = []
    out.extend(record_lines[:counts_idx])
    out.append(_rewrite_v30_counts_line(record_lines[counts_idx],len(new_atom_lines),len(new_bond_lines)))
    out.extend(record_lines[counts_idx + 1:atom_begin_idx + 1])
    out.extend(new_atom_lines)
    out.append(record_lines[atom_end_idx])
    out.extend(record_lines[atom_end_idx + 1:bond_begin_idx + 1])
    out.extend(new_bond_lines)
    out.append(record_lines[bond_end_idx])
    tail_start = bond_end_idx + 1
    if collection_begin_idx >= 0 and collection_end_idx >= 0 and collection_begin_idx < collection_end_idx:
        out.extend(record_lines[tail_start:collection_begin_idx + 1])
        for ln in record_lines[collection_begin_idx + 1:collection_end_idx]:
            remapped = _remap_v30_collection_line(ln,atom_map,bond_map)
            if remapped is None:
                continue
            out.append(remapped)
        out.append(record_lines[collection_end_idx])
        tail_start = collection_end_idx + 1
    out.extend(record_lines[tail_start:])
    return out,h_removed

# ---------------------------------------------------------------------------
# Ligand formal-charge helpers (special-case repair + total-charge extraction)
# ---------------------------------------------------------------------------
_V2000_CHARGE_CODE_TO_FORMAL: Dict[int, int] = {0: 0,
                                                1: 3,
                                                2: 2,
                                                3: 1,
                                                5: -1,
                                                6: -2,
                                                7: -3}

def _bond_order_value(order: int) -> int:
    """Map SDF bond order tokens to integer valence contributions."""
    try:
        o = int(order)
    except Exception:
        return 0
    if o in (1,2,3):
        return o
    #treat unknown / aromatic bond tokens conservatively
    return 0

def _parse_v2000_atom_charge_code(atom_line: str) -> int:
    """Parse the legacy MDL atom-block charge code from a V2000 atom line."""
    try:
        code = int(atom_line[36:39])
    except Exception:
        return 0
    return _V2000_CHARGE_CODE_TO_FORMAL.get(code,0)

def _zero_v2000_atom_charge_code(atom_line: str) -> str:
    """Normalise a V2000 atom line so that formal charges live only in M  CHG."""
    if len(atom_line) < 39:
        return atom_line
    return atom_line[:36] + "  0" + atom_line[39:]

def _parse_v2000_m_chg_lines(lines: List[str]) -> Dict[int, int]:
    """Parse V2000 M  CHG records into {atom_index: formal_charge}."""
    charges: Dict[int,int] = {}
    for line in lines:
        stripped = line.strip()
        if not stripped.startswith("M  CHG"):
            continue
        parts = stripped.split()
        if len(parts) < 3:
            continue
        try:
            n_pairs = int(parts[2])
        except Exception:
            continue
        raw_ints: List[int] = []
        try:
            for tok in parts[3:]:
                raw_ints.append(int(tok))
        except Exception:
            continue
        if len(raw_ints) < (2*n_pairs):
            continue
        for i in range(n_pairs):
            atom_idx = raw_ints[2*i]
            formal_charge = raw_ints[2*i + 1]
            charges[int(atom_idx)] = int(formal_charge)
    return charges

def _format_v2000_m_chg_lines(charges: Dict[int, int]) -> List[str]:
    """Build normalised V2000 M  CHG lines (max 8 atom/charge pairs per line)."""
    non_zero = sorted((int(atom_idx), int(formal_charge))
                      for atom_idx, formal_charge in charges.items()
                      if int(formal_charge) != 0)
    out: List[str] = []
    while non_zero:
        chunk = non_zero[:8]
        non_zero = non_zero[8:]
        line = f"M  CHG{len(chunk):>3}"
        for atom_idx,formal_charge in chunk:
            line += f"{atom_idx:>4}{formal_charge:>4}"
        out.append(line + "\n")
    return out

def _parse_v2000_record_graph(record_lines: List[str]) -> Tuple[List[str], List[Tuple[int, int, int]], Dict[int, int]]:
    """Extract elements, bonds and formal charges from a V2000 molblock."""
    if len(record_lines) < 5:
        raise ValueError("V2000 record too short")
    counts = record_lines[3]
    n_atoms,n_bonds = _parse_v2000_counts_line(counts)
    atom_start = 4
    atom_end = atom_start + n_atoms
    bond_start = atom_end
    bond_end = bond_start + n_bonds
    if len(record_lines) < bond_end:
        raise ValueError("V2000 record truncated before end of bond block")
    atom_lines = record_lines[atom_start:atom_end]
    bond_lines = record_lines[bond_start:bond_end]
    atoms: List[str] = []
    charges: Dict[int,int] = {}
    for atom_idx,atom_line in enumerate(atom_lines,start = 1):
        atoms.append(_extract_v2000_element(atom_line))
        atom_charge = _parse_v2000_atom_charge_code(atom_line)
        if atom_charge:
            charges[atom_idx] = atom_charge
    m_end_idx = _first_index(record_lines,lambda s: s.strip() == "M  END",start = bond_end)
    if m_end_idx < 0:
        pre_m_end = record_lines[bond_end:]
    else:
        pre_m_end = record_lines[bond_end:m_end_idx]
    charges.update(_parse_v2000_m_chg_lines(pre_m_end))
    bonds: List[Tuple[int,int,int]] = []
    for bond_line in bond_lines:
        try:
            a1 = int(bond_line[0:3])
            a2 = int(bond_line[3:6])
            order = int(bond_line[6:9])
        except Exception as e:
            raise ValueError(f"Failed to parse V2000 bond line: {bond_line!r}") from e
        bonds.append((a1,a2,order))
    return atoms,bonds,charges

def _rewrite_v2000_record_formal_charges(record_lines: List[str], charges: Dict[int, int]) -> List[str]:
    """Rewrite a V2000 record with normalised M  CHG lines."""
    if len(record_lines) < 5:
        raise ValueError("V2000 record too short")
    counts = record_lines[3]
    n_atoms,n_bonds = _parse_v2000_counts_line(counts)
    atom_start = 4
    atom_end = atom_start + n_atoms
    bond_start = atom_end
    bond_end = bond_start + n_bonds
    if len(record_lines) < bond_end:
        raise ValueError("V2000 record truncated before end of bond block")
    atom_lines = [
        _zero_v2000_atom_charge_code(line)
        for line in record_lines[atom_start:atom_end]
    ]
    bond_lines = record_lines[bond_start:bond_end]
    m_end_idx = _first_index(record_lines,lambda s: s.strip() == "M  END",start = bond_end)
    if m_end_idx < 0:
        pre_m_end = record_lines[bond_end:]
        tail_after_m_end: List[str] = []
        m_end_line = "M  END\n"
    else:
        pre_m_end = record_lines[bond_end:m_end_idx]
        m_end_line = record_lines[m_end_idx]
        tail_after_m_end = record_lines[m_end_idx + 1:]
    kept_pre_m_end: List[str] = []
    for line in pre_m_end:
        if line.strip().startswith("M  CHG"):
            continue
        kept_pre_m_end.append(line)
    out: List[str] = []
    out.extend(record_lines[:4])
    out.extend(atom_lines)
    out.extend(bond_lines)
    out.extend(kept_pre_m_end)
    out.extend(_format_v2000_m_chg_lines(charges))
    out.append(m_end_line)
    out.extend(tail_after_m_end)
    return out

def _parse_v3000_atom_charge(atom_line: str) -> int:
    """Parse CHG=<n> from a V3000 atom line."""
    for tok in _v30_body(atom_line).strip().split():
        if tok.startswith("CHG="):
            try:
                return int(tok.split("=", 1)[1])
            except Exception:
                return 0
    return 0

def _parse_v3000_record_graph(record_lines: List[str]) -> Tuple[List[str], List[Tuple[int, int, int]], Dict[int, int]]:
    """Extract elements, bonds and formal charges from a V3000 molblock."""
    if len(record_lines) < 5:
        raise ValueError("V3000 record too short")
    atom_begin_idx = _first_index(record_lines,lambda s: s.strip() == "M  V30 BEGIN ATOM")
    atom_end_idx = _first_index(record_lines,lambda s: s.strip() == "M  V30 END ATOM")
    bond_begin_idx = _first_index(record_lines,lambda s: s.strip() == "M  V30 BEGIN BOND")
    bond_end_idx = _first_index(record_lines,lambda s: s.strip() == "M  V30 END BOND")
    if min(atom_begin_idx,atom_end_idx,bond_begin_idx,bond_end_idx) < 0:
        raise ValueError("Failed to locate mandatory V3000 sections")
    atoms: List[str] = []
    charges: Dict[int,int] = {}
    for atom_line in record_lines[atom_begin_idx + 1:atom_end_idx]:
        toks = _v30_body(atom_line).strip().split()
        if len(toks) < 2:
            raise ValueError(f"Unexpected V3000 atom line: {atom_line!r}")
        atom_idx = int(toks[0])
        element = toks[1].upper()
        atoms.append(element)
        atom_charge = _parse_v3000_atom_charge(atom_line)
        if atom_charge:
            charges[atom_idx] = atom_charge
    bonds: List[Tuple[int,int,int]] = []
    for bond_line in record_lines[bond_begin_idx + 1:bond_end_idx]:
        toks = _v30_body(bond_line).strip().split()
        if len(toks) < 4:
            raise ValueError(f"Unexpected V3000 bond line: {bond_line!r}")
        order = int(toks[1])
        a1 = int(toks[2])
        a2 = int(toks[3])
        bonds.append((a1,a2,order))
    return atoms,bonds,charges

def _rewrite_v3000_atom_line_with_charge(atom_line: str, formal_charge: int) -> str:
    """Rewrite one V3000 atom line, normalising the CHG=<n> token."""
    toks = _v30_body(atom_line).strip().split()
    if len(toks) < 2:
        return atom_line
    kept = [tok for tok in toks[1:] if not tok.startswith("CHG=")]
    if int(formal_charge) != 0:
        kept.append(f"CHG={int(formal_charge)}")
    return _V30_PREFIX + " ".join([toks[0],*kept]) + "\n"

def _rewrite_v3000_record_formal_charges(record_lines: List[str], charges: Dict[int,int]) -> List[str]:
    """Rewrite a V3000 record with normalised CHG=<n> atom-line tokens."""
    atom_begin_idx = _first_index(record_lines,lambda s: s.strip() == "M  V30 BEGIN ATOM")
    atom_end_idx = _first_index(record_lines,lambda s: s.strip() == "M  V30 END ATOM")
    if min(atom_begin_idx,atom_end_idx) < 0:
        raise ValueError("Failed to locate V3000 atom block for charge rewrite")
    out: List[str] = []
    out.extend(record_lines[:atom_begin_idx + 1])
    for atom_line in record_lines[atom_begin_idx + 1:atom_end_idx]:
        atom_idx = int(_v30_body(atom_line).strip().split()[0])
        out.append(_rewrite_v3000_atom_line_with_charge(atom_line,charges.get(atom_idx,0)))
    out.extend(record_lines[atom_end_idx:])
    return out

def _build_adjacency(n_atoms: int, bonds: List[Tuple[int, int, int]]) -> Dict[int, List[Tuple[int, int, int]]]:
    """Build adjacency {atom_idx: [(nbr_idx, bond_order, bond_idx), ...]}."""
    adjacency: Dict[int,List[Tuple[int,int,int]]] = {i: [] for i in range(1,n_atoms + 1)}
    for bond_idx,(a1,a2,order) in enumerate(bonds,start = 1):
        adjacency[a1].append((a2,order,bond_idx))
        adjacency[a2].append((a1,order,bond_idx))
    return adjacency

def _rule_positive_nitrogen_valence4(atom_idx: int, atoms: List[str], adjacency: Dict[int, List[Tuple[int, int, int]]], charges: Dict[int, int]) -> Dict[int, int]:
    """Repair common ligand nitrogen motifs with bond-order sum 4.
       Covered examples:
         - N with 4 single bonds -> N(+1)
         - N with 2 singles + 1 double -> N(+1)
         - nitro / nitrate-like centres:
             N(=O)(-O)R
             N(=O)(-O)(-O)
             N(=O)R2"""
    if atoms[atom_idx - 1] != "N":
        return {}
    neigh = adjacency.get(atom_idx,[])
    if not neigh:
        return {}
    #avoid aromatic/delocalised ring nitrogens in the first implementation.
    if any(int(order) == 4 for _,order,_ in neigh):
        return {}
    bond_order_sum = sum(_bond_order_value(order) for _,order,_ in neigh)
    if bond_order_sum != 4:
        return {}
    updates: Dict[int,int] = {}
    current_n_charge = int(charges.get(atom_idx,0))
    if current_n_charge == 0:
        updates[atom_idx] = 1
    elif current_n_charge != 1:
        return {}
    for nbr_idx,order,_ in neigh:
        if atoms[nbr_idx - 1] != "O":
            continue
        if int(order) != 1:
            continue
        current_o_charge = int(charges.get(nbr_idx,0))
        if current_o_charge != 0:
            continue
        o_neigh = adjacency.get(nbr_idx,[])
        if len(o_neigh) == 1:
            updates[nbr_idx] = -1
    return updates

_SPECIAL_LIGAND_CHARGE_RULES = (("positive_nitrogen_valence4",_rule_positive_nitrogen_valence4),)

def _apply_special_ligand_charge_rules(atoms: List[str], bonds: List[Tuple[int, int, int]], charges: Dict[int, int]) -> Tuple[Dict[int, int], List[str], int]:
    """Apply modular special-case formal-charge rules to one molecular graph."""
    repaired = {int(k): int(v) for k,v in charges.items() if int(v) != 0}
    adjacency = _build_adjacency(len(atoms),bonds)
    atoms_changed: set = set()
    rule_hits: List[str] = []
    for atom_idx in range(1,len(atoms) + 1):
        for rule_name,rule_fn in _SPECIAL_LIGAND_CHARGE_RULES:
            updates = rule_fn(atom_idx,atoms,adjacency,repaired)
            if not updates:
                continue
            actually_changed = False
            for upd_idx,upd_charge in updates.items():
                upd_idx = int(upd_idx)
                upd_charge = int(upd_charge)
                if int(repaired.get(upd_idx,0)) == upd_charge:
                    continue
                repaired[upd_idx] = upd_charge
                atoms_changed.add(upd_idx)
                actually_changed = True
            if actually_changed:
                rule_hits.append(rule_name)
    return repaired,rule_hits,len(atoms_changed)

def repair_ligand_sdf_formal_charges(sdf_path: Path, logger: Optional[logging.Logger] = None) -> Dict[str, int]:
    """Repair special-case ligand formal charges directly in an SDF file.
       Behaviour:
         - preserve already assigned non-zero formal charges
         - fill in missing formal charges for selected graph motifs
         - rewrite charges in a normalised representation"""
    log = _get_logger(logger)
    path = Path(sdf_path)
    summary = {"records_processed": 0,
               "records_changed": 0,
               "atoms_changed": 0,
               "total_formal_charge": 0}
    try:
        if not path.exists():
            return summary
        lines = path.read_text().splitlines(True)
        records = _split_sdf_records(lines)
        out_lines: List[str] = []
        for record_lines,had_terminator in records:
            if (not record_lines) or all((not ln.strip()) for ln in record_lines):
                continue
            if len(record_lines) < 4:
                out_lines.extend(record_lines)
                if had_terminator:
                    out_lines.append("$$$$\n")
                continue
            counts_line = record_lines[3]
            fmt = None
            atoms: List[str] = []
            bonds: List[Tuple[int,int,int]] = []
            charges: Dict[int,int] = {}
            if "V3000" in counts_line:
                fmt = "V3000"
                atoms,bonds,charges = _parse_v3000_record_graph(record_lines)
            elif "V2000" in counts_line:
                fmt = "V2000"
                atoms,bonds,charges = _parse_v2000_record_graph(record_lines)
            else:
                out_lines.extend(record_lines)
                if had_terminator:
                    out_lines.append("$$$$\n")
                continue
            repaired,rule_hits,atoms_changed = _apply_special_ligand_charge_rules(atoms,bonds,charges)
            summary["records_processed"] += 1
            summary["atoms_changed"] += atoms_changed
            summary["total_formal_charge"] += sum(int(v) for v in repaired.values())
            if atoms_changed > 0:
                if fmt == "V3000":
                    new_record = _rewrite_v3000_record_formal_charges(record_lines,repaired)
                else:
                    new_record = _rewrite_v2000_record_formal_charges(record_lines,repaired)
                out_lines.extend(new_record)
                summary["records_changed"] += 1
                log.info("Repaired %d ligand formal-charge atom(s) in %s using rules: %s",atoms_changed,path.name,", ".join(sorted(set(rule_hits))) if rule_hits else "<?>")
            else:
                out_lines.extend(record_lines)
            if had_terminator:
                out_lines.append("$$$$\n")
        if summary["records_changed"] > 0:
            path.write_text("".join(out_lines))
    except Exception as e:
        log.warning("Failed to repair ligand formal charges in %s: %s",path,e)
    return summary

def get_sdf_total_formal_charge(sdf_path: Path, logger: Optional[logging.Logger] = None) -> Optional[int]:
    """Return the total formal charge encoded in an SDF, if parseable."""
    log = _get_logger(logger)
    path = Path(sdf_path)
    try:
        if not path.exists():
            return None
        lines = path.read_text().splitlines(True)
        records = _split_sdf_records(lines)
        total_charge = 0
        parsed_any = False
        for record_lines,_ in records:
            if (not record_lines) or all((not ln.strip()) for ln in record_lines):
                continue
            if len(record_lines) < 4:
                continue
            counts_line = record_lines[3]
            if "V3000" in counts_line:
                _,_,charges = _parse_v3000_record_graph(record_lines)
            elif "V2000" in counts_line:
                _,_,charges = _parse_v2000_record_graph(record_lines)
            else:
                continue
            total_charge += sum(int(v) for v in charges.values())
            parsed_any = True
        if not parsed_any:
            return None
        return int(total_charge)
    except Exception as e:
        log.debug("Failed to parse total formal charge from %s: %s",path,e)
        return None

# ---------------------------------------------------------------------------
# Public functions
# ---------------------------------------------------------------------------
def strip_h_from_pdb_lines(pdb_lines: List[str]) -> Tuple[List[str], int]:
    """Remove ATOM/HETATM hydrogen records from PDB lines.
       Detection is based on element column (77-78) when present, with a fallback
       to atom name (columns 13-16).
       Returns (filtered_lines, n_removed)."""
    out: List[str] = []
    removed = 0
    for ln in pdb_lines:
        rec = ln[0:6].strip()
        if rec in ("ATOM","HETATM"):
            elem = ln[76:78].strip() if len(ln) >= 78 else ""
            atname = ln[12:16].strip() if len(ln) >= 16 else ""
            is_h = (elem.upper() == "H")
            if (not elem) and atname.upper().startswith("H"):
                is_h = True
            if is_h:
                removed += 1
                continue
        out.append(ln)
    return out,removed

def strip_h_from_sdf_file(sdf_in: str, sdf_out: str) -> Tuple[int, int]:
    """Remove explicit H atoms from all molecules in an SDF without RDKit.
       Philosophy:
         - preserve original heavy-atom order
         - preserve original heavy-heavy bond orders/connectivity
         - remove only explicit H atom records and bonds touching them
         - remap common atom-indexed metadata records
         - preserve multi-record SDF structure and property/data blocks
       Supported:
         - V2000
         - common V3000 ligand molblocks"""
    in_path = Path(sdf_in)
    out_path = Path(sdf_out)
    if not in_path.exists():
        raise FileNotFoundError(f"SDF file not found: {sdf_in}")
    lines = in_path.read_text().splitlines(True)
    records = _split_sdf_records(lines)
    out_lines: List[str] = []
    n_written = 0
    h_removed_total = 0
    for record_lines,had_terminator in records:
        if (not record_lines) or all((not ln.strip()) for ln in record_lines):
            continue
        if len(record_lines) < 4:
            raise ValueError("SDF record too short to determine format")
        counts_line = record_lines[3]
        if "V3000" in counts_line:
            new_record,h_removed = _process_v3000_record(record_lines)
        elif "V2000" in counts_line:
            new_record,h_removed = _process_v2000_record(record_lines)
        else:
            raise ValueError(f"Unsupported/unknown SDF molblock format in record counts line: {counts_line!r}")
        out_lines.extend(new_record)
        if had_terminator:
            out_lines.append("$$$$\n")
        n_written += 1
        h_removed_total += h_removed
    out_path.write_text("".join(out_lines))
    return n_written,h_removed_total

def sanitize_ligand_file_for_protonation(ligand_path: str, output_folder: Path, keep_input_h: bool, logger: logging.Logger) -> str:
    """Sanitize ligand input file before ANY protonation backend.
       - If keep_input_h: return original path.
       - For SDF: write a sanitized copy into output_folder, optionally stripping
         explicit H atoms, then repair special-case ligand formal charges.
       - For PDB: nothing to do here (protein parsing already strips H); return original path."""
    lp = str(ligand_path)
    lower = lp.lower()
    if not lower.endswith(".sdf"):
        return ligand_path
    if keep_input_h:
        out_path = lp
    else:
        out_path = str(Path(output_folder) / (Path(lp).stem + "_noH.sdf"))
        try:
            n_mols,h_removed = strip_h_from_sdf_file(lp,out_path)
            logger.info("Stripped ~%d explicit H atoms from ligand SDF (%d mols) -> %s (use --keep-input-h to disable).",h_removed,n_mols,out_path)
        except Exception as e:
            logger.warning("Failed to sanitize ligand SDF (keeping original): %s",e)
            out_path = ligand_path
    repair_summary = repair_ligand_sdf_formal_charges(Path(out_path),logger = logger)
    if repair_summary.get("atoms_changed",0):
        logger.info("Applied special-case ligand charge repair to %s: atoms_changed = %d total_formal_charge = %d",
                    Path(out_path).name,
                    repair_summary.get("atoms_changed",0),
                    repair_summary.get("total_formal_charge",0))
    return out_path

def repair_native_sdf_h_connectivity(sdf_path: Path, logger: Optional[logging.Logger] = None) -> None:
    """Add missing explicit X-H bonds to a native/original SDF template.
       This is intentionally used only for original/native input processing,
       where explicit input hydrogens are preserved but the connectivity passed
       to WriteSDF(...) may still be a heavy-atom template."""
    log = _get_logger(logger)
    try:
        if not sdf_path.exists():
            return
        lines = sdf_path.read_text().splitlines(True)
        if len(lines) < 5:
            return
        counts = lines[3]
        if "V2000" not in counts:
            log.debug("Native SDF H-connectivity fix skipped for non-V2000 file: %s",sdf_path)
            return
        n_atoms = int(counts[0:3])
        n_bonds = int(counts[3:6])
        atom_start = 4
        bond_start = atom_start + n_atoms
        atoms_data = []
        for i in range(n_atoms):
            ln = lines[atom_start + i]
            x = float(ln[0:10])
            y = float(ln[10:20])
            z = float(ln[20:30])
            elem = (ln[31:34].strip() or ln[30:34].strip() or "?").upper()
            atoms_data.append((elem,x,y,z))
        bonded = {i: set() for i in range(1,n_atoms + 1)}
        existing_pairs = set()
        for i in range(n_bonds):
            ln = lines[bond_start + i]
            a1 = int(ln[0:3])
            a2 = int(ln[3:6])
            pair = tuple(sorted((a1,a2)))
            existing_pairs.add(pair)
            bonded[a1].add(a2)
            bonded[a2].add(a1)
        r_h = AD._covalent_radius_angstrom("H") or 0.31
        new_bonds: List[str] = []
        for h_idx,(elem_h,xh,yh,zh) in enumerate(atoms_data,start = 1):
            if elem_h != "H":
                continue
            if bonded[h_idx]:
                continue
            best_idx = None
            best_d2 = None
            for heavy_idx,(elem_x,xx,yx,zx) in enumerate(atoms_data,start = 1):
                if heavy_idx == h_idx:
                    continue
                if elem_x == "H":
                    continue
                r_x = AD._covalent_radius_angstrom(elem_x) or 0.76
                thr = 1.25*(r_h + r_x)
                d2 = (xh - xx)**2 + (yh - yx)**2 + (zh - zx)**2
                if d2 > thr*thr:
                    continue
                if (best_d2 is None) or (d2 < best_d2):
                    best_d2 = d2
                    best_idx = heavy_idx
            if best_idx is None:
                continue
            pair = tuple(sorted((h_idx,best_idx)))
            if pair in existing_pairs:
                continue
            new_bonds.append(f"{h_idx:3d}{best_idx:3d}{1:3d}  0  0  0  0\n")
            existing_pairs.add(pair)
            bonded[h_idx].add(best_idx)
            bonded[best_idx].add(h_idx)
        if not new_bonds:
            return
        lines[3] = f"{n_atoms:3d}{n_bonds + len(new_bonds):3d}{counts[6:]}"
        lines[bond_start + n_bonds: bond_start + n_bonds] = new_bonds
        sdf_path.write_text("".join(lines))
        log.info("Added %d missing X-H bond(s) to native/original SDF template: %s",len(new_bonds),sdf_path.name)
    except Exception as e:
        log.warning("Failed to add missing X-H bonds to native/original SDF template %s: %s",sdf_path,e)

def repair_prepared_sdf_h_connectivity(sdf_path: Path, logger: Optional[logging.Logger] = None) -> None:
    """Repair missing explicit X-H bonds in a prepared ligand/residue SDF
       immediately before the final geometry/connectivity validation.
       This intentionally reuses the same geometric X-H assignment logic used by
       repair_native_sdf_h_connectivity(), but is applied generically to any
       prepared SDF that will actually be used downstream."""
    path = Path(sdf_path)
    if path.suffix.lower() != ".sdf":
        return
    repair_native_sdf_h_connectivity(path,logger = logger)

def _read_xyz_atoms_and_geometry(xyz_path: Path) -> Tuple[List[str], List[Tuple[float, float, float]]]:
    """Read a simple XYZ file into (elements, geometry)."""
    path = Path(xyz_path)
    if not path.exists():
        raise FileNotFoundError(f"XYZ file not found: {path}")
    raw = path.read_text().splitlines()
    if len(raw) < 3:
        raise ValueError(f"XYZ file too short: {path}")
    try:
        n_atoms = int(raw[0].strip())
    except Exception as e:
        raise ValueError(f"Failed to parse XYZ atom count from {path}") from e
    body = raw[2:2 + n_atoms]
    if len(body) < n_atoms:
        raise ValueError(f"XYZ file truncated: expected {n_atoms} atoms in {path}")
    atoms: List[str] = []
    geom: List[Tuple[float,float,float]] = []
    for ln in body:
        parts = ln.split()
        if len(parts) < 4:
            raise ValueError(f"Malformed XYZ atom line in {path}: {ln!r}")
        elem = str(parts[0]).strip().upper()
        x = float(parts[1])
        y = float(parts[2])
        z = float(parts[3])
        atoms.append(elem)
        geom.append((x,y,z))
    return atoms,geom

def _read_pdb_atoms_and_geometry(pdb_path: Path) -> Tuple[List[str], List[Tuple[float, float, float]]]:
    """Read a simple PDB file into (elements, geometry)."""
    path = Path(pdb_path)
    if not path.exists():
        raise FileNotFoundError(f"PDB file not found: {path}")
    atoms: List[str] = []
    geom: List[Tuple[float,float,float]] = []
    for ln in path.read_text().splitlines():
        if not ln.startswith(("ATOM","HETATM")):
            continue
        elem = (ln[76:78].strip() if len(ln) >= 78 else "").upper()
        if not elem:
            atname = ln[12:16].strip().upper() if len(ln) >= 16 else ""
            elem = atname[:1]
        try:
            x = float(ln[30:38])
            y = float(ln[38:46])
            z = float(ln[46:54])
        except Exception as e:
            raise ValueError(f"Malformed PDB coordinates in {path}: {ln!r}") from e
        atoms.append(elem)
        geom.append((x,y,z))
    if not atoms:
        raise ValueError(f"No ATOM/HETATM records found in {path}")
    return atoms,geom

def _parse_v2000_record_atoms_and_geometry(record_lines: List[str]) -> Tuple[List[str], List[Tuple[float, float, float]]]:
    if len(record_lines) < 5:
        raise ValueError("V2000 record too short")
    counts = record_lines[3]
    n_atoms,_ = _parse_v2000_counts_line(counts)
    atom_start = 4
    atom_end = atom_start + n_atoms
    if len(record_lines) < atom_end:
        raise ValueError("V2000 record truncated before end of atom block")
    atoms: List[str] = []
    geom: List[Tuple[float,float,float]] = []
    for atom_line in record_lines[atom_start:atom_end]:
        try:
            x = float(atom_line[0:10])
            y = float(atom_line[10:20])
            z = float(atom_line[20:30])
        except Exception as e:
            raise ValueError(f"Failed to parse V2000 atom coordinates: {atom_line!r}") from e
        atoms.append(_extract_v2000_element(atom_line))
        geom.append((x,y,z))
    return atoms,geom

def _rewrite_v2000_atom_line_coords(atom_line: str, x: float, y: float, z: float) -> str:
    rest = atom_line[30:] if len(atom_line) > 30 else "\n"
    return f"{x:10.4f}{y:10.4f}{z:10.4f}{rest}"

def _parse_v3000_record_atoms_and_geometry(record_lines: List[str]) -> Tuple[List[str], List[Tuple[float, float, float]]]:
    if len(record_lines) < 5:
        raise ValueError("V3000 record too short")
    atom_begin_idx = _first_index(record_lines,lambda s: s.strip() == "M  V30 BEGIN ATOM")
    atom_end_idx = _first_index(record_lines,lambda s: s.strip() == "M  V30 END ATOM")
    if min(atom_begin_idx,atom_end_idx) < 0:
        raise ValueError("Failed to locate V3000 atom block")
    atoms: List[str] = []
    geom: List[Tuple[float,float,float]] = []
    for atom_line in record_lines[atom_begin_idx + 1:atom_end_idx]:
        toks = _v30_body(atom_line).strip().split()
        if len(toks) < 6:
            raise ValueError(f"Unexpected V3000 atom line: {atom_line!r}")
        atoms.append(toks[1].upper())
        geom.append((float(toks[3]),float(toks[4]),float(toks[5])))
    return atoms,geom

def _rewrite_v3000_atom_line_coords(atom_line: str, x: float, y: float, z: float) -> str:
    toks = _v30_body(atom_line).strip().split()
    if len(toks) < 6:
        return atom_line
    toks[3] = f"{float(x):.4f}"
    toks[4] = f"{float(y):.4f}"
    toks[5] = f"{float(z):.4f}"
    return _V30_PREFIX + " ".join(toks) + "\n"

def _build_geometry_adjacency(atoms: List[str], geom: List[Tuple[float, float, float]], scale: float = 1.30) -> Dict[int, List[int]]:
    adjacency: Dict[int,List[int]] = {i: [] for i in range(1,len(atoms) + 1)}
    for i in range(len(atoms)):
        ri = AD._covalent_radius_angstrom(atoms[i])
        if ri is None:
            continue
        xi,yi,zi = geom[i]
        for j in range(i + 1,len(atoms)):
            rj = AD._covalent_radius_angstrom(atoms[j])
            if rj is None:
                continue
            xj,yj,zj = geom[j]
            thr = float(scale)*(ri + rj)
            d2 = (xi - xj)**2 + (yi - yj)**2 + (zi - zj)**2
            if d2 <= thr*thr:
                adjacency[i + 1].append(j + 1)
                adjacency[j + 1].append(i + 1)
    return adjacency

def _adjacency_neighbor_indices(adjacency: Dict[int, List], idx: int) -> List[int]:
    """Normalise adjacency payloads to plain neighbour indices.
       _build_geometry_adjacency() stores [nbr_idx, ...], while _build_adjacency()
       stores [(nbr_idx, bond_order, bond_idx), ...]. This helper lets the
       validation/fixup routines work with both representations."""
    nbrs: List[int] = []
    for item in adjacency.get(idx,[]):
        if isinstance(item,tuple):
            if not item:
                continue
            nbrs.append(int(item[0]))
        else:
            nbrs.append(int(item))
    return nbrs

def _connected_components_from_adjacency(adjacency: Dict[int, List]) -> List[List[int]]:
    components: List[List[int]] = []
    seen = set()
    for start in adjacency:
        if start in seen:
            continue
        stack = [start]
        comp: List[int] = []
        seen.add(start)
        while stack:
            node = stack.pop()
            comp.append(node)
            for nbr in _adjacency_neighbor_indices(adjacency,node):
                if nbr in seen:
                    continue
                seen.add(nbr)
                stack.append(nbr)
        components.append(sorted(comp))
    return components

def _validate_connectivity_from_graph(atoms: List[str], adjacency: Dict[int, List], path_label: str) -> Dict[str, object]:
    n_atoms = len(atoms)
    isolated = [idx for idx in adjacency if len(_adjacency_neighbor_indices(adjacency,idx)) == 0]
    components = _connected_components_from_adjacency(adjacency)
    if n_atoms <= 1:
        return {"ok": True,"n_atoms": n_atoms,"isolated_indices": isolated,"n_components": len(components),"components": components}
    if isolated:
        labels = ", ".join(f"{atoms[idx - 1]}{idx}" for idx in isolated[:10])
        raise ValueError(f"Detected isolated atom(s) in prepared structure {path_label}: {labels}")
    if len(components) > 1:
        sizes = ", ".join(str(len(comp)) for comp in components)
        raise ValueError(f"Detected disconnected prepared structure {path_label}: "
                         f"{len(components)} components (sizes: {sizes})")
    return {"ok": True,"n_atoms": n_atoms,"isolated_indices": isolated,"n_components": len(components),"components": components}

def validate_no_isolated_atoms(structure_path: Path, logger: Optional[logging.Logger] = None) -> Dict[str, object]:
    """Fail if a prepared structure contains isolated atoms or disconnected components."""
    log = _get_logger(logger)
    path = Path(structure_path)
    suffix = path.suffix.lower()
    if suffix == ".xyz":
        print("REPLACE THIS WITH GEOMETRY.py xyz")
        atoms,geom = _read_xyz_atoms_and_geometry(path)
        adjacency = _build_geometry_adjacency(atoms,geom)
        result = _validate_connectivity_from_graph(atoms,adjacency,str(path))
        log.debug("Connectivity validation OK for %s: n_atoms = %d n_components = %d",path.name,result.get("n_atoms",0),result.get("n_components",0))
        return result
    if suffix == ".pdb":
        print("REPLACE THIS WITH GEOMETRY.py")
        atoms,geom = _read_pdb_atoms_and_geometry(path)
        adjacency = _build_geometry_adjacency(atoms,geom)
        result = _validate_connectivity_from_graph(atoms,adjacency,str(path))
        log.debug(
            "Connectivity validation OK for %s: n_atoms = %d n_components = %d",path.name,result.get("n_atoms",0),result.get("n_components",0))
        return result
    if suffix == ".sdf":
        lines = path.read_text().splitlines(True)
        records = _split_sdf_records(lines)
        processed = 0
        for rec_idx,(record_lines,_) in enumerate(records,start = 1):
            if (not record_lines) or all((not ln.strip()) for ln in record_lines):
                continue
            if len(record_lines) < 4:
                raise ValueError(f"SDF record too short while validating connectivity: {path}")
            counts_line = record_lines[3]
            if "V3000" in counts_line:
                atoms,bonds,_ = _parse_v3000_record_graph(record_lines)
            elif "V2000" in counts_line:
                atoms,bonds,_ = _parse_v2000_record_graph(record_lines)
            else:
                raise ValueError(f"Unsupported SDF format while validating connectivity: {path}")
            adjacency = _build_adjacency(len(atoms),bonds)
            _validate_connectivity_from_graph(atoms,adjacency,f"{path} [record {rec_idx}]")
            processed += 1
        log.debug("Connectivity validation OK for %s: records = %d",path.name,processed)
        return {"ok": True,"records_processed": processed}
    raise ValueError(f"Unsupported prepared-structure format for connectivity validation: {path}")

def _ideal_x_h_distance(elem_x: str) -> float:
    r_h = AD._covalent_radius_angstrom("H") or 0.31
    r_x = AD._covalent_radius_angstrom(elem_x) or 0.76
    return float(r_h + r_x)

def _is_bad_h_distance(actual: float, ideal: float) -> bool:
    try:
        d = float(actual)
    except Exception:
        return True
    if d <= 1.0e-8:
        return True
    if d < max(0.25,0.40*float(ideal)):
        return True
    if d > max(2.20,1.60*float(ideal)):
        return True
    return False

def _reposition_bonded_hydrogen(heavy_xyz: Tuple[float, float, float], hydrogen_xyz: Tuple[float, float, float], ideal_distance: float) -> Tuple[Tuple[float, float, float], bool]:
    xh,yh,zh = hydrogen_xyz
    xx,yx,zx = heavy_xyz
    vx = xh - xx
    vy = yh - yx
    vz = zh - zx
    norm2 = vx*vx + vy*vy + vz*vz
    if norm2 <= 1.0e-12:
        return (xx + ideal_distance,yx,zx),False
    norm = norm2**0.5
    scale = float(ideal_distance)/norm
    return (xx + vx*scale,yx + vy*scale,zx + vz*scale),True

#FMCM: this should go to Geometry.py
def _distance_xyz(a_xyz: Tuple[float, float, float], b_xyz: Tuple[float, float, float]) -> float:
    ax,ay,az = a_xyz
    bx,by,bz = b_xyz
    return ((ax - bx)**2 + (ay - by)**2 + (az - bz)**2)**0.5

def _is_plausible_x_h_contact(actual: float, ideal: float) -> bool:
    """Loose plausibility check for a geometrically reasonable X-H contact."""
    try:
        d = float(actual)
    except Exception:
        return False
    return d <= max(1.25,1.20*float(ideal))

def _find_nearest_heavy_atom(atoms: List[str], geom: List[Tuple[float, float, float]], h_idx: int) -> Tuple[Optional[int], Optional[float]]:
    best_idx: Optional[int] = None
    best_d = float("inf")
    if h_idx < 1 or h_idx > len(atoms):
        return None,None
    h_xyz = geom[h_idx - 1]
    for idx,elem in enumerate(atoms,start = 1):
        if idx == h_idx:
            continue
        if str(elem).strip().upper() == "H":
            continue
        d = _distance_xyz(h_xyz,geom[idx - 1])
        if d < best_d:
            best_d = d
            best_idx = idx
    if best_idx is None:
        return None,None
    return best_idx,float(best_d)

def _count_h_connectivity_shift_warnings_in_record(atoms: List[str], geom: List[Tuple[float, float, float]], adjacency: Dict[int, List], path: Path, rec_idx: int, logger: Optional[logging.Logger] = None) -> int:
    """Warn when post-QM geometry suggests an H atom is now closer to a different
       heavy atom than the one encoded in the SDF connectivity.
       This is intentionally warning-only for now: connectivity is preserved."""
    log = _get_logger(logger)
    warnings_count = 0
    for h_idx,elem in enumerate(atoms,start = 1):
        if str(elem).strip().upper() != "H":
            continue
        neigh = _adjacency_neighbor_indices(adjacency,h_idx)
        if len(neigh) != 1:
            continue
        bonded_idx = neigh[0]
        bonded_elem = atoms[bonded_idx - 1]
        if str(bonded_elem).strip().upper() == "H":
            continue
        nearest_idx,nearest_d = _find_nearest_heavy_atom(atoms,geom,h_idx)
        if nearest_idx is None or nearest_d is None:
            continue
        if nearest_idx == bonded_idx:
            continue
        bonded_d = _distance_xyz(geom[h_idx - 1],geom[bonded_idx - 1])
        bonded_ideal = _ideal_x_h_distance(bonded_elem)
        nearest_elem = atoms[nearest_idx - 1]
        nearest_ideal = _ideal_x_h_distance(nearest_elem)
        bonded_stretched = bonded_d >= max(1.35,1.15*float(bonded_ideal))
        nearest_plausible = _is_plausible_x_h_contact(nearest_d,nearest_ideal)
        clear_preference = (bonded_d - nearest_d) >= 0.25
        if not (bonded_stretched and nearest_plausible and clear_preference):
            continue
        warnings_count += 1
        log.warning("Post-QM H-connectivity warning in %s [record %d]: "
                    "H%d is bonded to %s%d in the SDF, but geometry places it closer to %s%d "
                    "(%0.3f Å vs bonded %0.3f Å). Connectivity was kept unchanged.",
                    path.name,rec_idx,h_idx,bonded_elem,bonded_idx,nearest_elem,nearest_idx,nearest_d,bonded_d)
    return warnings_count

def warn_sdf_h_connectivity_changes(sdf_path: Path, logger: Optional[logging.Logger] = None) -> Dict[str, int]:
    """Inspect an SDF and emit warning(s) when an H atom appears geometrically
       reassigned after optimisation, while the bond block still preserves the
       original connectivity."""
    path = Path(sdf_path)
    summary = {"records_processed": 0,"connectivity_shift_warnings": 0}
    if not path.exists():
        return summary
    lines = path.read_text().splitlines(True)
    records = _split_sdf_records(lines)
    for rec_idx,(record_lines,_) in enumerate(records,start = 1):
        if (not record_lines) or all((not ln.strip()) for ln in record_lines):
            continue
        if len(record_lines) < 4:
            raise ValueError(f"SDF record too short while checking H connectivity shifts: {path}")
        counts_line = record_lines[3]
        summary["records_processed"] += 1
        if "V3000" in counts_line:
            atoms,geom = _parse_v3000_record_atoms_and_geometry(record_lines)
            bonds = _parse_v3000_record_graph(record_lines)[1]
        elif "V2000" in counts_line:
            atoms,geom = _parse_v2000_record_atoms_and_geometry(record_lines)
            bonds = _parse_v2000_record_graph(record_lines)[1]
        else:
            raise ValueError(f"Unsupported SDF format while checking H connectivity shifts: {path}")
        adjacency = _build_adjacency(len(atoms),bonds)
        summary["connectivity_shift_warnings"] += _count_h_connectivity_shift_warnings_in_record(atoms,geom,adjacency,path,rec_idx,logger = logger)
    return summary

def check_and_fix_sdf_h_distances(sdf_path: Path, fix: bool = True, logger: Optional[logging.Logger] = None) -> Dict[str, int]:
    """Check bonded H distances in an SDF; optionally rewrite clearly wrong coordinates."""
    log = _get_logger(logger)
    path = Path(sdf_path)
    summary = {"records_processed": 0,"bad_h_distances": 0,"h_repositioned": 0}
    if not path.exists():
        return summary
    lines = path.read_text().splitlines(True)
    records = _split_sdf_records(lines)
    out_lines: List[str] = []
    changed_any = False
    for rec_idx,(record_lines,had_terminator) in enumerate(records,start = 1):
        if (not record_lines) or all((not ln.strip()) for ln in record_lines):
            continue
        if len(record_lines) < 4:
            raise ValueError(f"SDF record too short while checking H distances: {path}")
        counts_line = record_lines[3]
        summary["records_processed"] += 1
        if "V3000" in counts_line:
            atoms,geom = _parse_v3000_record_atoms_and_geometry(record_lines)
            bonds = _parse_v3000_record_graph(record_lines)[1]
            adjacency = _build_adjacency(len(atoms), bonds)
            atom_begin_idx = _first_index(record_lines,lambda s: s.strip() == "M  V30 BEGIN ATOM")
            atom_end_idx = _first_index(record_lines,lambda s: s.strip() == "M  V30 END ATOM")
            atom_lines = list(record_lines[atom_begin_idx + 1:atom_end_idx])
            local_changed = False
            for h_idx,elem in enumerate(atoms,start = 1):
                if elem != "H":
                    continue
                neigh = _adjacency_neighbor_indices(adjacency,h_idx)
                if len(neigh) != 1:
                    continue
                heavy_idx = neigh[0]
                if atoms[heavy_idx - 1] == "H":
                    continue
                ideal = _ideal_x_h_distance(atoms[heavy_idx - 1])
                #FMCM: why not use the distance function written?
                hx,hy,hz = geom[h_idx - 1]
                xx,yx,zx = geom[heavy_idx - 1]
                actual = ((hx - xx)**2 + (hy - yx)**2 + (hz - zx)**2)**0.5
                if not _is_bad_h_distance(actual,ideal):
                    continue
                summary["bad_h_distances"] += 1
                msg = (f"Detected suspicious X-H distance in {path.name} [record {rec_idx}]: "
                       f"H{h_idx}-{atoms[heavy_idx - 1]}{heavy_idx} = {actual:.3f} Å "
                       f"(target ~{ideal:.3f} Å)")
                if fix:
                    new_xyz,used_vector = _reposition_bonded_hydrogen(geom[heavy_idx - 1],geom[h_idx - 1],ideal)
                    geom[h_idx - 1] = new_xyz
                    atom_lines[h_idx - 1] = _rewrite_v3000_atom_line_coords(atom_lines[h_idx - 1],*new_xyz)
                    summary["h_repositioned"] += 1
                    local_changed = True
                    log.warning("%s; repositioned H atom coordinates%s.",msg," using existing bond vector" if used_vector else " with axis fallback")
                else:
                    log.warning("%s; native/original mode keeps the input coordinates unchanged.",msg)
            if local_changed:
                changed_any = True
                new_record = list(record_lines)
                new_record[atom_begin_idx + 1:atom_end_idx] = atom_lines
                out_lines.extend(new_record)
            else:
                out_lines.extend(record_lines)
        elif "V2000" in counts_line:
            atoms,geom = _parse_v2000_record_atoms_and_geometry(record_lines)
            bonds = _parse_v2000_record_graph(record_lines)[1]
            adjacency = _build_adjacency(len(atoms),bonds)
            n_atoms,_ = _parse_v2000_counts_line(counts_line)
            atom_start = 4
            atom_end = atom_start + n_atoms
            atom_lines = list(record_lines[atom_start:atom_end])
            local_changed = False
            for h_idx,elem in enumerate(atoms,start = 1):
                if elem != "H":
                    continue
                neigh = _adjacency_neighbor_indices(adjacency,h_idx)
                if len(neigh) != 1:
                    continue
                heavy_idx = neigh[0]
                if atoms[heavy_idx - 1] == "H":
                    continue
                ideal = _ideal_x_h_distance(atoms[heavy_idx - 1])
                hx,hy,hz = geom[h_idx - 1]
                xx,yx,zx = geom[heavy_idx - 1]
                actual = ((hx - xx)**2 + (hy - yx)**2 + (hz - zx)**2)**0.5
                if not _is_bad_h_distance(actual,ideal):
                    continue
                summary["bad_h_distances"] += 1
                msg = (f"Detected suspicious X-H distance in {path.name} [record {rec_idx}]: "
                       f"H{h_idx}-{atoms[heavy_idx - 1]}{heavy_idx} = {actual:.3f} Å "
                       f"(target ~{ideal:.3f} Å)")
                if fix:
                    new_xyz,used_vector = _reposition_bonded_hydrogen(geom[heavy_idx - 1],geom[h_idx - 1],ideal)
                    geom[h_idx - 1] = new_xyz
                    atom_lines[h_idx - 1] = _rewrite_v2000_atom_line_coords(atom_lines[h_idx - 1],*new_xyz)
                    summary["h_repositioned"] += 1
                    local_changed = True
                    log.warning("%s; repositioned H atom coordinates%s.",msg," using existing bond vector" if used_vector else " with axis fallback")
                else:
                    log.warning("%s; native/original mode keeps the input coordinates unchanged.",msg)
            if local_changed:
                changed_any = True
                new_record = list(record_lines)
                new_record[atom_start:atom_end] = atom_lines
                out_lines.extend(new_record)
            else:
                out_lines.extend(record_lines)
        else:
            raise ValueError(f"Unsupported SDF format while checking H distances: {path}")
        if had_terminator:
            out_lines.append("$$$$\n")
    if fix and changed_any:
        path.write_text("".join(out_lines))
    return summary

def validate_prepared_structure(structure_path: Path, logger: Optional[logging.Logger] = None, fix_sdf_h: bool = True, warn_h_connectivity_change: bool = False) -> Dict[str, object]:
    """Validate a prepared structure immediately before QM calculations."""
    path = Path(structure_path)
    result = {"connectivity": validate_no_isolated_atoms(path,logger = logger),"h_distance": None,"h_connectivity_warning": None}
    if path.suffix.lower() == ".sdf":
        result["h_distance"] = check_and_fix_sdf_h_distances(path,fix = bool(fix_sdf_h),logger = logger)
        h_summary = result.get("h_distance") or {}
        if h_summary.get("bad_h_distances", 0):
            if fix_sdf_h:
                _get_logger(logger).info("SDF H-distance sanitization summary for %s: bad = %d fixed = %d",path.name,h_summary.get("bad_h_distances",0),h_summary.get("h_repositioned",0))
            else:
                _get_logger(logger).info("SDF H-distance check summary for %s: bad = %d fixed = 0 (native/original mode)",path.name,h_summary.get("bad_h_distances",0))
        if warn_h_connectivity_change:
            result["h_connectivity_warning"] = warn_sdf_h_connectivity_changes(path,logger = logger)
            hc_summary = result.get("h_connectivity_warning") or {}
            if hc_summary.get("connectivity_shift_warnings",0):
                _get_logger(logger).info("Post-QM H-connectivity warning summary for %s: warnings = %d",path.name,hc_summary.get("connectivity_shift_warnings",0))
    return result

def _format_v2000_atom_line(x: float, y: float, z: float, elem: str) -> str:
    """Create a minimal V2000 atom line."""
    el = (str(elem).strip().upper() or "?")[:3]
    return (f"{x:10.4f}{y:10.4f}{z:10.4f} "
            f"{el:<3}"
            " 0  0  0  0  0  0  0  0  0  0  0  0\n")

def _insert_sdf_property_blocks(lines: List[str], props: Dict[str, object]) -> List[str]:
    """Insert SDF property blocks before $$$$ (or append if missing)."""
    if not props:
        return lines
    prop_blocks: List[str] = []
    for key,value in props.items():
        prop_blocks.append(f">  <{key}>\n{value}\n")
    prop_text = "".join(prop_blocks)
    out = list(lines)
    try:
        idx = next(i for i,ln in enumerate(out) if ln.strip() == "$$$$")
        out.insert(idx,prop_text)
    except StopIteration:
        if out and (not out[-1].endswith("\n")):
            out[-1] += "\n"
        out.append("\n")
        out.append(prop_text)
    return out

def write_template_sdf_with_xyz(template_sdf: Path, xyz_file: Path, out_sdf: Path, props: Optional[Dict[str, object]] = None, logger: Optional[logging.Logger] = None) -> None:
    """Write an SDF by taking connectivity from a V2000 template SDF and replacing coordinates with those from an XYZ file.
       Native/original extension:
         - if XYZ contains additional H atoms beyond the template,
           append those H atoms to the molblock
         - infer their X-H bonds geometrically using the same covalent-radius
           philosophy used by repair_native_sdf_h_connectivity()
       Constraints:
         - V2000 template only
         - heavy-atom sequence must be preserved
         - any extra XYZ atoms must be hydrogens"""
    log = _get_logger(logger)
    template_path = Path(template_sdf)
    xyz_path = Path(xyz_file)
    out_path = Path(out_sdf)
    props = props or {}
    if not template_path.exists():
        raise FileNotFoundError(f"Template SDF not found: {template_path}")
    if not xyz_path.exists():
        raise FileNotFoundError(f"XYZ file not found: {xyz_path}")
    xyz_atoms,xyz_geom = _read_xyz_atoms_and_geometry(xyz_path)
    if not xyz_atoms or not xyz_geom:
        raise ValueError(f"Empty XYZ geometry: {xyz_path}")
    lines = template_path.read_text().splitlines(True)
    if len(lines) < 5:
        raise ValueError(f"Template SDF too short / invalid: {template_path}")
    counts = lines[3]
    if "V2000" not in counts:
        raise ValueError(f"Template SDF must be V2000 for native/original augmentation: {template_path}")
    n_atoms,n_bonds = _parse_v2000_counts_line(counts)
    atom_start = 4
    atom_end = atom_start + n_atoms
    bond_start = atom_end
    bond_end = bond_start + n_bonds
    if len(lines) < bond_end:
        raise ValueError(f"Template SDF truncated before bond block: {template_path}")
    atom_lines = list(lines[atom_start:atom_end])
    bond_lines = list(lines[bond_start:bond_end])
    template_atoms = [_extract_v2000_element(ln) for ln in atom_lines]
    if len(xyz_atoms) < n_atoms:
        raise ValueError(f"XYZ has fewer atoms than template: template={n_atoms} xyz={len(xyz_atoms)} "
                         f"({template_path} vs {xyz_path})")
    #map template atoms onto the XYZ sequence while allowing extra XYZ hydrogens.
    template_to_xyz: List[int] = []
    xyz_idx = 0
    for t_elem in template_atoms:
        found = None
        while xyz_idx < len(xyz_atoms):
            x_elem = str(xyz_atoms[xyz_idx]).strip().upper()
            if x_elem == t_elem:
                found = xyz_idx
                xyz_idx += 1
                break
            if x_elem == "H":
                xyz_idx += 1
                continue
            raise ValueError("Heavy-atom order mismatch between template SDF and XYZ "
                             f"({template_path} vs {xyz_path}); encountered XYZ atom {x_elem} while matching {t_elem}")
        if found is None:
            raise ValueError(f"Failed to align template atom sequence to XYZ ({template_path} vs {xyz_path})")
        template_to_xyz.append(found)
    mapped_xyz = set(template_to_xyz)
    extra_h_indices = [i for i in range(len(xyz_atoms)) if i not in mapped_xyz]
    if any(str(xyz_atoms[i]).strip().upper() != "H" for i in extra_h_indices):
        bad = [xyz_atoms[i] for i in extra_h_indices if str(xyz_atoms[i]).strip().upper() != "H"]
        raise ValueError("Template/XYZ mismatch contains extra non-hydrogen atoms, which is unsupported: "
                         f"{bad} ({template_path} vs {xyz_path})")
    #rewrite coordinates for the template atoms.
    new_atom_lines: List[str] = []
    atoms_data: List[Tuple[str,float,float,float]] = []
    for i,old_line in enumerate(atom_lines):
        x,y,z = xyz_geom[template_to_xyz[i]]
        rest = old_line[30:] if len(old_line) > 30 else "\n"
        new_atom_lines.append(f"{x:10.4f}{y:10.4f}{z:10.4f}{rest}")
        #FMCM: is it needed to convert to float?
        atoms_data.append((template_atoms[i],float(x),float(y),float(z)))
    existing_pairs = set()
    for ln in bond_lines:
        a1 = int(ln[0:3])
        a2 = int(ln[3:6])
        existing_pairs.add(tuple(sorted((a1,a2))))
    new_bonds: List[str] = []
    r_h = AD._covalent_radius_angstrom("H") or 0.31
    for xyz_i in extra_h_indices:
        xh,yh,zh = xyz_geom[xyz_i]
        new_h_idx = len(atoms_data) + 1
        best_idx = None
        best_d2 = None
        for heavy_idx,(elem_x,xx,yx,zx) in enumerate(atoms_data,start = 1):
            if str(elem_x).strip().upper() == "H":
                continue
            r_x = AD._covalent_radius_angstrom(elem_x) or 0.76
            thr = 1.25*(r_h + r_x)
            d2 = (xh - xx)**2 + (yh - yx)**2 + (zh - zx)**2
            if d2 > thr*thr:
                continue
            if (best_d2 is None) or (d2 < best_d2):
                best_d2 = d2
                best_idx = heavy_idx
        if best_idx is None:
            raise ValueError("Failed to assign an added hydrogen from XYZ to a heavy atom in the template "
                            f"({template_path} vs {xyz_path})")
        pair = tuple(sorted((new_h_idx,best_idx)))
        if pair in existing_pairs:
            raise ValueError(f"Unexpected duplicate X-H bond while augmenting template: {pair} "
                             f"({template_path} vs {xyz_path})")
        new_atom_lines.append(_format_v2000_atom_line(xh,yh,zh,"H"))
        new_bonds.append(f"{new_h_idx:3d}{best_idx:3d}{1:3d}  0  0  0  0\n")
        existing_pairs.add(pair)
        atoms_data.append(("H",float(xh),float(yh),float(zh)))
    out_lines: List[str] = []
    out_lines.extend(lines[:4])
    out_lines[3] = _rebuild_v2000_counts_line(lines[3],n_atoms = len(new_atom_lines),n_bonds = len(bond_lines) + len(new_bonds))
    out_lines.extend(new_atom_lines)
    out_lines.extend(bond_lines)
    out_lines.extend(new_bonds)
    out_lines.extend(lines[bond_end:])
    out_lines = _insert_sdf_property_blocks(out_lines,props)
    out_path.parent.mkdir(parents = True,exist_ok = True)
    out_path.write_text("".join(out_lines))
    repair_prepared_sdf_h_connectivity(out_path,logger = log)
    validate_prepared_structure(out_path,logger = log,fix_sdf_h = True,warn_h_connectivity_change = True)
    if extra_h_indices:
        log.info("Augmented native/original SDF template %s with %d additional H atom(s) from %s and %d inferred bond(s).",
                  template_path.name,len(extra_h_indices),xyz_path.name,len(new_bonds))

class InputSanitizer:
    """Central input sanitization utilities:
         - strip explicit H atoms from input containers (unless --keep-input-h)
         - keep heavy-atom connectivity intact
         - repair selected ligand formal-charge edge cases in SDF inputs"""

    @staticmethod
    def strip_h_from_pdb_lines(pdb_lines: List[str]) -> Tuple[List[str], int]:
        return strip_h_from_pdb_lines(pdb_lines)

    @staticmethod
    def strip_h_from_sdf_file(sdf_in: str, sdf_out: str) -> Tuple[int, int]:
        return strip_h_from_sdf_file(sdf_in,sdf_out)

    @staticmethod
    def sanitize_ligand_file_for_protonation(ligand_path: str, output_folder: Path, keep_input_h: bool, logger: logging.Logger) -> str:
        return sanitize_ligand_file_for_protonation(ligand_path = ligand_path,output_folder = output_folder,keep_input_h = keep_input_h,logger = logger)

    @staticmethod
    def repair_ligand_sdf_formal_charges(sdf_path: Path, logger: Optional[logging.Logger] = None) -> Dict[str, int]:
        return repair_ligand_sdf_formal_charges(sdf_path = sdf_path,logger = logger)

    @staticmethod
    def get_sdf_total_formal_charge(sdf_path: Path, logger: Optional[logging.Logger] = None) -> Optional[int]:
        return get_sdf_total_formal_charge(sdf_path = sdf_path,logger = logger)

    @staticmethod
    def write_template_sdf_with_xyz(template_sdf: Path, xyz_file: Path, out_sdf: Path, props: Optional[Dict[str, object]] = None, logger: Optional[logging.Logger] = None) -> None:
        return write_template_sdf_with_xyz(template_sdf = template_sdf,xyz_file = xyz_file,out_sdf = out_sdf,props = props,logger = logger)

    @staticmethod
    def validate_no_isolated_atoms(structure_path: Path, logger: Optional[logging.Logger] = None) -> Dict[str, object]:
        return validate_no_isolated_atoms(structure_path = structure_path,logger = logger)

    @staticmethod
    def check_and_fix_sdf_h_distances(sdf_path: Path, fix: bool = True, logger: Optional[logging.Logger] = None) -> Dict[str, int]:
        return check_and_fix_sdf_h_distances(sdf_path = sdf_path,fix = fix,logger = logger)

    @staticmethod
    def repair_prepared_sdf_h_connectivity(sdf_path: Path, logger: Optional[logging.Logger] = None) -> None:
        return repair_prepared_sdf_h_connectivity(sdf_path = sdf_path,logger = logger)
  
    @staticmethod
    def warn_sdf_h_connectivity_changes(sdf_path: Path, logger: Optional[logging.Logger] = None) -> Dict[str, int]:
        return warn_sdf_h_connectivity_changes(sdf_path = sdf_path,logger = logger)

    @staticmethod
    def validate_prepared_structure(structure_path: Path, logger: Optional[logging.Logger] = None, fix_sdf_h: bool = True, warn_h_connectivity_change: bool = False) -> Dict[str, object]:
        return validate_prepared_structure(structure_path = structure_path,logger = logger,fix_sdf_h = fix_sdf_h,warn_h_connectivity_change = warn_h_connectivity_change)

__all__ = [
    "InputSanitizer",
    "get_sdf_total_formal_charge",
    "repair_ligand_sdf_formal_charges",
    "repair_native_sdf_h_connectivity",
    "sanitize_ligand_file_for_protonation",
    "strip_h_from_pdb_lines",
    "strip_h_from_sdf_file",
    "repair_prepared_sdf_h_connectivity",
    "validate_no_isolated_atoms",
    "check_and_fix_sdf_h_distances",
    "warn_sdf_h_connectivity_changes",
    "validate_prepared_structure"
]