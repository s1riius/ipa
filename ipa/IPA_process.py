# ------------------------------------------------------------
# Script to process IPA data
# ------------------------------------------------------------

from __future__ import annotations
import argparse
import logging
import re
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle

#config
STAT_LINE_RE = re.compile(r"^\s*(?P<key>[^:]+?)\s*:\s*(?P<val>[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*(?P<unit>.*)\s*$")

DEFAULT_STATS_KEYS = ["Deformation Energy",
                      "Distortion",
                      "Energy Chemical State",
                      "IPA Interaction Energy",
                      "IPA Binding Energy"]

#currently restricted to 50 simultaneous curves in the interaction plot.
MAX_CURVES_FOR_INTERACTION_PLOT = 50

#CSS palette with required starting order:
CSS_PALETTE = ["black","blue","red","green","darkviolet","darkorange","deeppink",
               "brown","darkcyan","gold","magenta","navy","olive","teal",
               "crimson","indigo","limegreen","orchid","peru","royalblue",
               "salmon","sienna","slateblue","turquoise","tomato",
               "mediumvioletred","darkgreen","darkred","darkblue","firebrick",
               "forestgreen","steelblue","mediumseagreen","mediumorchid",
               "cadetblue","chocolate","coral","cornflowerblue","darkgoldenrod",
               "darkmagenta","darkslateblue","dodgerblue","hotpink","khaki",
               "lightseagreen","maroon","mediumblue","mediumspringgreen",
               "midnightblue"]

#IO file parsing
def parse_io_file(io_path: Path, default_outdir: Path) -> Tuple[List[Path], Path, str, bool]:
    """Input file lines (non-empty, non-comment) can be:
      i <path>          -> input system folder (each i-line corresponds to one system)
      o <path>          -> output directory (last one wins)
      env intersect     -> environment mode (default: intersect); the data is presented for the intersection of residues in the analyses, i.e., for the residues present in all analyses
      env unite         -> environment mode (unite); the data is presented for the union of residues in the analyses, i.e., for the residues appearing in at least one analysis
      first             -> choose FIRST chemical state (override "best", which is the default)
    Also accepted shorthand lines:
      intersect
      unite
    Returns: (input_folders, outdir, env_mode, use_first_state)"""
    in_folders: List[Path] = []
    outdir = default_outdir
    env_mode = "intersect"
    use_first_state = False
    for raw in io_path.read_text(encoding = "utf-8",errors = "replace").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        low = line.lower().strip()
        if low == "first":
            use_first_state = True
            continue
        parts = line.split(maxsplit = 1)
        tag0 = parts[0].strip().lower()
        rest_raw = parts[1].strip() if len(parts) == 2 else ""
        rest = rest_raw.lower()
        if tag0 in ("intersect","unite"):
            env_mode = tag0
            continue
        if tag0 == "env":
            if rest not in ("intersect","unite"):
                raise ValueError(f"Bad env line (expected 'env intersect' or 'env unite'): {raw!r}")
            env_mode = rest
            continue
        if tag0 == "i":
            if not rest_raw:
                raise ValueError(f"Bad line (expected 'i <path>'): {raw!r}")
            in_folders.append(Path(rest_raw).expanduser())
            continue
        if tag0 == "o":
            if not rest_raw:
                raise ValueError(f"Bad line (expected 'o <path>'): {raw!r}")
            outdir = Path(rest_raw).expanduser()
            continue
        raise ValueError(f"Unknown line in input file: {raw!r}")
    return in_folders, outdir, env_mode, use_first_state

#stats.txt parsing (multi-state)
@dataclass
class ChemicalState:
    idx: int
    label: str
    metrics: Dict[str,float]
    units: Dict[str,str]

def parse_stats_txt_multi(stats_path: Path) -> List[ChemicalState]:
    """stats.txt contains multiple chemical states (blocks).
    Each block starts with a label line ending in ':' and then key/value lines.
    Returns chemical states in file order, indexed from 1 ... N."""
    lines = stats_path.read_text(encoding = "utf-8",errors = "replace").splitlines()
    states: List[ChemicalState] = []
    cur_label: Optional[str] = None
    cur_metrics: Dict[str,float] = {}
    cur_units: Dict[str,str] = {}

    def flush() -> None:
        nonlocal cur_label,cur_metrics,cur_units
        if cur_label is None: return
        states.append(ChemicalState(idx = len(states) + 1,label = cur_label,metrics = cur_metrics,units = cur_units))
        cur_label = None
        cur_metrics = {}
        cur_units = {}

    for raw in lines:
        s = raw.strip()
        if not s: continue
        #block header
        if s.endswith(":"):
            flush()
            cur_label = s[:-1].strip() or f"state_{len(states) + 1}"
            continue
        m = STAT_LINE_RE.match(raw)
        if m and cur_label is not None:
            key = m.group("key").strip()
            val = float(m.group("val"))
            unit = (m.group("unit") or "").strip()
            cur_metrics[key] = val
            cur_units[key] = unit
    flush()
    return states

#interaction_energies.csv parsing
def parse_interaction_csv_multi(csv_path: Path) -> Dict[int, pd.DataFrame]:
    """Parses multi-block ';' files like:
         label:
         header;...
         data...
         label2:
         header;...
         data...
       If no blocks detected, parses as a single table (comma first, then ';') and returns {1: df}."""
    text = csv_path.read_text(encoding = "utf-8",errors = "replace")
    lines = text.splitlines()
    block_starts: List[int] = []
    for i in range(len(lines) - 1):
        if lines[i].strip().endswith(":") and (";" in lines[i + 1]):
            block_starts.append(i)
    if not block_starts:
        try:
            df = pd.read_csv(StringIO(text))
            return {1:df}
        except Exception:
            df = pd.read_csv(StringIO(text),sep = ";")
            return {1:df}
    out: Dict[int,pd.DataFrame] = {}
    for b_idx,start in enumerate(block_starts):
        end = block_starts[b_idx + 1] if b_idx + 1 < len(block_starts) else len(lines)
        header = lines[start + 1].strip()
        data_lines = [ln for ln in lines[start + 2 : end] if ln.strip()]
        if not data_lines: continue
        buf = StringIO()
        buf.write(header + "\n")
        buf.write("\n".join(data_lines) + "\n")
        buf.seek(0)
        df = pd.read_csv(buf,sep = ";")
        out[len(out) + 1] = df
    if not out:
        df = pd.read_csv(csv_path,sep = ";",skiprows = 1)
        return {1: df}
    return out

#main processing
@dataclass
class FolderResult:
    system_id: str
    folder: str
    interactions_by_state: Dict[int,pd.DataFrame]
    interactions_raw: pd.DataFrame
    chemical_states: List[ChemicalState]

def process_folder(folder: Path,system_id: str) -> Optional[FolderResult]:
    inter_path = folder / "interaction_energies.csv"
    stats_path = folder / "stats.txt"
    if not inter_path.exists():
        logging.warning("Missing %s",inter_path)
        return None
    if not stats_path.exists():
        logging.warning("Missing %s",stats_path)
        return None
    chemical_states = parse_stats_txt_multi(stats_path)
    if not chemical_states:
        logging.warning("No chemical states parsed from %s",stats_path)
        return None
    interactions_by_state = parse_interaction_csv_multi(inter_path)
    if not interactions_by_state:
        logging.warning("No interaction data parsed from %s",inter_path)
        return None
    #concatenate into raw (and carry __chem_state_idx when multi-block)
    if len(interactions_by_state) > 1:
        parts = []
        for k,df in interactions_by_state.items():
            d2 = df.copy()
            d2["__chem_state_idx"] = k
            parts.append(d2)
        interactions_raw = pd.concat(parts,ignore_index = True)
    else:
        interactions_raw = list(interactions_by_state.values())[0].copy()
    #add identifiers to each per-state df
    for k,df in list(interactions_by_state.items()):
        df2 = df.copy()
        df2["system_id"] = system_id
        df2["folder"] = folder.name
        interactions_by_state[k] = df2
    interactions_raw["system_id"] = system_id
    interactions_raw["folder"] = folder.name
    return FolderResult(system_id = system_id,folder = folder.name,interactions_by_state = interactions_by_state,interactions_raw = interactions_raw,chemical_states = chemical_states)

#plot helpers
def _get_colors(n: int) -> List[str]:
    if n > MAX_CURVES_FOR_INTERACTION_PLOT:
        raise ValueError(f"Too many curves ({n}). Current restriction is {MAX_CURVES_FOR_INTERACTION_PLOT}.")
    if n > len(CSS_PALETTE):
        raise ValueError(f"CSS palette only has {len(CSS_PALETTE)} colors; requested {n}. Extend CSS_PALETTE.")
    return CSS_PALETTE[:n]

def _infer_residue_column_fallback(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    for c in cols:
        if "res" in str(c).lower():
            return c
    return cols[0]

def _get_residue_labels(df: pd.DataFrame) -> pd.Series:
    if "res+resn" in df.columns: return df["res+resn"].astype(str)
    if "res" in df.columns and "resn" in df.columns: return (df["res"].astype(str) + df["resn"].astype(str)).astype(str)
    col = _infer_residue_column_fallback(df)
    logging.warning("Could not find 'res+resn' (or 'res'+'resn'). Falling back to column: %s", col)
    return df[col].astype(str)

def _default_prot_state_for_res(resname3: str) -> Optional[str]:
    r = resname3.upper()
    if r == "LYS": return "+"
    elif r in ("GLU","ASP"): return "-"
    elif r == "ARG": return "+"
    elif r == "HIS": return None
    return None

def _get_base_and_prot_state(df: pd.DataFrame) -> Tuple[pd.Series, pd.Series]:
    base = _get_residue_labels(df).astype(str)
    state_cols = ["protonation","prot","state","charge","chg","q","res_state"]
    col = next((c for c in state_cols if c in df.columns),None)
    if col is None:
        return base,pd.Series(["__nostate__"]*len(df),index = df.index)
    st = df[col].astype(str).str.strip()
    st = st.replace({"nan": "","None": "","NA": "","N/A": ""})
    st = st.where(st.isin(["+","-","0"]),other = "__nostate__")
    return base,st

def _get_flipped_flags(df: pd.DataFrame) -> pd.Series:
    flip_cols = ["flip","flipped","is_flipped","flp"]
    col = next((c for c in flip_cols if c in df.columns),None)
    if col is None:
        return pd.Series([False]*len(df),index = df.index)
    s = df[col]
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False)
    s2 = s.astype(str).str.strip().str.lower()
    true_set = {"1","true","t","yes","y","flipped","flip"}
    return s2.isin(true_set)

def _choose_best_state_idx_from_stats(r: FolderResult, use_first_state: bool) -> int:
    """Chooses:
      - state 1 if 'first' was set in input file
      - otherwise, the state with the LOWEST (most negative) IPA Binding Energy
        If missing, falls back to 1."""
    if use_first_state:
        return 1
    best_idx = 1
    best_val: Optional[float] = None
    for st in r.chemical_states:
        v = st.metrics.get("IPA Binding Energy",None)
        if v is None: continue
        if best_val is None or v < best_val:
            best_val = float(v)
            best_idx = st.idx
    return best_idx

def _find_chem_state_col(df: pd.DataFrame) -> Optional[str]:
    candidates = ["__chem_state_idx",
                  "chemical_state","chem_state","chemical_state_index","state_index",
                  "ligand_state","lig_state","chemstate","cs","state"]
    for c in candidates:
        if c in df.columns:
            return c
    return None

def _filter_df_to_state(df: pd.DataFrame, state_idx: int) -> pd.DataFrame:
    col = _find_chem_state_col(df)
    if col is None: return df
    sn = pd.to_numeric(df[col],errors = "coerce")
    if sn.notna().any():
        out = df.loc[sn == state_idx].copy()
        return out if not out.empty else df
    ss = df[col].astype(str).str.strip().str.lower()
    want = str(state_idx)
    mask = (ss == want) | ss.str.contains(rf"\b{re.escape(want)}\b",regex = True)
    out = df.loc[mask].copy()
    return out if not out.empty else df

def _split_df_by_state(df: pd.DataFrame) -> Dict[int, pd.DataFrame]:
    col = _find_chem_state_col(df)
    if col is None: return {1: df}
    sn = pd.to_numeric(df[col],errors = "coerce")
    if sn.notna().any():
        out: Dict[int,pd.DataFrame] = {}
        for sidx in sorted({int(x) for x in sn.dropna().unique().tolist()}):
            part = df.loc[sn == sidx].copy()
            if not part.empty:
                out[sidx] = part
        return out if out else {1: df}
    ss = df[col].astype(str).str.strip()
    extracted = ss.str.extract(r"(\d+)")[0]
    sn2 = pd.to_numeric(extracted,errors = "coerce")
    if sn2.notna().any():
        out: Dict[int,pd.DataFrame] = {}
        for sidx in sorted({int(x) for x in sn2.dropna().unique().tolist()}):
            part = df.loc[sn2 == sidx].copy()
            if not part.empty:
                out[sidx] = part
        return out if out else {1: df}
    return {1: df}

#stats plots: stacked by system, not spread across x-axis
def make_stats_per_key_plots_multi(results: List[FolderResult], outdir: Path, keys: List[str]) -> None:
    """Chemical states are plotted vertically per system (colored by state index).
       Writes a CSV per stats plot with columns: system_id, state_idx, value."""
    if not results:
        logging.warning("No results for stats plots.")
        return
    system_order = [r.system_id for r in results]
    #determine maximum number of states for consistent coloring
    max_state = 1
    for r in results:
        max_state = max(max_state,max((st.idx for st in r.chemical_states),default = 1))
    state_colors = _get_colors(min(max_state,len(CSS_PALETTE)))
    st2c = {i + 1: state_colors[i] for i in range(len(state_colors))}
    #fixed x positions per system
    x_pos = {sid: i + 1 for i, sid in enumerate(system_order)}
    x_ticks = [x_pos[sid] for sid in system_order]
    for key in keys:
        rows: List[Dict[str,object]] = []
        xs: List[float] = []
        ys: List[float] = []
        cs: List[str] = []
        labels_for_legend: Dict[int,str] = {}
        for r in results:
            for st in r.chemical_states:
                if key not in st.metrics: continue
                val = float(st.metrics[key])
                rows.append({"system_id": r.system_id,"state_idx": st.idx,"value": val})
                xs.append(float(x_pos[r.system_id]))
                ys.append(val)
                cs.append(st2c.get(st.idx,CSS_PALETTE[-1]))
                labels_for_legend[st.idx] = f"state {st.idx}"
        if not rows:
            logging.warning("Stats key not found anywhere (skipping plot): %s",key)
            continue
        #write CSV for this plot
        safe = re.sub(r"[^A-Za-z0-9_.-]+","_",key).strip("_")
        df_out = pd.DataFrame(rows).sort_values(["system_id","state_idx"])
        df_out.to_csv(outdir / f"stats_{safe}.csv",index = False)
        #plot (stacked points per system; color by state)
        fig,ax = plt.subplots(figsize = (max(10,len(system_order)),4.8))
        ax.scatter(xs,ys,c = cs,s = 55)
        #make a clean legend: one handle per state_idx that appears
        handles = []
        shown_states = sorted(set(df_out["state_idx"].tolist()))
        for sidx in shown_states:
            handles.append(plt.Line2D([0],[0],marker = "o",linestyle = "",color = st2c.get(sidx,CSS_PALETTE[-1]),label = f"state {sidx}"))
        ax.legend(handles = handles,ncol = min(6,len(handles)),fontsize = "small")
        ax.set_title(key)
        ax.set_xlabel("System",fontsize = 14,fontweight = "bold")
        ax.set_ylabel("Energy (kcal/mol)",fontsize = 14,fontweight = "bold")
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(system_order,rotation = 0)
        fig.tight_layout()
        fig.savefig(outdir / f"stats_{safe}.png",dpi = 200)
        plt.close(fig)

#interaction plot
def make_interactions_plot(results: List[FolderResult], out_png: Path, env_mode: str, use_first_state: bool, energy_col: str = "Eint(kcal/mol)") -> None:
    """Behavior:
         - If 1 system folder: plot ALL chemical states for that folder.
         - If >1 system folder: plot ONLY ONE state per system:
             * first if 'first' is in io file
             * else best (lowest IPA Binding Energy) from stats.txt
       env_mode:
         - intersect: only plot residues present in ALL curves (no fake zeros)
         - unite: plot union; missing residues get 0 interaction energy
       Also writes a CSV for the interaction plot (same name as png but .csv),
       with one row per residue and one column per curve label."""
    if not results: raise ValueError("No results to plot.")
    curves: List[Tuple[str,pd.DataFrame]] = []
    #decide curves
    if len(results) == 1:
        r = results[0]
        #plot all states available
        if len(r.interactions_by_state) > 1:
            by_state = dict(sorted(r.interactions_by_state.items(),key = lambda kv: kv[0]))
        else:
            by_state = _split_df_by_state(r.interactions_raw)
        state_indices = sorted(by_state.keys())
        if len(state_indices) > MAX_CURVES_FOR_INTERACTION_PLOT:
            raise ValueError(f"Too many chemical states ({len(state_indices)}) to plot at once; "
                             f"limit is {MAX_CURVES_FOR_INTERACTION_PLOT}.")
        for sidx in state_indices:
            curves.append((f"{r.system_id}({sidx})",by_state[sidx]))
    else:
        if len(results) > MAX_CURVES_FOR_INTERACTION_PLOT:
            raise ValueError(f"Too many systems ({len(results)}) to plot at once; limit is {MAX_CURVES_FOR_INTERACTION_PLOT}.")
        for r in results:
            chosen = _choose_best_state_idx_from_stats(r,use_first_state = use_first_state)
            #prefer explicit per-state df if we truly have multi-block
            if chosen in r.interactions_by_state and len(r.interactions_by_state) > 1:
                df_state = r.interactions_by_state[chosen]
            else:
                #otherwise filter from raw by state column if present; if not present, raw is single-state
                df_state = _filter_df_to_state(r.interactions_raw,chosen)
            curves.append((f"{r.system_id}({chosen})",df_state))
    colors = _get_colors(len(curves))
    #per-curve maps and PRESENCE sets
    per_curve_states: Dict[str,Dict[str,Dict[str,float]]] = {}
    per_curve_flipped: Dict[str,Dict[str,float]] = {}
    presence_sets: List[set] = []
    for (label,df) in curves:
        if energy_col not in df.columns:
            raise ValueError(f"Column {energy_col!r} not found for curve {label}")
        base,prot = _get_base_and_prot_state(df)
        flipped = _get_flipped_flags(df)
        energies = pd.to_numeric(df[energy_col],errors = "coerce").fillna(0.0)
        present = set(str(b).strip() for b in base.tolist() if str(b).strip() and str(b).strip().lower() != "nan")
        presence_sets.append(present)
        states_map: Dict[str,Dict[str,float]] = {}
        flipped_map: Dict[str,float] = {}
        for b,st,e,fl in zip(base.tolist(),prot.tolist(),energies.tolist(),flipped.tolist()):
            b = str(b).strip()
            if not b or b.lower() == "nan": continue
            st = str(st).strip() if st else "__nostate__"
            if st not in ["+","-","0"]:
                st = "__nostate__"
            states_map.setdefault(b,{})
            states_map[b][st] = float(states_map[b].get(st,0.0) + float(e))
            if fl:
                flipped_map[b] = float(flipped_map.get(b,0.0) + float(e))
        per_curve_states[label] = states_map
        per_curve_flipped[label] = flipped_map
    #determine residues to plot
    if env_mode == "intersect":
        bases_final = sorted(set.intersection(*presence_sets)) if presence_sets else []
    elif env_mode == "unite":
        bases_final = sorted(set.union(*presence_sets)) if presence_sets else []
    else:
        raise ValueError(f"Unknown env_mode: {env_mode!r}")
    if not bases_final:
        raise ValueError(f"No residues to plot after applying env mode: {env_mode}")
    x_numeric = list(range(1,len(bases_final) + 1))
    fig,ax = plt.subplots(figsize = (max(10,len(bases_final)*0.15),6))
    #CSV table we will fill: one row per residue, columns per curve
    csv_table = pd.DataFrame({"residue": bases_final})
    #precompute y-range including detached & flipped so nothing spills out
    all_y: List[float] = []
    cache: List[Dict[str,object]] = []
    for (label,_),color in zip(curves,colors):
        m = per_curve_states[label]
        mf = per_curve_flipped[label]
        y_line: List[float] = []
        det0_x: List[int] = []
        det0_y: List[float] = []
        detp_x: List[int] = []
        detp_y: List[float] = []
        detf_x: List[int] = []
        detf_y: List[float] = []
        for i,base in enumerate(bases_final,start = 1):
            states = m.get(base,{})
            #unite: missing -> 0
            #intersect: missing -> NaN (shouldn't happen if intersection is correct, but safe)
            if not states:
                if env_mode == "unite":
                    y_line.append(0.0)
                    all_y.append(0.0)
                else:
                    y_line.append(float("nan"))
                continue
            resname3 = base[:3].upper()
            default = _default_prot_state_for_res(resname3)
            #choose connected-line protonation state
            if resname3 == "HIS":
                if "0" in states:
                    line_state = "0"
                elif "__nostate__" in states:
                    line_state = "__nostate__"
                else:
                    line_state = list(states.keys())[0]
            else:
                if default is not None and default in states:
                    line_state = default
                elif "__nostate__" in states:
                    line_state = "__nostate__"
                else:
                    pref = ["+","-","0"]
                    line_state = next((p for p in pref if p in states),list(states.keys())[0])
            yv = float(states.get(line_state,0.0))
            y_line.append(yv)
            all_y.append(yv)
            #detached 0-state: hollow circle
            if resname3 == "LYS" and "0" in states and line_state != "0":
                det0_x.append(i)
                det0_y.append(float(states["0"]))
                all_y.append(float(states["0"]))
            if resname3 in ("GLU","ASP") and "0" in states and line_state != "0":
                det0_x.append(i)
                det0_y.append(float(states["0"]))
                all_y.append(float(states["0"]))
            #detached HIS '+'
            if resname3 == "HIS" and "+" in states and line_state != "+":
                detp_x.append(i)
                detp_y.append(float(states["+"]))
                all_y.append(float(states["+"]))
            #flipped
            if base in mf:
                detf_x.append(i)
                detf_y.append(float(mf[base]))
                all_y.append(float(mf[base]))
        cache.append({"label": label,
                      "color": color,
                      "y_line": y_line,
                      "det0_x": det0_x,
                      "det0_y": det0_y,
                      "detp_x": detp_x,
                      "detp_y": detp_y,
                      "detf_x": detf_x,
                      "detf_y": detf_y})
        #add curve to CSV table
        csv_table[label] = y_line
    #write interaction CSV (same basename as png)
    out_csv = out_png.with_suffix(".csv")
    csv_table.to_csv(out_csv,index = False)
    if all_y:
        y_min = min(all_y)
        y_max = max(all_y)
        span = y_max - y_min
        pad = max(0.5,0.05*span)
        ax.set_ylim(y_min - pad,y_max + pad)
    #plot curves
    for item in cache:
        label = str(item["label"])
        color = str(item["color"])
        y_line = item["y_line"]
        ax.plot(x_numeric,y_line,label = label,color = color,marker = "o",markersize = 6,linewidth = 1.5)
        det0_x = item["det0_x"]
        det0_y = item["det0_y"]
        if det0_x:
            ax.scatter(det0_x,det0_y,marker = "o",facecolors = "none",edgecolors = color,linewidths = 1.5,s = 45,zorder = 6)
        detp_x = item["detp_x"]
        detp_y = item["detp_y"]
        if detp_x:
            ax.scatter(detp_x,detp_y,marker = "+",color = color,linewidths = 1.8,s = 70,zorder = 7)
        detf_x = item["detf_x"]
        detf_y = item["detf_y"]
        if detf_x:
            for xi,yi in zip(detf_x,detf_y):
                ax.text(xi,yi,"f",color = color,fontsize = 10,fontweight = "bold",ha = "center",va = "center",zorder = 8,clip_on = True)
    ax.set_xlabel("Residue",fontsize = 14,fontweight = "bold")
    ax.set_ylabel("Energy (kcal/mol)",fontsize = 14,fontweight = "bold")
    ax.set_xticks(x_numeric)
    ax.set_xticklabels(bases_final,rotation = 90)
    ax.legend(ncol = min(5,len(curves)),fontsize = "small")
    fig.tight_layout()
    fig.savefig(out_png,dpi = 200)
    plt.close(fig)

def make_interactions_heatmap(results: List[FolderResult], out_png: Path, env_mode: str, use_first_state: bool, energy_col: str = "Eint(kcal/mol)") -> None:
    """Heatmap:
         - rows: residues (res+resn)
         - columns: curves (systems or chem states, using the same selection logic as make_interactions_plot)
         - values: chosen line-state interaction energy per residue (same logic as line plot)
         - env_mode:
             * intersect -> only residues present in all curves
             * unite     -> union; missing residues -> 0
       Also writes CSV next to PNG."""
    if not results:
        raise ValueError("No results to plot.")
    #decide curves (same as make_interactions_plot)
    curves: List[Tuple[str,pd.DataFrame]] = []
    if len(results) == 1:
        r = results[0]
        if len(r.interactions_by_state) > 1:
            by_state = dict(sorted(r.interactions_by_state.items(),key = lambda kv: kv[0]))
        else:
            by_state = _split_df_by_state(r.interactions_raw)
        for sidx in sorted(by_state.keys()):
            curves.append((f"{r.system_id}({sidx})",by_state[sidx]))
    else:
        for r in results:
            chosen = _choose_best_state_idx_from_stats(r,use_first_state = use_first_state)
            if chosen in r.interactions_by_state and len(r.interactions_by_state) > 1:
                df_state = r.interactions_by_state[chosen]
            else:
                df_state = _filter_df_to_state(r.interactions_raw,chosen)
            curves.append((f"{r.system_id}({chosen})",df_state))
    if len(curves) > MAX_CURVES_FOR_INTERACTION_PLOT:
        raise ValueError(f"Too many curves ({len(curves)}) for heatmap; limit is {MAX_CURVES_FOR_INTERACTION_PLOT}.")
    #build per-curve residue->states map + presence sets
    per_curve_states: Dict[str,Dict[str,Dict[str,float]]] = {}
    presence_sets: List[set] = []
    for (label,df) in curves:
        if energy_col not in df.columns:
            raise ValueError(f"Column {energy_col!r} not found for curve {label}")
        base,prot = _get_base_and_prot_state(df)
        energies = pd.to_numeric(df[energy_col],errors = "coerce").fillna(0.0)
        present = set(str(b).strip()
                      for b in base.tolist()
                      if str(b).strip() and str(b).strip().lower() != "nan")
        presence_sets.append(present)
        states_map: Dict[str,Dict[str,float]] = {}
        for b,st,e in zip(base.tolist(),prot.tolist(),energies.tolist()):
            b = str(b).strip()
            if not b or b.lower() == "nan": continue
            st = str(st).strip() if st else "__nostate__"
            if st not in ["+","-","0"]:
                st = "__nostate__"
            states_map.setdefault(b,{})
            states_map[b][st] = float(states_map[b].get(st,0.0) + float(e))
        per_curve_states[label] = states_map
    #determine residues (rows) according to env_mode
    if env_mode == "intersect":
        residues = sorted(set.intersection(*presence_sets)) if presence_sets else []
    elif env_mode == "unite":
        residues = sorted(set.union(*presence_sets)) if presence_sets else []
    else:
        raise ValueError(f"Unknown env_mode: {env_mode!r}")
    if not residues:
        raise ValueError(f"No residues to plot after applying env mode: {env_mode}")
    #build matrix: rows = residues, cols = curves
    col_labels = [lab for (lab,_) in curves]
    mat = np.zeros((len(residues),len(col_labels)),dtype = float)
    for j,label in enumerate(col_labels):
        m = per_curve_states[label]
        for i,base in enumerate(residues):
            states = m.get(base,{})
            if not states:
                #unite -> 0 already; intersect shouldn't happen, but keep 0
                mat[i,j] = 0.0
                continue
            resname3 = base[:3].upper()
            default = _default_prot_state_for_res(resname3)
            #choose connected-line protonation state (same rule as line plot)
            if resname3 == "HIS":
                if "0" in states: line_state = "0"
                elif "__nostate__" in states: line_state = "__nostate__"
                else: line_state = list(states.keys())[0]
            else:
                if default is not None and default in states: line_state = default
                elif "__nostate__" in states: line_state = "__nostate__"
                else:
                    pref = ["+","-","0"]
                    line_state = next((p for p in pref if p in states),list(states.keys())[0])
            mat[i,j] = float(states.get(line_state,0.0))
    #write CSV (wide matrix: residue + one column per curve)
    df_csv = pd.DataFrame(mat,columns = col_labels)
    df_csv.insert(0,"residue",residues)
    out_csv = out_png.with_suffix(".csv")
    df_csv.to_csv(out_csv,index = False)
    #plot heatmap
    fig,ax = plt.subplots(figsize = (max(10,0.45*len(col_labels)),max(6,0.22*len(residues))))
    #make color scale symmetric around 0
    max_abs = np.nanmax(np.abs(mat))
    norm = mcolors.TwoSlopeNorm(vmin = -max_abs,vcenter = 0.0,vmax = max_abs)
    im = ax.imshow(mat,
                   aspect = "auto",
                   interpolation = "nearest",
                   cmap = "bwr",                    #blue–white–red
                   norm = norm)                     #centered at 0
    ax.set_title("Interaction energy heatmap")
    ax.set_xlabel("Inhibitor / curve",fontsize = 14,fontweight = "bold")
    ax.set_ylabel("Residue",fontsize = 14,fontweight = "bold")
    ax.set_xticks(range(len(col_labels)))
    ax.set_xticklabels(col_labels,rotation = 45,ha = "right")
    ax.set_yticks(range(len(residues)))
    ax.set_yticklabels(residues)
    cbar = fig.colorbar(im,ax = ax)
    cbar.set_label("Energy (kcal/mol)",fontsize = 12,fontweight = "bold")
    fig.tight_layout()
    fig.savefig(out_png,dpi = 200)
    plt.close(fig)

def compute_heatmap_row_stats(df_heat: pd.DataFrame) -> pd.DataFrame:
    """df_heat must be the heatmap matrix in WIDE format:
         - first column: 'residue'
         - other columns: curves/inhibitors
         - cell values: energies (floats)
       Returns a DataFrame with:
         residue, n, mean, median, min, max, std"""
    if "residue" not in df_heat.columns:
        raise ValueError("Heatmap DataFrame must contain a 'residue' column.")
    value_cols = [c for c in df_heat.columns if c != "residue"]
    if not value_cols:
        raise ValueError("Heatmap DataFrame must contain at least one curve/inhibitor column.")
    vals = df_heat[value_cols].apply(pd.to_numeric,errors = "coerce")
    stats = pd.DataFrame({"residue": df_heat["residue"].astype(str),
                          "n": vals.notna().sum(axis = 1).astype(int),
                          "mean": vals.mean(axis = 1,skipna = True),
                          "median": vals.median(axis = 1,skipna = True),
                          "min": vals.min(axis = 1,skipna = True),
                          "max": vals.max(axis = 1,skipna = True),
                          "std": vals.std(axis = 1,ddof = 1,skipna = True)})            #sample std; change to ddof = 0 if you want population std
    return stats

def write_heatmap_stats(heatmap_csv: Path, out_csv: Path) -> pd.DataFrame:
    """Reads the heatmap CSV (wide), computes row-wise stats, writes them to out_csv.
       Returns the stats DataFrame."""
    df_heat = pd.read_csv(heatmap_csv)
    stats = compute_heatmap_row_stats(df_heat)
    stats.to_csv(out_csv,index = False)
    return stats

def plot_residue_candlestick(stats_df: pd.DataFrame, out_png: Path, title: str = "Residue interaction statistics") -> None:
    """Candlestick-style plot per residue with SYSTEMATIC styling:
         - rectangle spans min..max
             * outline: dark blue
             * fill: light grey
         - median: black horizontal line inside rectangle
         - mean: black 'x' marker
         - std: black whiskers for mean ± std (clipped to [min, max]) + black caps
       Expects stats_df columns:
         residue, mean, median, min, max, std"""
    required = {"residue","mean","median","min","max","std"}
    missing = required - set(stats_df.columns)
    if missing: raise ValueError(f"stats_df missing columns: {sorted(missing)}")
    residues = stats_df["residue"].astype(str).tolist()
    x = np.arange(1,len(residues) + 1)
    vmin = pd.to_numeric(stats_df["min"],errors = "coerce").to_numpy()
    vmax = pd.to_numeric(stats_df["max"],errors = "coerce").to_numpy()
    vmed = pd.to_numeric(stats_df["median"],errors = "coerce").to_numpy()
    vmean = pd.to_numeric(stats_df["mean"],errors = "coerce").to_numpy()
    vstd = pd.to_numeric(stats_df["std"],errors = "coerce").fillna(0.0).to_numpy()
    #whiskers for mean ± std, clipped to [min,max]
    wlow = np.maximum(vmin,vmean - vstd)
    whigh = np.minimum(vmax,vmean + vstd)
    #styling knobs
    rect_edge = "#27408B"
    rect_face = "0.90"            #light grey (string grayscale)
    black = "black"
    box_width = 0.65
    cap_half = 0.18
    fig,ax = plt.subplots(figsize = (max(12,0.35*len(residues)),6))
    all_y = []
    for xi,lo,hi,med,mu,wl,wh in zip(x,vmin,vmax,vmed,vmean,wlow,whigh):
        if np.isnan(lo) or np.isnan(hi):
            continue
        #rectangle min..max (filled lwith "rect_face", and with "rect_edge" outline)
        rect = Rectangle((xi - box_width/2,lo),box_width,hi - lo,facecolor = rect_face,edgecolor = rect_edge,linewidth = 1.4)
        ax.add_patch(rect)
        #median line (black)
        if not np.isnan(med):
            ax.hlines(med,xi - box_width/2,xi + box_width/2,colors = black,linewidth = 1.8,zorder = 5)
        #mean as black 'x'
        if not np.isnan(mu):
            ax.plot([xi],[mu],marker = "x",markersize = 7,linestyle = "None",color = black,markeredgewidth = 1.6,zorder = 6)
        #std whiskers (black), clipped to [min,max], with caps
        if not (np.isnan(wl) or np.isnan(wh)):
            ax.vlines(xi,wl,wh,colors = black,linewidth = 1.4,zorder = 6)
            ax.hlines([wl,wh],xi - cap_half,xi + cap_half,colors = black,linewidth = 1.4,zorder = 6)
        all_y.extend([lo,hi,med,mu,wl,wh])
    #y limits with padding
    if all_y:
        y_min = np.nanmin(all_y)
        y_max = np.nanmax(all_y)
        span = y_max - y_min
        pad = max(0.5,0.05*span)
        ax.set_ylim(y_min - pad,y_max + pad)
    ax.set_title(title)
    ax.set_xlabel("Residue",fontsize = 14,fontweight = "bold")
    ax.set_ylabel("Energy (kcal/mol)",fontsize = 14,fontweight = "bold")
    ax.set_xticks(x)
    ax.set_xticklabels(residues,rotation = 90)
    ax.grid(True,axis = "y",linewidth = 0.6,alpha = 0.4)
    fig.tight_layout()
    fig.savefig(out_png,dpi = 200)
    plt.close(fig)

#CLI
def main() -> int:
    ap = argparse.ArgumentParser(description = "Read analysis folders and plot results.")
    ap.add_argument("io_txt",type = Path,help = "input file with lines: 'i <folder>', optional 'o <outdir>', optional 'env intersect|unite', optional 'first'")
    ap.add_argument("--default-outdir",type = Path,default = Path("out"),help = "default output directory")
    ap.add_argument("--stats-keys",nargs = "*",default = DEFAULT_STATS_KEYS,help = "stats.txt keys to plot")
    ap.add_argument("--verbose",action = "store_true",help = "verbose logging")
    args = ap.parse_args()
    logging.basicConfig(level = logging.DEBUG if args.verbose else logging.INFO,format = "%(levelname)s: %(message)s")
    in_folders,outdir,env_mode,use_first_state = parse_io_file(args.io_txt,args.default_outdir)
    if not in_folders:
        logging.error("No input folders found (no 'i <path>' lines) in %s",args.io_txt)
        return 2
    outdir.mkdir(parents = True,exist_ok = True)
    logging.info("Output directory: %s",outdir)
    logging.info("Interaction environment mode: %s",env_mode)
    logging.info("Interaction state selection: %s","FIRST" if use_first_state else "BEST (lowest binding energy)")
    logging.info("Interaction curve limit: %d",MAX_CURVES_FOR_INTERACTION_PLOT)
    results: List[FolderResult] = []
    for idx,folder in enumerate(in_folders,start = 1):
        system_id = f"S{idx}"
        r = process_folder(folder,system_id)
        if r is not None:
            results.append(r)
    if not results:
        logging.error("No valid folders processed (missing files or unparsable content).")
        return 2
    #heatmap (PNG + CSV)
    heat_png = outdir / f"interaction_heatmap_{env_mode}.png"
    make_interactions_heatmap(results = results,out_png = heat_png,env_mode = env_mode,use_first_state = use_first_state)
    logging.info("Wrote: %s",heat_png)
    logging.info("Wrote: %s",heat_png.with_suffix(".csv"))
    heat_png = outdir / f"interaction_heatmap_{env_mode}.png"
    heat_csv = heat_png.with_suffix(".csv")
    stats_csv = outdir / f"interaction_heatmap_{env_mode}_stats.csv"
    stats_png = outdir / f"interaction_heatmap_{env_mode}_candlestick.png"
    stats_df = write_heatmap_stats(heat_csv,stats_csv)
    plot_residue_candlestick(stats_df = stats_df,out_png = stats_png,title = f"Residue interaction stats ({env_mode})")
    #stats plots (PNG + CSV per key)
    make_stats_per_key_plots_multi(results,outdir,args.stats_keys)
    logging.info("Wrote stats plots (and CSVs) into: %s",outdir)
    #interaction plot (PNG + CSV)
    inter_plot = outdir / f"interaction_energies_{env_mode}.png"
    make_interactions_plot(results = results,out_png = inter_plot,env_mode = env_mode,use_first_state = use_first_state)
    logging.info("Wrote: %s",inter_plot)
    logging.info("Wrote: %s",inter_plot.with_suffix(".csv"))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
