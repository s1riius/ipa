# ------------------------------------------------------------
# Utilities to convert Gaussian .cub volumetric density grids into .mrc maps.
# ------------------------------------------------------------

#designed for In-Pocket Analysis:
#- In-pocket runs may produce many "*_dens.cub" files under <output>/geometries/.
#- This module converts all .cub files to .mrc and archives (moves) the original .cub
#  files into a dedicated folder under the analysis output directory.
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional,Tuple,List
import shutil

import numpy as np
import mrcfile
import util.UnitConversion as UC

BOHR_TO_ANGSTROM = UC.Bohr2Angstrom

@dataclass(frozen = True)
class CubeGrid:
    origin_ang: Tuple[float,float,float]
    nsteps: Tuple[int,int,int]
    step_vecs_ang: Tuple[Tuple[float,float,float],
                         Tuple[float,float,float],
                         Tuple[float,float,float]]
    data: np.ndarray                 #shape (nx, ny, nz), float64
    comment: str = ""

def _split_fields(line: str) -> List[str]:
    #cube files can have variable whitespace; splitting on whitespace is safest.
    return [f for f in line.strip().split() if f]

def read_cube(path: Path) -> CubeGrid:
    """Parse a Gaussian cube file (orthogonal or general axes are accepted; however
       this converter currently assumes an orthogonal grid when setting MRC header
       voxel sizes/cell dimensions).
       Data ordering follows cube convention: values are listed with z varying fastest,
       then y, then x; we store as grid[x, y, z]"""
    path = Path(path)
    lines = path.read_text().splitlines()
    if len(lines) < 10:
        raise ValueError(f"Cube file looks too short: {path}")
    comment = "\n".join(lines[:2]).strip()
    #line 3: natoms, origin (in Bohr)
    l3 = _split_fields(lines[2])
    natoms = int(float(l3[0]))
    origin_bohr = (float(l3[1]),float(l3[2]),float(l3[3]))
    origin_ang = tuple(v*BOHR_TO_ANGSTROM for v in origin_bohr)
    #lines 4-6: nsteps and axis vectors (in Bohr)
    axes = []
    nsteps = []
    for i in range(3):
        li = _split_fields(lines[3 + i])
        n = int(float(li[0]))
        vec_bohr = (float(li[1]),float(li[2]),float(li[3]))
        vec_ang = tuple(v*BOHR_TO_ANGSTROM for v in vec_bohr)
        nsteps.append(n)
        axes.append(vec_ang)
    nsteps_t = (nsteps[0],nsteps[1],nsteps[2])
    step_vecs_ang = (axes[0],axes[1],axes[2])
    #skip atom block: natoms lines starting at line 7 (0-based index 6)
    #natoms can be negative in some cube variants; abs() is standard.
    nat = abs(natoms)
    data_start = 6 + nat
    if data_start >= len(lines):
        raise ValueError(f"Cube file atom block extends past EOF: {path}")
    #read all scalar values from remaining lines
    values: List[float] = []
    for line in lines[data_start:]:
        fields = _split_fields(line)
        if not fields:
            continue
        values.extend(float(x) for x in fields)
    nx,ny,nz = nsteps_t
    expected = nx*ny*nz
    if len(values) != expected:
        raise ValueError(f"Cube grid size mismatch for {path}: expected {expected} values "
                         f"({nx}*{ny}*{nz}), got {len(values)}")
    arr = np.asarray(values,dtype = np.float64)
    #cube order: z fastest, then y, then x => reshape (nx, ny, nz) with C-order
    grid = arr.reshape((nx,ny,nz),order = "C")
    return CubeGrid(origin_ang = origin_ang,nsteps = nsteps_t,step_vecs_ang = step_vecs_ang,data = grid,comment = comment)

def _axis_step_and_orthogonality(step_vecs_ang) -> Tuple[Tuple[float, float, float], bool]:
    #return voxel sizes along x,y,z as magnitudes (Angstrom), and whether axes are orthogonal.
    vecs = [np.array(v,dtype = np.float64) for v in step_vecs_ang]
    mags = [float(np.linalg.norm(v)) for v in vecs]
    #orthogonality check: dot products close to 0 compared to magnitudes
    orth = True
    for i in range(3):
        for j in range(i + 1,3):
            dij = float(np.dot(vecs[i],vecs[j]))
            denom = (mags[i]*mags[j]) if mags[i] and mags[j] else 1.0
            if abs(dij) > 1e-6*denom:
                orth = False
    return (mags[0],mags[1],mags[2]),orth

def cube_to_mrc( cube_path: Path, mrc_path: Path, *, scale: Optional[float] = None, dtype: str = "float32", overwrite: bool = True) -> None:
    """Convert a single .cub to .mrc.
       - If scale is None: store float32 (recommended, lossless-ish).
       - If scale is a number: store (scale * density) cast to int16 by default
         (set dtype="int16" explicitly if you want)."""
    cube = read_cube(Path(cube_path))
    voxel_sizes,orthogonal = _axis_step_and_orthogonality(cube.step_vecs_ang)
    if not orthogonal:
        #still write an MRC, but header voxel size/cella will only reflect magnitudes
        #(many viewers assume orthogonal anyway).
        pass
    data = cube.data
    if scale is not None:
        data = (scale*data)
    np_dtype = np.dtype(dtype)
    data_out = data.astype(np_dtype,copy = False)
    data_out = np.transpose(data_out,(2,1,0))
    mrc_path = Path(mrc_path)
    mrc_path.parent.mkdir(parents = True,exist_ok = True)
    with mrcfile.new(str(mrc_path), overwrite=overwrite) as mrc:
        mrc.set_data(data_out)
        #header metadata (best-effort)
        nx,ny,nz = cube.nsteps
        cella = (nx*voxel_sizes[0],ny*voxel_sizes[1],nz*voxel_sizes[2])
        mrc.header.cella = cella
        mrc.voxel_size = voxel_sizes        #(x, y, z) in Angstrom
        #cube origin is the coordinate of the first voxel corner (usually).
        #we store it in the MRC origin field so that maps overlay in viewers.
        mrc.header.origin = cube.origin_ang
        #map axis order:
        #mrcfile defaults to mapc = 1, mapr = 2, maps = 3 (x,y,z); we keep that
        mrc.header.mapc = 1
        mrc.header.mapr = 2
        mrc.header.maps = 3

def convert_cub_tree_to_mrc(cub_dir: Path, mrc_dir: Path, *, glob_pattern: str = "**/*.cub", scale: Optional[float] = None, dtype: str = "float32", overwrite: bool = True, move_cub_to_archive: bool = True) -> int:
    """Convert all .cub files under cub_dir to .mrc.
       IMPORTANT (IPA behavior):
       - .mrc files are written under mrc_dir (NOT next to the .cub),
         preserving the relative directory structure under cub_dir.
       - If move_cub_to_archive is True, cubes are moved into mrc_dir as an archive
         preserving relative structure (legacy mode / standalone usage)."""
    cub_dir = Path(cub_dir)
    mrc_dir = Path(mrc_dir)
    if not cub_dir.exists():
        return 0
    cub_files = [p for p in cub_dir.glob(glob_pattern) if p.is_file() and p.suffix == ".cub"]
    if not cub_files:
        return 0
    converted = 0
    for cub in cub_files:
        rel = cub.relative_to(cub_dir)
        #write MRC into mrc_dir, mirroring the cub relative structure
        mrc_path = (mrc_dir / rel).with_suffix(".mrc")
        mrc_path.parent.mkdir(parents = True,exist_ok = True)
        cube_to_mrc(cub,mrc_path,scale = scale,dtype = dtype,overwrite = overwrite)
        converted += 1
        if move_cub_to_archive:
            #archive the cub into mrc_dir as well (legacy behavior)
            dest = (mrc_dir / rel)
            dest.parent.mkdir(parents = True,exist_ok = True)
            if dest.exists():
                dest.unlink()
            shutil.move(str(cub),str(dest))
    return converted

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description = "Convert .cub files to .mrc (IPA utility).")
    ap.add_argument("--cub-folder",required = True,help = "Path to folder with cub files")
    ap.add_argument("--mrc-folder",required = True,help = "Path where mrc files will be saved")
    ap.add_argument("--scale",type = float,default = None,help = "Optional scaling factor (e.g., 600) before writing.")
    ap.add_argument("--dtype",default = "float32",help = "Output dtype (float32 recommended; int16 if using scale).")
    ap.add_argument("--no-move",action = "store_true",help = "Do not move .cub files to archive.")
    args = ap.parse_args()
    n = convert_cub_tree_to_mrc(Path(args.cub_folder),Path(args.mrc_folder),scale = args.scale,dtype = args.dtype,move_cub_to_archive = (not args.no_move))
    print(f"Converted {n} cube file(s).")
