# ------------------------------------------------------------
# Functions to determine molecular connectivity and topology
# ------------------------------------------------------------

from typing import List,Optional
import util.AtomicData as AD

def _nearest_index_by_element_and_coord(elements: List[str], coords: List, target_elem: str, target_coord, max_dist: float = 1.2) -> Optional[int]:
    """Find the closest atom index (matching element) to target_coord within max_dist (Å)."""
    te = str(target_elem).strip().upper()
    tx,ty,tz = (float(target_coord[0]),float(target_coord[1]),float(target_coord[2]))
    max_d2 = float(max_dist)**2
    best_i = None
    best_d2 = float("inf")
    for i,el in enumerate(elements or []):
        if str(el).strip().upper() != te:
            continue
        try:
            x,y,z = coords[i]
            dx = float(x) - tx
            dy = float(y) - ty
            dz = float(z) - tz
            d2 = dx*dx + dy*dy + dz*dz
            if d2 < best_d2:
                best_d2 = d2
                best_i = i
        except (TypeError,ValueError,IndexError): continue
    if best_i is None or best_d2 > max_d2:
        return None
    return int(best_i)

def _count_covalent_neighbors(elements: List[str], coords: List, idx: int, scale: float = 1.3) -> Optional[int]:
    """Count neighbours using a covalent-radius distance criterion (no bonds required).
       We use squared distances to avoid sqrt, and compare against:
           d^2 < ( scale * (r_i + r_j) )^2
       where r are covalent radii in Å."""
    if idx is None:
        return None
    idx = int(idx)
    if idx < 0 or idx >= len(elements):
        return None
    ri = AD._covalent_radius_angstrom(elements[idx])
    if ri is None:
        return None
    xi,yi,zi = coords[idx]
    xi = xi
    yi = yi
    zi = zi
    s = scale
    n = 0
    for j,elj in enumerate(elements):
        if j == idx:
            continue
        rj = AD._covalent_radius_angstrom(elj)
        if rj is None:
            continue
        try:
            xj,yj,zj = coords[j]
            dx = xj - xi
            dy = yj - yi
            dz = zj - zi
            d2 = dx*dx + dy*dy + dz*dz
            thr = s*(ri + rj)
            if d2 < (thr*thr):
                n += 1
        except (TypeError,ValueError,IndexError):
            continue
    return int(n)

__all__ = [
    "_nearest_index_by_element_and_coord",
    "_count_covalent_neighbors"
]