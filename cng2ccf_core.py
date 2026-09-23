"""
CNG -> CCF SWC back-conversion with automatic scale reconciliation.

Pipeline (unchanged from the original, plus Step 0):
  0. Detect the scale factor between Source and CNG (e.g. nm vs um, voxel vs um)
     and convert the Source to CNG units (microns).
  1. Find 3 anchor points in each tree: root, furthest tip, furthest branch point
     (fallback: midpoint along the root->tip path).
  2. SVD (Kabsch) rotation between the anchor sets.
  3. Apply the inverse rotation to the CNG tree and translate it to the Source root.
Output is always in microns (CNG units).

Contact: NeuroMorpho.Org / Sumit Nanda (snanda@mednet.ucla.edu)
"""
import os
import re
import numpy as np

# --------------------------------------------------------------------------
# I/O
# --------------------------------------------------------------------------
def read_swc(file_path):
    """Return (points array Nx7, header lines)."""
    points, header = [], []
    with open(file_path, "r", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                header.append(s)
                continue
            p = s.split()
            if len(p) < 7:
                continue
            points.append([int(float(p[0])), int(float(p[1])),
                           float(p[2]), float(p[3]), float(p[4]),
                           float(p[5]), int(float(p[6]))])
    if not points:
        raise ValueError(f"No SWC nodes found in {os.path.basename(file_path)}")
    return np.array(points, dtype=float), header


def header_units(header):
    """Look for a declared unit in the header (e.g. navis/FlyWire: 'units': '1 nanometer')."""
    txt = " ".join(header).lower()
    if re.search(r"nanomet|\bnm\b", txt):
        return "nm"
    if re.search(r"micromet|micron|\bum\b|µm", txt):
        return "um"
    return None


def root_index(points):
    idx = np.where(points[:, 6] == -1)[0]
    return int(idx[0]) if len(idx) else 0


# --------------------------------------------------------------------------
# Step 0: scale detection
# --------------------------------------------------------------------------
def _size_metrics(points):
    """Rotation/translation-invariant size measures, robust to small edits."""
    xyz = points[:, 2:5]
    root = xyz[root_index(points)]
    max_reach = np.linalg.norm(xyz - root, axis=1).max()           # furthest point from root
    rg = np.sqrt(((xyz - xyz.mean(axis=0)) ** 2).sum(axis=1).mean())  # radius of gyration
    id2i = {int(n): i for i, n in enumerate(points[:, 0])}
    par = [id2i.get(int(p), -1) for p in points[:, 6]]
    cable = sum(np.linalg.norm(xyz[i] - xyz[j]) for i, j in enumerate(par) if j >= 0)
    return np.array([max_reach, rg, cable])


def detect_scale(source_points, cng_points, header=None,
                 same_tol=0.10, snap_tol=0.03):
    """
    Returns (factor, info) where source_um = source / factor.
    - within +-same_tol of 1       -> 1 (same units)
    - within +-snap_tol of 10^k    -> exactly 10^k (nm/um, mm/um, ...)
    - otherwise                    -> measured value (e.g. voxel -> um, isotropic)
    """
    ms, mc = _size_metrics(source_points), _size_metrics(cng_points)
    valid = mc > 0
    ratios = ms[valid] / mc[valid]
    raw = float(np.median(ratios))
    spread = float(ratios.max() / ratios.min()) if len(ratios) > 1 else 1.0

    if abs(raw - 1) <= same_tol:
        factor, kind = 1.0, "same units"
    else:
        k = round(np.log10(raw))
        if k != 0 and abs(raw / 10 ** k - 1) <= snap_tol:
            factor = float(10 ** k)
            kind = {3: "nanometers -> microns", -3: "millimeters -> microns"}.get(k, f"power of ten (10^{k})")
        else:
            factor, kind = raw, "non-standard (e.g. voxel -> micron), isotropic"

    declared = header_units(header or [])
    info = {"raw_ratio": raw, "factor": factor, "kind": kind,
            "metric_ratios": ratios.tolist(), "metric_spread": spread,
            "declared_units": declared, "warnings": []}
    if declared == "nm" and factor != 1000:
        info["warnings"].append(f"Header says nanometers but measured ratio is {raw:.4g}.")
    if spread > 1.15:
        info["warnings"].append(f"Size metrics disagree (spread {spread:.2f}x): trees may differ "
                                "substantially or scaling may be anisotropic; check the output.")
    return factor, info


# --------------------------------------------------------------------------
# Steps 1-3 (original method)
# --------------------------------------------------------------------------
def find_key_points(points):
    ri = root_index(points)
    root = points[ri, 2:5]
    all_ids, parent_ids = set(points[:, 0]), set(points[:, 6])

    tip_ids = all_ids - parent_ids
    tips = points[np.isin(points[:, 0], list(tip_ids))]
    order = np.argsort(np.linalg.norm(tips[:, 2:5] - root, axis=1))[::-1]
    furthest_tip = tips[order[0], 2:5]

    counts = {}
    for pid in points[:, 6]:
        counts[pid] = counts.get(pid, 0) + 1
    branch_ids = [pid for pid, c in counts.items() if c > 1 and pid != -1]
    branches = points[np.isin(points[:, 0], branch_ids)]

    if branches.shape[0] > 0:
        furthest_branch = branches[np.argmax(np.linalg.norm(branches[:, 2:5] - root, axis=1)), 2:5]
    else:
        id2i = {int(r[0]): i for i, r in enumerate(points)}
        path, cur = [], int(tips[order[0], 0])
        while cur != -1 and cur in id2i:
            i = id2i[cur]
            path.append(points[i, 2:5])
            cur = int(points[i, 6])
        furthest_branch = path[len(path) // 2] if path else root
    return root, furthest_tip, furthest_branch


def compute_svd_rotation(source, target):
    cs, ct = source - source.mean(axis=0), target - target.mean(axis=0)
    U, _, Vt = np.linalg.svd(cs.T @ ct)
    Rm = Vt.T @ U.T
    if np.linalg.det(Rm) < 0:
        Vt[-1, :] *= -1
        Rm = Vt.T @ U.T
    return Rm


def apply_transformation(points, rotation_matrix, translation):
    out = points.copy()
    out[:, 2:5] = points[:, 2:5] @ rotation_matrix.T + translation
    return out


def refine_rotation_icp(cng_xyz, src_xyz, R_init, trim=0.80, max_iter=30,
                        max_angle_deg=15.0, min_gain=0.05, log=print):
    """
    Fine-tune the back-rotation using ALL nodes instead of only the 3 anchors.

    Trimmed ICP: nearest-neighbour correspondence CNG->Source, Kabsch fit on the
    closest `trim` fraction (so CNG soma points and edited branches can't drag the
    fit), iterate. Rotation is about the root, which stays exactly on the Source
    root - no translation is introduced.

    Safety: the refined rotation is returned ONLY if it improves the median
    node-to-node distance by at least `min_gain` and rotates by less than
    `max_angle_deg`. Otherwise the original anchor-based rotation is returned
    unchanged.
    Both inputs must be root-centred. Returns (R, stats).
    """
    from scipy.spatial import cKDTree

    tree = cKDTree(src_xyz)
    med0 = float(np.median(tree.query(cng_xyz @ R_init.T)[0]))
    R, med = R_init.copy(), med0
    for _ in range(max_iter):
        d, idx = tree.query(cng_xyz @ R.T)
        keep = d <= np.quantile(d, trim)
        P, Q = cng_xyz[keep] @ R.T, src_xyz[idx[keep]]
        R_new = compute_svd_rotation(P, Q) @ R
        med_new = float(np.median(tree.query(cng_xyz @ R_new.T)[0]))
        if med_new >= med * (1 - 1e-4):
            break
        R, med = R_new, med_new

    angle = float(np.degrees(np.arccos(np.clip((np.trace(R @ R_init.T) - 1) / 2, -1, 1))))
    stats = {"median_before_um": med0, "median_after_um": med, "angle_deg": angle, "accepted": False}
    if angle > max_angle_deg:
        log(f"  [Note] refinement rejected: {angle:.1f} deg is too large a correction; "
            "keeping the anchor-based rotation.")
        return R_init, stats
    if med > med0 * (1 - min_gain):
        log(f"  [Note] refinement made no useful difference; keeping the anchor-based rotation.")
        return R_init, stats
    stats["accepted"] = True
    log(f"  refined: median node distance {med0:.3f} -> {med:.3f} um "
        f"(extra rotation {angle:.2f} deg)")
    return R, stats


# --------------------------------------------------------------------------
# One pair
# --------------------------------------------------------------------------
def _write_swc(path, points, lines):
    with open(path, "w", newline="\n") as f:
        for l in lines:
            f.write(l + "\n")
        f.write("# For questions and further development, contact NeuroMorpho.Org and "
                "Sumit Nanda at snanda@mednet.ucla.edu\n")
        np.savetxt(f, points, fmt="%d %d %.6f %.6f %.6f %.6f %d")


def convert_pair(source_path, cng_path, output_path, scale_override=None, log=print,
                 original_res_path=None, refine=True):
    src, src_header = read_swc(source_path)
    cng, _ = read_swc(cng_path)

    # Step 0: bring Source into CNG units (microns)
    if scale_override:
        factor, info = float(scale_override), {"kind": "user override", "raw_ratio": float("nan"), "warnings": []}
    else:
        factor, info = detect_scale(src, cng, src_header)
    src[:, 2:5] /= factor
    src[:, 5] /= factor
    log(f"  scale: /{factor:g}  ({info['kind']}; measured ratio {info['raw_ratio']:.5g})")
    for w in info["warnings"]:
        log(f"  [Warning] {w}")

    # Center both trees on their roots
    root_offset = src[root_index(src), 2:5].copy()
    src[:, 2:5] -= root_offset
    cng_root = cng[root_index(cng), 2:5].copy()
    cng_c = cng.copy()
    cng_c[:, 2:5] -= cng_root

    key_s = np.vstack(find_key_points(src))
    key_c = np.vstack(find_key_points(cng_c))
    Rm = compute_svd_rotation(key_s, key_c)
    R_back = np.linalg.inv(Rm)

    # Quality: anchor residual after back-rotation (microns)
    resid = np.linalg.norm(key_c @ R_back.T - key_s, axis=1)
    log(f"  anchor residuals (um): " + ", ".join(f"{r:.2f}" for r in resid))

    # Optional Step 2b: whole-tree refinement of the same rotation (root stays pinned)
    refine_stats = None
    if refine:
        R_back, refine_stats = refine_rotation_icp(cng_c[:, 2:5], src[:, 2:5], R_back, log=log)

    out = apply_transformation(cng_c, R_back, root_offset)

    provenance = [f"# Source: {os.path.basename(source_path)} | CNG: {os.path.basename(cng_path)}",
                  f"# Source scale factor detected: {factor:g} source units per micron ({info['kind']})"]
    if refine_stats and refine_stats["accepted"]:
        provenance.append(f"# Rotation refined on all nodes: median node distance "
                          f"{refine_stats['median_before_um']:.3f} -> {refine_stats['median_after_um']:.3f} um "
                          f"(+{refine_stats['angle_deg']:.2f} deg)")
    _write_swc(output_path,
               out,
               ["# CNG SWC file rotated back to the CCF space (units: micrometers)."] + provenance)

    # Optional second copy in the Source's original resolution (only if scaling was needed)
    original_output = None
    if original_res_path and factor != 1.0:
        scaled = out.copy()
        scaled[:, 2:6] *= factor          # X, Y, Z and radius
        _write_swc(original_res_path,
                   scaled,
                   ["# CNG SWC file rotated back to the CCF space, rescaled to the "
                    "original Source resolution.",
                    f"# Coordinates AND radii multiplied by {factor:g} (1 micron = {factor:g} source units)."]
                   + provenance)
        original_output = original_res_path

    return {"factor": factor, "kind": info["kind"], "warnings": info["warnings"],
            "anchor_residuals_um": resid.tolist(), "output": output_path,
            "original_output": original_output, "refine": refine_stats}


# --------------------------------------------------------------------------
# Pairing & batch
# --------------------------------------------------------------------------
_CNG_SUFFIX = re.compile(r"([._-]?cng)$", re.IGNORECASE)

def pair_key(filename):
    stem = os.path.splitext(os.path.basename(filename))[0]
    return _CNG_SUFFIX.sub("", stem).lower()


def output_stem(cng_filename):
    stem = os.path.splitext(os.path.basename(cng_filename))[0]
    return _CNG_SUFFIX.sub("", stem)


def pair_files(source_files, cng_files):
    """Match by name (X.swc <-> X_CNG.swc / X.CNG.swc); fall back to sorted order."""
    s_map = {pair_key(f): f for f in source_files}
    c_map = {pair_key(f): f for f in cng_files}
    common = sorted(set(s_map) & set(c_map))
    if len(common) == len(source_files) == len(cng_files):
        return [(s_map[k], c_map[k]) for k in common], []
    if len(source_files) == len(cng_files) and not common:
        return list(zip(sorted(source_files), sorted(cng_files))), ["Names did not match; paired by sorted order."]
    unmatched = sorted(set(s_map) ^ set(c_map))
    return [(s_map[k], c_map[k]) for k in common], [f"Unmatched: {', '.join(unmatched)}"] if unmatched else []


def process_swc_files(input_source_folder, input_cng_folder, output_ccf_folder,
                      scale_override=None, log=print, original_res_folder=None, refine=True):
    """original_res_folder: where to also write copies rescaled to the Source's original
    resolution. The folder is created only when at least one pair actually needs scaling.
    Defaults to '<output_ccf_folder>_OriginalResolution'; pass False to disable."""
    os.makedirs(output_ccf_folder, exist_ok=True)
    if original_res_folder is None:
        original_res_folder = output_ccf_folder.rstrip("/\\") + "_OriginalResolution"
    src = [os.path.join(input_source_folder, f) for f in os.listdir(input_source_folder) if f.lower().endswith(".swc")]
    cng = [os.path.join(input_cng_folder, f) for f in os.listdir(input_cng_folder) if f.lower().endswith(".swc")]
    pairs, notes = pair_files(src, cng)
    for n in notes:
        log(f"[Note] {n}")
    results = []
    for s, c in pairs:
        stem = output_stem(c)
        out = os.path.join(output_ccf_folder, f"{stem}.CNG.CCF.swc")
        orig = None
        if original_res_folder:
            os.makedirs(original_res_folder, exist_ok=True)  # removed below if unused
            orig = os.path.join(original_res_folder, f"{stem}.CNG.CCF.OriginalResolution.swc")
        log(f"{os.path.basename(s)}  +  {os.path.basename(c)}")
        try:
            r = convert_pair(s, c, out, scale_override, log, original_res_path=orig, refine=refine)
            results.append(r)
            log(f"  saved: {out}")
            if r["original_output"]:
                log(f"  saved: {r['original_output']}  (original resolution)")
        except Exception as e:
            log(f"  [Error] skipped: {e}")
    if original_res_folder and os.path.isdir(original_res_folder) and not os.listdir(original_res_folder):
        os.rmdir(original_res_folder)  # no pair needed scaling
        log("[Note] No scaling was needed; no OriginalResolution output written.")
    return results


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Rotate CNG SWC files back to CCF/source space (output in microns).")
    ap.add_argument("--source", default="InputSource")
    ap.add_argument("--cng", default="InputCNG")
    ap.add_argument("--out", default="OutputCCF")
    ap.add_argument("--scale", type=float, default=None,
                    help="Force Source/CNG unit ratio (e.g. 1000 for nm->um). Default: auto-detect.")
    ap.add_argument("--no-original-resolution", action="store_true",
                    help="Skip the extra <out>_OriginalResolution copies.")
    ap.add_argument("--no-refine", action="store_true",
                    help="Use the 3-anchor rotation only, without whole-tree refinement.")
    a = ap.parse_args()
    process_swc_files(a.source, a.cng, a.out, a.scale,
                      original_res_folder=False if a.no_original_resolution else None,
                      refine=not a.no_refine)
