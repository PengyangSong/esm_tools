#!/usr/bin/env python3
"""
After depth@node.out and cavity_depth@node.out exist: ensure at least N vertical layer
midpoints lie between ice-shelf base and ocean bed by deepening bed (depth@node) where needed.

Uses depth_zlev.out (zbar interfaces) like aux3d / fesom2_read_mesh:
  zmid[k] = 0.5 * (zbar[k] + zbar[k+1]), length len(zbar)-1.

  Typical FESOM files use negative z (0 at surface, more negative downward), e.g. zmid[0]≈-2.5,
  zmid[-1]≈-6000. If zmid[-1] < zmid[0], the script uses the same z convention for bed as
  depth@node (modify_nodhn: negative ocean depth) and maps cavity_depth (positive metres
  below surface) to z as -nbase before comparing to zmid.

  If zmid increases with depth (positive-down midpoints), comparisons use positive-down
  bathymetry: nbed_pd = -nbed when depth@node is negative-ocean, unless
  --depth_at_node_positive_down.

Per-node logic (from user recipe, with guards):
  - nbase = cavity_depth@node (shelf base, positive down from surface in m)
  - nbed  = depth@node (negative z / ocean depth unless positive-down legacy)
  - By default, any node with cavity_depth > cavity_floor is considered.
  - --thin_cavity_only: nbase shallower than the magnitude of the uppermost layer mid
    (nbase < -zmid[0] if zmid negative, else nbase < zmid[0]).
  - base_level = first j with nbase_cmp > zmid[j]; bed_level = first j with nbed_cmp > zmid[j]
    (nbase_cmp = -nbase in negative-z grids; nbed_cmp = nbed when z matches depth@node).
  - If bed_level - base_level < min_layer_span: set bed to zmid[base_level + bed_index_offset] - eps

Reads depth_zlev.out from submesh/global/depth_zlev.out if present, else submesh/depth_zlev.out.
Rewrites depth@node.out and rebuilds aux3d.out (same zlev content + depth@node).

Requires: numpy
"""

from __future__ import annotations

import argparse
import os
from typing import Optional

import numpy as np


def read_nod2d_n(path):
    with open(path) as f:
        n = int(f.readline().strip())
    return n


def read_column(path):
    a = np.loadtxt(path, dtype=np.float64)
    if a.ndim > 1:
        a = a[:, 0]
    return a


def read_depth_zlev(path):
    """First line: nz (number of zbar values); next nz lines: zbar."""
    with open(path) as f:
        nz = int(f.readline().strip())
        zbar = np.array([float(f.readline()) for _ in range(nz)], dtype=np.float64)
    if zbar.size < 2:
        raise ValueError(f"need at least 2 zbar levels in {path}")
    zmid = 0.5 * (zbar[:-1] + zbar[1:])
    return zbar, zmid


def zmid_is_fesom_negative_z(zmid: np.ndarray) -> bool:
    """True when z becomes more negative with index (standard FESOM z, 0 at surface)."""
    return bool(zmid[-1] < zmid[0])


def first_index_where_depth_gt_zmid(depth: float, zmid: np.ndarray) -> Optional[int]:
    """First j such that depth > zmid[j] (user recipe)."""
    w = np.where(depth > zmid)[0]
    if w.size == 0:
        return None
    return int(w[0])


def to_positive_down_bed(nd: float, positive_down_file: bool) -> float:
    """Map depth@node to positive-down bathymetry depth for zmid layer logic."""
    if positive_down_file:
        return float(nd)
    nd = float(nd)
    return -nd if nd < 0.0 else nd


def from_positive_down_bed(nd_pd_new: float, nd_file_was_negative: bool) -> float:
    """Write back depth@node after deepening in positive-down space."""
    if nd_file_was_negative:
        return -float(nd_pd_new)
    return float(nd_pd_new)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--submesh_dir", required=True)
    p.add_argument("--min_layer_span", type=int, default=3,
                   help="require bed_level - base_level >= this (default 3)")
    p.add_argument("--bed_index_offset", type=int, default=2,
                   help="new bed level index = base_level + this when spanning too small (user: 2)")
    p.add_argument("--depth_adjust_eps", type=float, default=0.1,
                   help="new nbed = zmid[bed_level_new] - this (m)")
    p.add_argument(
        "--thin_cavity_only",
        action="store_true",
        help="only nodes with cavity_depth < zmid[0] (shallow shelf base vs first layer mid)",
    )
    p.add_argument("--cavity_floor", type=float, default=0.0,
                   help="treat cavity_depth <= this as no cavity")
    p.add_argument(
        "--depth_zlev",
        default=None,
        help="path to depth_zlev.out (default: <submesh>/global/depth_zlev.out or <submesh>/depth_zlev.out)",
    )
    p.add_argument(
        "--depth_at_node_positive_down",
        action="store_true",
        help="depth@node is positive down (skip flip); default expects negative ocean depth like modify_nodhn",
    )
    p.add_argument(
        "--zmid_convention",
        choices=("auto", "negative_z", "positive_down"),
        default="auto",
        help="vertical coordinate of zmid: FESOM negative z (default auto: zmid[-1]<zmid[0]) or positive-down m",
    )
    args = p.parse_args()

    sub = args.submesh_dir.rstrip("/") + "/"
    nod2d = os.path.join(sub, "nod2d.out")
    d_node = os.path.join(sub, "depth@node.out")
    d_cav = os.path.join(sub, "cavity_depth@node.out")

    for path in (nod2d, d_node, d_cav):
        if not os.path.isfile(path):
            raise FileNotFoundError(path)

    if args.depth_zlev:
        dz_path = args.depth_zlev
    else:
        dz_global = os.path.join(sub, "global", "depth_zlev.out")
        dz_root = os.path.join(sub, "depth_zlev.out")
        if os.path.isfile(dz_global):
            dz_path = dz_global
        elif os.path.isfile(dz_root):
            dz_path = dz_root
        else:
            raise FileNotFoundError(
                f"depth_zlev.out not found: try {dz_global} or {dz_root}"
            )

    mesh_dim_n = read_nod2d_n(nod2d)
    zbar, zmid = read_depth_zlev(dz_path)
    with open(dz_path) as _f:
        depth_zlev_text = _f.read()
    nbed = read_column(d_node).copy()
    nbase = read_column(d_cav)

    if nbed.shape[0] != mesh_dim_n or nbase.shape[0] != mesh_dim_n:
        raise ValueError("depth@node / cavity_depth length mismatch vs nod2d")

    if args.zmid_convention == "auto":
        z_neg = zmid_is_fesom_negative_z(zmid)
    elif args.zmid_convention == "negative_z":
        z_neg = True
    else:
        z_neg = False

    z_top_for_thin = float(-zmid[0]) if z_neg else float(zmid[0])

    n_fix = 0
    pos_down = args.depth_at_node_positive_down
    for i in range(mesh_dim_n):
        nb = float(nbase[i])
        nd = float(nbed[i])

        if nb <= args.cavity_floor:
            continue

        if args.thin_cavity_only:
            if not (nb < z_top_for_thin):
                continue

        if z_neg:
            nbase_cmp = -nb
            nbed_cmp = nd
            base_level = first_index_where_depth_gt_zmid(nbase_cmp, zmid)
            bed_level = first_index_where_depth_gt_zmid(nbed_cmp, zmid)
        else:
            nd_pd = to_positive_down_bed(nd, pos_down)
            base_level = first_index_where_depth_gt_zmid(nb, zmid)
            bed_level = first_index_where_depth_gt_zmid(nd_pd, zmid)

        if base_level is None:
            base_level = 0
        if bed_level is None:
            bed_level = len(zmid) - 1

        if bed_level - base_level >= args.min_layer_span:
            continue

        bed_level_new = base_level + args.bed_index_offset
        if bed_level_new >= len(zmid):
            bed_level_new = len(zmid) - 1
        if bed_level_new <= base_level:
            bed_level_new = min(base_level + args.min_layer_span, len(zmid) - 1)

        nbed_new = float(zmid[bed_level_new] - args.depth_adjust_eps)
        if z_neg:
            if nbed_new < nd:
                nbed[i] = nbed_new
                n_fix += 1
        elif nbed_new > nd_pd:
            file_was_neg = (not pos_down) and nd < 0.0
            nbed[i] = from_positive_down_bed(nbed_new, file_was_neg)
            n_fix += 1

    np.savetxt(d_node, nbed, fmt="%f")
    with open(d_node) as f_n:
        depth_node_text = f_n.read()
    aux3d_path = os.path.join(sub, "aux3d.out")
    with open(aux3d_path, "w") as out:
        out.write(depth_zlev_text)
        out.write(depth_node_text)

    _zconv = "negative_z" if z_neg else "positive_down"
    print(
        f"wetdry_submesh_cavity_min_layers: nodes={mesh_dim_n} zmid_levels={len(zmid)} "
        f"zmid_conv={_zconv} deepened_bed={n_fix} (min_layer_span={args.min_layer_span})"
    )


if __name__ == "__main__":
    main()
