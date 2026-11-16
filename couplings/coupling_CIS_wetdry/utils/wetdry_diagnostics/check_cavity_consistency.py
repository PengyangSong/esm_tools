#!/usr/bin/env python3
"""Check node/element cavity-layer consistency around selected nodes."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


def counted(path: Path, dtype=int):
    with path.open() as handle:
        count = int(handle.readline())
    return np.loadtxt(path, skiprows=1, max_rows=count, dtype=dtype)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("blowup", type=Path)
    parser.add_argument("mesh", type=Path)
    parser.add_argument("--node", type=int, nargs="+", required=True)
    args = parser.parse_args()

    nodes = counted(args.mesh / "nod2d.out", dtype=float)
    elements = counted(args.mesh / "elem2d.out", dtype=int) - 1
    elvl = np.loadtxt(args.mesh / "elvls.out", dtype=int) - 1
    ecav = np.loadtxt(args.mesh / "cavity_elvls.out", dtype=int) - 1
    nlvl = np.loadtxt(args.mesh / "nlvls.out", dtype=int) - 1
    ncav = np.loadtxt(args.mesh / "cavity_nlvls.out", dtype=int) - 1

    with Dataset(args.blowup) as nc:
        helem = np.asarray(nc.variables["helem"][0])
        u = np.asarray(nc.variables["u"][0])
        v = np.asarray(nc.variables["v"][0])
        rhs = np.asarray(nc.variables["ssh_rhs"][0])

    target = {item - 1 for item in args.node}
    selected = [i for i, tri in enumerate(elements) if any(n in tri for n in target)]
    print(f"Selected elements: {len(selected)}")
    print("element nodes cavity_el raw node-cavity(raw) active-thickness zero-active max|uv| max|uv_rhs|")

    for e in selected:
        tri = elements[e]
        node_cav = ncav[tri]
        # FESOM output uses one-based layer counts; convert back for display.
        ecav_raw, elvl_raw = int(ecav[e] + 1), int(elvl[e] + 1)
        node_cav_raw = ",".join(str(int(item + 1)) for item in node_cav)
        start, stop = max(int(ecav[e]), 0), max(int(elvl[e]), 0)
        thickness = helem[e, start:stop]
        vel = np.hypot(u[e, start:stop], v[e, start:stop])
        # The blow-up file does not contain u_rhs/v_rhs under one common name;
        # use ssh_rhs as the node-level indicator and report element velocity.
        zero_active = int(np.sum(thickness <= 0)) if thickness.size else 0
        mismatch = ecav[e] < node_cav.min() or ecav[e] > node_cav.max()
        print(
            f"{e + 1:7d} {','.join(str(int(n + 1)) for n in tri):>18s} "
            f"{ecav_raw:4d}/{elvl_raw:<4d} {node_cav_raw:>18s} "
            f"{thickness.min() if thickness.size else np.nan: .4e}/"
            f"{thickness.max() if thickness.size else np.nan: .4e} "
            f"{zero_active:10d} {vel.max() if vel.size else np.nan: .4e} "
            f"mismatch={mismatch}"
        )

    print("\nSelected node ssh_rhs:")
    for node in sorted(target):
        print(f"  node {node + 1}: lon={nodes[node,1]:.8f} "
              f"lat={nodes[node,2]:.8f} ssh_rhs={rhs[node]:.6e} "
              f"cavity/nlvls(raw)={ncav[node] + 1}/{nlvl[node] + 1}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
