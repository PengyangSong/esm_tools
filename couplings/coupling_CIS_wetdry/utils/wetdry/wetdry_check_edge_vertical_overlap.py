#!/usr/bin/env python3
"""Validate vertical-layer overlap across every edge of a wet/dry mesh.

This is deliberately a validation step.  It does not alter the mesh.  FESOM
can receive very large transports when neighboring columns have no common
active layer, so the caller should repair the bathymetry or reject the mesh
before running ``fesom_ini.x``.
"""

from __future__ import annotations

import argparse
import os

import numpy as np


def read_column(path: str) -> np.ndarray:
    values = np.loadtxt(path, dtype=float)
    return values[:, 0] if values.ndim > 1 else values


def read_mesh(submesh: str):
    with open(os.path.join(submesh, "nod2d.out")) as f:
        nnode = int(f.readline().strip())
    with open(os.path.join(submesh, "elem2d.out")) as f:
        nelem = int(f.readline().strip())
        elem = np.loadtxt(f, dtype=int, max_rows=nelem)
    elem = np.atleast_2d(elem) - 1
    if elem.shape != (nelem, 3):
        raise ValueError(f"elem2d.out has shape {elem.shape}, expected ({nelem}, 3)")
    return nnode, elem


def read_zmid(path: str) -> np.ndarray:
    with open(path) as f:
        nz = int(f.readline().strip())
        zbar = np.array([float(f.readline()) for _ in range(nz)])
    if zbar.size < 2:
        raise ValueError(f"need at least two levels in {path}")
    return 0.5 * (zbar[:-1] + zbar[1:])


def first_layer(depth: float, zmid: np.ndarray) -> int:
    found = np.where(depth > zmid)[0]
    return int(found[0]) if found.size else len(zmid) - 1


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--submesh_dir", required=True)
    parser.add_argument("--depth_zlev", default=None)
    parser.add_argument("--max_report", type=int, default=20)
    parser.add_argument("--depth_at_node_positive_down", action="store_true")
    parser.add_argument("--zmid_convention", choices=("auto", "negative_z", "positive_down"), default="auto")
    args = parser.parse_args()

    sub = args.submesh_dir.rstrip(os.sep)
    nnode, elem = read_mesh(sub)
    depth = read_column(os.path.join(sub, "depth@node.out"))
    cavity = read_column(os.path.join(sub, "cavity_depth@node.out"))
    zpath = args.depth_zlev
    if zpath is None:
        for candidate in (os.path.join(sub, "global", "depth_zlev.out"), os.path.join(sub, "depth_zlev.out")):
            if os.path.isfile(candidate):
                zpath = candidate
                break
    if zpath is None:
        raise FileNotFoundError("depth_zlev.out not found")
    zmid = read_zmid(zpath)

    if len(depth) != nnode or len(cavity) != nnode:
        raise ValueError("depth@node/cavity_depth length does not match nod2d.out")
    if args.zmid_convention == "auto":
        negative_z = bool(zmid[-1] < zmid[0])
    else:
        negative_z = args.zmid_convention == "negative_z"

    # Convert both boundaries to the vertical-coordinate convention used by zmid.
    if negative_z:
        base_coord = -cavity
        bed_coord = depth
    else:
        base_coord = cavity
        bed_coord = depth if args.depth_at_node_positive_down else np.where(depth < 0, -depth, depth)

    first = np.array([first_layer(float(x), zmid) for x in base_coord], dtype=int)
    last = np.array([first_layer(float(x), zmid) for x in bed_coord], dtype=int)
    if np.any(last < first):
        bad_nodes = np.where(last < first)[0]
        raise ValueError(f"{len(bad_nodes)} nodes have bed above cavity base; first={first[bad_nodes[:5]]}, last={last[bad_nodes[:5]]}")

    edges = np.vstack((elem[:, [0, 1]], elem[:, [1, 2]], elem[:, [2, 0]]))
    edges = np.unique(np.sort(edges, axis=1), axis=0)
    overlap = np.minimum(last[edges[:, 0]], last[edges[:, 1]]) - np.maximum(first[edges[:, 0]], first[edges[:, 1]]) + 1
    bad = np.where(overlap < 1)[0]

    print(f"wetdry_check_edge_vertical_overlap: nodes={nnode} edges={len(edges)} bad_edges={len(bad)}")
    for idx in bad[: max(0, args.max_report)]:
        a, b = edges[idx]
        print(f"  edge {a + 1} {b + 1}: node_a=[{first[a]},{last[a]}] node_b=[{first[b]},{last[b]}]")
    if len(bad):
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
