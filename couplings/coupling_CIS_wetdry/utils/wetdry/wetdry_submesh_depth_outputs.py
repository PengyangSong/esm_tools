#!/usr/bin/env python3
"""
After submesh_mapping: map maxmesh bedrock and PISM shelf to submesh nodes and write
depth@node.out, cavity_depth@node.out, then build aux3d.out (depth_zlev + depth@node).

  depth@node:          same as modify_nodhn.py: -nodhn from submesh nodhn.out when that
                       file exists and matches the submesh node count (negative = ocean
                       depth, more negative = deeper). Otherwise -bedrock_sub mapped from
                       QGLB (legacy fallback).
  cavity_depth@node:   shelf_sub     (ice-shelf base level on maxmesh, mapped)

Applies the same minimum-depth clamp as modify_nodhn.py (default -20 m).

Requires: numpy
"""

from __future__ import annotations

import argparse
import os
import shutil

import numpy as np


def read_map_nod(path):
    m = np.loadtxt(path, dtype=np.int64)
    if m.ndim > 1:
        m = m[:, 0]
    return m


def read_column(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(path)
    a = np.loadtxt(path, dtype=np.float64)
    if a.ndim > 1:
        a = a[:, 0]
    return a


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--submesh_dir", required=True)
    p.add_argument("--qglb_global", required=True, help="QGLB global (bedrock_topography.out)")
    p.add_argument("--maxmesh_shelf_file", required=True, help="wetdry_shelf_on_maxmesh.out")
    p.add_argument("--bedrock_file", default="bedrock_topography.out")
    p.add_argument("--map_nod", default="map_nod.out")
    p.add_argument("--min_depth", type=float, default=-20.0,
                   help="clamp depth@node (same convention as modify_nodhn)")
    p.add_argument(
        "--nodhn_file",
        default="nodhn.out",
        help="basename under submesh_dir; if present and length matches submesh, depth@node=-nodhn",
    )
    p.add_argument(
        "--no_nodhn",
        action="store_true",
        help="force depth@node from mapped bedrock only (ignore nodhn.out)",
    )
    args = p.parse_args()

    sub = args.submesh_dir.rstrip("/") + "/"
    gdir = args.qglb_global.rstrip("/")

    map_path = os.path.join(sub, args.map_nod)
    map_nod = read_map_nod(map_path)
    n_sub = map_nod.shape[0]

    bedrock_max = read_column(os.path.join(gdir, args.bedrock_file))
    shelf_max = read_column(args.maxmesh_shelf_file)
    if shelf_max.shape[0] != bedrock_max.shape[0]:
        raise ValueError("shelf maxmesh file length mismatch vs bedrock")

    # map_nod is 1-based maxmesh indices
    j = map_nod.astype(np.int64) - 1
    if np.any(j < 0) or np.any(j >= bedrock_max.shape[0]):
        raise ValueError("map_nod indices out of range for maxmesh fields")

    bed_sub = bedrock_max[j]
    shelf_sub = shelf_max[j]
    if bed_sub.shape[0] != n_sub or shelf_sub.shape[0] != n_sub:
        raise ValueError("mapped array size mismatch")

    nodhn_path = os.path.join(sub, args.nodhn_file)
    use_nodhn = (not args.no_nodhn) and os.path.isfile(nodhn_path)
    if use_nodhn:
        nodhn_sub = read_column(nodhn_path)
        if nodhn_sub.shape[0] != n_sub:
            print(
                f"wetdry_submesh_depth_outputs: {nodhn_path} length {nodhn_sub.shape[0]} "
                f"!= submesh {n_sub} — using mapped bedrock for depth@node",
                flush=True,
            )
            use_nodhn = False

    if use_nodhn:
        depth_node = (-1.0) * nodhn_sub
    else:
        depth_node = (-1.0) * bed_sub
    depth_node = np.where(depth_node > args.min_depth, args.min_depth, depth_node)

    depth_cav = shelf_sub

    cwd = os.getcwd()
    try:
        os.chdir(sub)
        np.savetxt("depth@node.out", depth_node, fmt="%f")
        np.savetxt("cavity_depth@node.out", depth_cav, fmt="%f")

        dz = os.path.join("global", "depth_zlev.out")
        if not os.path.isfile(dz):
            raise FileNotFoundError(f"missing {os.path.join(sub, dz)}")
        shutil.copyfile(dz, "depth_zlev.out")
        # aux3d: vertical structure + surface depth (same as cutting pipeline)
        with open("depth_zlev.out") as f_z, open("depth@node.out") as f_n, open(
            "aux3d.out", "w"
        ) as out:
            out.write(f_z.read())
            out.write(f_n.read())
    finally:
        os.chdir(cwd)

    _src = args.nodhn_file if use_nodhn else "mapped_bedrock"
    print(
        f"wetdry_submesh_depth_outputs: submesh nodes={n_sub} depth@node from {_src} "
        f"(negative ocean, clamp min={args.min_depth}) — wrote depth@node.out, "
        f"cavity_depth@node.out, aux3d.out under {sub}"
    )


if __name__ == "__main__":
    main()
