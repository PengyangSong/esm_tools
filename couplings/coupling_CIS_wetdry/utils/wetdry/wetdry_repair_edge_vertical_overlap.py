#!/usr/bin/env python3
"""Repair missing active-layer overlap between neighboring wet/dry columns.

This follows the FESOM--ISSM strategy: for every bad edge, deepen the
shallower column's bed to the bottom of the first layer reached by the deeper
column.  The shelf-base field is not changed.  The repair is bounded and is
followed by a complete re-check; unresolved edges are fatal.
"""

from __future__ import annotations

import argparse
import os

import numpy as np


def read_column(path):
    a = np.loadtxt(path, dtype=float)
    return a[:, 0] if a.ndim > 1 else a


def read_mesh(sub):
    with open(os.path.join(sub, "nod2d.out")) as f:
        n = int(f.readline().strip())
    with open(os.path.join(sub, "elem2d.out")) as f:
        ne = int(f.readline().strip())
        elem = np.loadtxt(f, dtype=int, max_rows=ne).reshape((-1, 3)) - 1
    edges = np.vstack((elem[:, [0, 1]], elem[:, [1, 2]], elem[:, [2, 0]]))
    return n, np.unique(np.sort(edges, axis=1), axis=0)


def read_zmid(path):
    with open(path) as f:
        nz = int(f.readline().strip())
        zbar = np.array([float(f.readline()) for _ in range(nz)])
    return 0.5 * (zbar[:-1] + zbar[1:])


def layer_index(value, zmid):
    found = np.where(value > zmid)[0]
    return int(found[0]) if found.size else len(zmid) - 1


def active_layers(depth, cavity, zmid, positive_down_file=False, zconv="auto"):
    if zconv == "auto":
        negative_z = bool(zmid[-1] < zmid[0])
    else:
        negative_z = zconv == "negative_z"
    if negative_z:
        base = -cavity
        bed = depth
    else:
        base = cavity
        bed = depth if positive_down_file else np.where(depth < 0.0, -depth, depth)
    first = np.array([layer_index(x, zmid) for x in base], dtype=int)
    last = np.array([layer_index(x, zmid) for x in bed], dtype=int)
    return first, last, negative_z


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--submesh_dir", required=True)
    p.add_argument("--depth_zlev", default=None)
    p.add_argument("--depth_at_node_positive_down", action="store_true")
    p.add_argument("--zmid_convention", choices=("auto", "negative_z", "positive_down"), default="auto")
    p.add_argument("--depth_adjust_eps", type=float, default=0.1)
    p.add_argument("--max_bed_deepening", type=float, default=500.0,
                   help="maximum permitted deepening per node in metres")
    p.add_argument("--max_passes", type=int, default=10)
    args = p.parse_args()

    sub = args.submesh_dir.rstrip(os.sep)
    depth_path = os.path.join(sub, "depth@node.out")
    cavity_path = os.path.join(sub, "cavity_depth@node.out")
    zpath = args.depth_zlev
    if zpath is None:
        for candidate in (os.path.join(sub, "global", "depth_zlev.out"), os.path.join(sub, "depth_zlev.out")):
            if os.path.isfile(candidate):
                zpath = candidate
                break
    if not zpath:
        raise FileNotFoundError("depth_zlev.out not found")

    nnode, edges = read_mesh(sub)
    depth = read_column(depth_path).copy()
    cavity = read_column(cavity_path)
    zmid = read_zmid(zpath)
    if len(depth) != nnode or len(cavity) != nnode:
        raise ValueError("depth@node/cavity_depth length mismatch")

    before = depth.copy()
    first, last, negative_z = active_layers(
        depth, cavity, zmid, args.depth_at_node_positive_down, args.zmid_convention
    )
    repairs = []
    for _pass in range(args.max_passes):
        overlap = np.minimum(last[edges[:, 0]], last[edges[:, 1]]) - np.maximum(
            first[edges[:, 0]], first[edges[:, 1]]
        ) + 1
        bad = np.where(overlap < 1)[0]
        if not len(bad):
            break
        changed = False
        for bi in bad:
            a, b = edges[bi]
            # Select the shallower bed; in negative-z coordinates this is the
            # larger value, while in positive-down coordinates it is smaller.
            if negative_z:
                shallow, deep = (a, b) if depth[a] >= depth[b] else (b, a)
                target_level = min(first[deep] + 1, len(zmid) - 1)
                new_depth = float(zmid[target_level] - args.depth_adjust_eps)
                delta = max(0.0, depth[shallow] - new_depth)
                if delta <= 0.0:
                    continue
                if delta > args.max_bed_deepening:
                    raise RuntimeError(
                        f"edge {shallow + 1}-{deep + 1} requires {delta:.3f} m deepening "
                        f"> max_bed_deepening={args.max_bed_deepening}"
                    )
                depth[shallow] = new_depth
            else:
                depth_pd = depth if args.depth_at_node_positive_down else np.where(depth < 0.0, -depth, depth)
                shallow, deep = (a, b) if depth_pd[a] <= depth_pd[b] else (b, a)
                target_level = min(first[deep] + 1, len(zmid) - 1)
                new_pd = float(zmid[target_level] + args.depth_adjust_eps)
                delta = max(0.0, new_pd - depth_pd[shallow])
                if delta <= 0.0:
                    continue
                if delta > args.max_bed_deepening:
                    raise RuntimeError(
                        f"edge {shallow + 1}-{deep + 1} requires {delta:.3f} m deepening "
                        f"> max_bed_deepening={args.max_bed_deepening}"
                    )
                depth[shallow] = new_pd if args.depth_at_node_positive_down else -new_pd
            repairs.append((shallow + 1, deep + 1, delta))
            changed = True
        if not changed:
            break
        first, last, _ = active_layers(
            depth, cavity, zmid, args.depth_at_node_positive_down, args.zmid_convention
        )
    else:
        raise RuntimeError(f"edge-overlap repair exceeded max_passes={args.max_passes}")

    overlap = np.minimum(last[edges[:, 0]], last[edges[:, 1]]) - np.maximum(
        first[edges[:, 0]], first[edges[:, 1]]
    ) + 1
    bad = np.where(overlap < 1)[0]
    if len(bad):
        raise RuntimeError(f"edge-overlap repair left {len(bad)} bad edges")

    np.savetxt(depth_path + ".before_edge_repair", before, fmt="%f")
    np.savetxt(depth_path, depth, fmt="%f")
    # fesom_ini.x consumes aux3d.out, so keep it synchronized with the
    # repaired depth@node field rather than leaving the pre-repair bathymetry.
    with open(zpath) as fz, open(os.path.join(sub, "aux3d.out"), "w") as fa:
        fa.write(fz.read())
        with open(depth_path) as fd:
            fa.write(fd.read())
    with open(os.path.join(sub, "edge_overlap_repair_diagnostics.txt"), "w") as f:
        f.write(f"repaired_edges={len(repairs)}\n")
        f.write(f"changed_nodes={int(np.count_nonzero(depth != before))}\n")
        f.write(f"max_deepening_m={float(np.max(np.abs(depth - before))):.8e}\n")
        f.write(f"remaining_bad_edges={len(bad)}\n")
        for a, b, delta in repairs:
            f.write(f"edge_shallow_node={a} edge_deep_node={b} deepening_m={delta:.8e}\n")
    print(
        f"wetdry_repair_edge_vertical_overlap: repaired_edges={len(repairs)} "
        f"changed_nodes={np.count_nonzero(depth != before)} "
        f"max_deepening_m={np.max(np.abs(depth - before)):.3f} remaining_bad_edges=0"
    )


if __name__ == "__main__":
    main()
