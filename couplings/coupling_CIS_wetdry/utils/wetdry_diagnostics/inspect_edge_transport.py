#!/usr/bin/env python3
"""Inspect edge transport near an FESOM blow-up cluster.

The transport estimate uses a local lon/lat tangent-plane geometry. It is
diagnostic, not a replacement for FESOM's exact spherical edge metrics.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


def counted_table(path: Path, dtype=int):
    with path.open() as handle:
        count = int(handle.readline())
    return np.loadtxt(path, skiprows=1, max_rows=count, dtype=dtype)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("blowup", type=Path)
    parser.add_argument("mesh", type=Path)
    parser.add_argument("--node", type=int, nargs="+", required=True,
                        help="1-based node IDs to inspect")
    parser.add_argument("--limit", type=int, default=30)
    args = parser.parse_args()

    nodes = counted_table(args.mesh / "nod2d.out", dtype=float)
    elements = counted_table(args.mesh / "elem2d.out", dtype=int) - 1
    with (args.mesh / "edgenum.out").open() as handle:
        edge_count = int(handle.readline())
    edges = np.loadtxt(args.mesh / "edges.out", max_rows=edge_count, dtype=int) - 1
    edge_tri = np.loadtxt(args.mesh / "edge_tri.out", max_rows=edge_count, dtype=int) - 1
    elvl = np.loadtxt(args.mesh / "elvls.out", dtype=int) - 1
    ecav = np.loadtxt(args.mesh / "cavity_elvls.out", dtype=int) - 1

    # Local tangent-plane coordinates in metres. This preserves the sign of
    # the edge-normal estimate while avoiding a dependency on cartopy/pyfesom.
    lon = np.deg2rad(nodes[:, 1])
    lat = np.deg2rad(nodes[:, 2])
    radius = 6.371e6
    x = radius * lon * np.cos(np.mean(lat))
    y = radius * lat
    centroids = np.column_stack((x[elements].mean(axis=1), y[elements].mean(axis=1)))

    target = {item - 1 for item in args.node}
    selected = []
    for edge_id, pair in enumerate(edges):
        if any(node in pair for node in target):
            selected.append(edge_id)

    with Dataset(args.blowup) as nc:
        u = np.asarray(nc.variables["u"][0])
        v = np.asarray(nc.variables["v"][0])
        helem = np.asarray(nc.variables["helem"][0])
        ssh_rhs = np.asarray(nc.variables["ssh_rhs"][0])

    records = []
    for edge_id in selected:
        n1, n2 = edges[edge_id]
        e1, e2 = edge_tri[edge_id]
        midpoint = 0.5 * (np.array([x[n1], y[n1]]) + np.array([x[n2], y[n2]]))
        fluxes = []
        elem_ids = [e1] if e2 < 0 else [e1, e2]
        for elem_id in elem_ids:
            delta = midpoint - centroids[elem_id]
            # FESOM's term is (v*delta_x - u*delta_y)*helem.
            levels = max(int(elvl[elem_id]), 0)
            flux = np.sum((v[elem_id, :levels] * delta[0]
                           - u[elem_id, :levels] * delta[1])
                          * helem[elem_id, :levels])
            fluxes.append(float(flux))
        mismatch = int(max(ecav[elem_ids]) - min(ecav[elem_ids])) if elem_ids else 0
        records.append((max(abs(value) for value in fluxes), edge_id, n1, n2,
                        elem_ids, fluxes, mismatch))

    records.sort(reverse=True)
    print(f"Selected edges: {len(records)}")
    print("Approximate transport uses local tangent-plane edge metrics.")
    print("edge  nodes       elements       cavity-start(raw)  mismatch  flux contributions")
    for magnitude, edge_id, n1, n2, elem_ids, fluxes, mismatch in records[:args.limit]:
        raw_cavity = [int(ecav[e] + 1) for e in elem_ids]
        print(f"{edge_id + 1:5d}  {n1 + 1:6d},{n2 + 1:<6d} "
              f"{','.join(str(e + 1) for e in elem_ids):14s} "
              f"{','.join(map(str, raw_cavity)):17s} {mismatch:8d} "
              f"{','.join(f'{v:.6e}' for v in fluxes)}")

    print("\nNode ssh_rhs values:")
    for node in sorted(target):
        print(f"  node {node + 1}: lon={nodes[node,1]:.8f} "
              f"lat={nodes[node,2]:.8f} ssh_rhs={ssh_rhs[node]:.6e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
