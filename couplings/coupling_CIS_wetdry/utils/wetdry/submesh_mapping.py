#!/usr/bin/env python
"""
submesh_mapping: mapping from submesh to maxmesh — writes map_nod.out / map_elem.out
between a submesh and its parent maxmesh (QGLB).

The submesh is derived from the maxmesh via relax cut_off + SFC resort.
Since the submesh namelist uses no smoothing and no edge swapping,
node coordinates are preserved exactly. We exploit this by matching
Cartesian (x,y,z) coordinates via cKDTree on the unit sphere.

Usage:
    python submesh_mapping.py --maxmesh_dir <path> --submesh_dir <path>
"""

import os
import numpy as np
from scipy.spatial import cKDTree


Rearth = 6371.0  # km


def geo2cart(glon, glat):
    x = Rearth * np.cos(np.radians(glat)) * np.cos(np.radians(glon))
    y = Rearth * np.cos(np.radians(glat)) * np.sin(np.radians(glon))
    z = Rearth * np.sin(np.radians(glat))
    return x, y, z


def read_nod2d(mesh_dir):
    fpath = os.path.join(mesh_dir, "nod2d.out")
    with open(fpath) as f:
        n = int(f.readline().strip())
        data = np.loadtxt(f, max_rows=n)
    lon = data[:, 1]
    lat = data[:, 2]
    return n, lon, lat


def read_elem2d(mesh_dir):
    fpath = os.path.join(mesh_dir, "elem2d.out")
    with open(fpath) as f:
        n = int(f.readline().strip())
        data = np.loadtxt(f, dtype=int, max_rows=n)
    return n, data


def build_node_mapping(max_lon, max_lat, sub_lon, sub_lat, tol_km=0.1):
    """
    For each submesh node, find the nearest maxmesh node via cKDTree
    in 3D Cartesian space. Returns map_nod (1-based maxmesh indices).
    """
    max_x, max_y, max_z = geo2cart(max_lon, max_lat)
    sub_x, sub_y, sub_z = geo2cart(sub_lon, sub_lat)

    tree = cKDTree(np.column_stack([max_x, max_y, max_z]))
    dists, indices = tree.query(np.column_stack([sub_x, sub_y, sub_z]))

    max_dist = np.max(dists)
    mean_dist = np.mean(dists)
    print(f"  Node matching (Cartesian): max dist = {max_dist:.4e} km, "
          f"mean = {mean_dist:.4e} km")

    n_bad = np.sum(dists > tol_km)
    if n_bad > 0:
        print(f"  WARNING: {n_bad} nodes exceed tolerance {tol_km} km")
        worst = np.argsort(dists)[-min(5, n_bad):]
        for w in worst:
            print(f"      node {w+1}: dist={dists[w]:.6e} km")
    else:
        print(f"  All nodes within tolerance {tol_km} km. OK.")

    map_nod = indices + 1  # 1-based
    return map_nod


def build_elem_mapping(max_elem, sub_elem, map_nod):
    """
    For each submesh element, map its 3 nodes to maxmesh indices,
    then find the matching maxmesh element via node-set lookup.
    """
    n_max_elem = max_elem.shape[0]
    n_sub_elem = sub_elem.shape[0]

    max_elem_set = {}
    for i in range(n_max_elem):
        key = frozenset(max_elem[i, :].tolist())
        max_elem_set[key] = i + 1  # 1-based

    map_elem = np.zeros(n_sub_elem, dtype=int)
    n_matched = 0
    n_fallback = 0

    for i in range(n_sub_elem):
        sub_nodes = sub_elem[i, :]
        mapped_nodes = map_nod[sub_nodes - 1]  # sub_elem is 1-based
        key = frozenset(mapped_nodes.tolist())

        if key in max_elem_set:
            map_elem[i] = max_elem_set[key]
            n_matched += 1
        else:
            map_elem[i] = -1
            n_fallback += 1

    print(f"  Element matching: {n_matched}/{n_sub_elem} exact matches")
    if n_fallback > 0:
        print(f"  WARNING: {n_fallback} elements have no exact match "
              f"(likely due to edge swapping in maxmesh)")
        print(f"  Attempting centroid-based fallback...")
        map_elem = _fallback_elem_centroid(
            max_elem, map_nod, map_elem)

    return map_elem


def _fallback_elem_centroid(max_elem, map_nod, map_elem):
    """
    Fallback for unmatched elements: match by Cartesian centroid proximity.
    """
    unmatched = np.where(map_elem == -1)[0]
    if len(unmatched) == 0:
        return map_elem

    maxmesh_dir = os.environ.get("maxmesh_dir", ".")
    _, max_lon, max_lat = read_nod2d(maxmesh_dir)
    max_x, max_y, max_z = geo2cart(max_lon, max_lat)

    n_max_elem = max_elem.shape[0]
    cx = np.zeros(n_max_elem)
    cy = np.zeros(n_max_elem)
    cz = np.zeros(n_max_elem)
    for i in range(n_max_elem):
        nds = max_elem[i, :] - 1
        cx[i] = np.mean(max_x[nds])
        cy[i] = np.mean(max_y[nds])
        cz[i] = np.mean(max_z[nds])

    tree = cKDTree(np.column_stack([cx, cy, cz]))

    for idx in unmatched:
        mapped = map_nod[idx] - 1  # not used; need sub_elem
        # re-derive from stored info -- use maxmesh coords of mapped nodes
        # (this function is only called from build_elem_mapping context)
        pass

    # Simpler: re-read submesh and compute centroids
    submesh_dir = os.environ.get("submesh_dir", ".")
    _, sub_lon, sub_lat = read_nod2d(submesh_dir)
    sub_x, sub_y, sub_z = geo2cart(sub_lon, sub_lat)

    _, sub_elem_data = read_elem2d(submesh_dir)
    for idx in unmatched:
        nds = sub_elem_data[idx, :] - 1
        qx = np.mean(sub_x[nds])
        qy = np.mean(sub_y[nds])
        qz = np.mean(sub_z[nds])
        _, nearest = tree.query([qx, qy, qz])
        map_elem[idx] = nearest + 1

    print(f"  Fallback: assigned {len(unmatched)} elements via centroid matching")
    return map_elem


def write_mapping(fpath, mapping):
    with open(fpath, "w") as f:
        for val in mapping:
            f.write(f"{val}\n")


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Generate mesh mapping files")
    parser.add_argument("--maxmesh_dir", required=True,
                        help="Path to maxmesh (QGLB) directory")
    parser.add_argument("--submesh_dir", required=True,
                        help="Path to submesh directory")
    parser.add_argument("--tol", type=float, default=0.1,
                        help="Coordinate matching tolerance in km (default: 0.1)")

    args = parser.parse_args()
    maxmesh_dir = args.maxmesh_dir
    submesh_dir = args.submesh_dir

    os.environ["maxmesh_dir"] = maxmesh_dir
    os.environ["submesh_dir"] = submesh_dir

    print(f"Generating mesh mapping:")
    print(f"  maxmesh: {maxmesh_dir}")
    print(f"  submesh: {submesh_dir}")

    # Read meshes
    n_max, max_lon, max_lat = read_nod2d(maxmesh_dir)
    n_sub, sub_lon, sub_lat = read_nod2d(submesh_dir)
    print(f"  maxmesh: {n_max} nodes")
    print(f"  submesh: {n_sub} nodes")

    n_max_e, max_elem = read_elem2d(maxmesh_dir)
    n_sub_e, sub_elem = read_elem2d(submesh_dir)
    print(f"  maxmesh: {n_max_e} elements")
    print(f"  submesh: {n_sub_e} elements")

    # Build node mapping in Cartesian space
    print("\nBuilding node mapping...")
    map_nod = build_node_mapping(max_lon, max_lat, sub_lon, sub_lat,
                                 tol_km=args.tol)

    # Validate uniqueness
    unique_mapped = len(np.unique(map_nod))
    if unique_mapped != n_sub:
        print(f"  WARNING: {n_sub - unique_mapped} duplicate mappings detected!")
    else:
        print(f"  All {n_sub} node mappings are unique. OK.")

    # Build element mapping
    print("\nBuilding element mapping...")
    map_elem = build_elem_mapping(max_elem, sub_elem, map_nod)

    # Write output
    out_nod = os.path.join(submesh_dir, "map_nod.out")
    out_elem = os.path.join(submesh_dir, "map_elem.out")
    write_mapping(out_nod, map_nod)
    write_mapping(out_elem, map_elem)
    print(f"\nWritten: {out_nod}")
    print(f"Written: {out_elem}")
    print("Done.")


if __name__ == "__main__":
    main()
