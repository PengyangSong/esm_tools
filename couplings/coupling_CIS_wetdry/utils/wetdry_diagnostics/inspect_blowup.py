#!/usr/bin/env python3
"""Inspect a FESOM blow-up file without loading unnecessary model fields."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


def read_nodes(mesh_dir: Path):
    node_file = mesh_dir / "nod2d.out"
    with node_file.open() as handle:
        count = int(handle.readline())
    nodes = np.loadtxt(node_file, skiprows=1, max_rows=count)
    return nodes


def read_elements(mesh_dir: Path):
    elem_file = mesh_dir / "elem2d.out"
    with elem_file.open() as handle:
        count = int(handle.readline())
    return np.loadtxt(elem_file, skiprows=1, max_rows=count, dtype=int) - 1


def neighbours(elements: np.ndarray, node: int) -> set[int]:
    result: set[int] = set()
    for triangle in elements:
        if node in triangle:
            result.update(int(item) for item in triangle)
    return result


def field(nc: Dataset, name: str) -> np.ndarray | None:
    variable = nc.variables.get(name)
    if variable is None:
        return None
    return np.asarray(variable[0])


def scalar_summary(name: str, values: np.ndarray) -> str:
    finite = np.isfinite(values)
    if not finite.any():
        return f"{name}: no finite values"
    finite_values = values[finite]
    return (
        f"{name}: min={finite_values.min():.6e} "
        f"max={finite_values.max():.6e} "
        f"absmax={np.max(np.abs(finite_values)):.6e} "
        f"nonfinite={values.size - finite.sum()}"
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("blowup", type=Path)
    parser.add_argument("mesh", type=Path)
    parser.add_argument("--radius", type=int, default=1)
    parser.add_argument("--threshold", type=float, default=1.0e6)
    args = parser.parse_args()

    nodes = read_nodes(args.mesh)
    elements = read_elements(args.mesh)
    with Dataset(args.blowup) as nc:
        print("Dimensions:", ", ".join(f"{k}={len(v)}" for k, v in nc.dimensions.items()))
        print("Global field summaries:")
        for name in (
            "ssh_rhs", "ssh_rhs_old", "d_eta", "hbar", "eta_n",
            "temp", "salt", "hnode", "w", "cfl_z",
        ):
            values = field(nc, name)
            if values is not None:
                print("  " + scalar_summary(name, values))

        temp = field(nc, "temp")
        if temp is None:
            raise SystemExit("The blow-up file has no temp variable")
        # FESOM deliberately stores NaN in inactive/cavity/below-bed tracer
        # cells.  Do not classify those masked cells as unstable; classify a
        # node only when a finite tracer value exceeds the threshold.
        bad = np.nanmax(np.abs(temp), axis=-1) > args.threshold
        bad_nodes = np.flatnonzero(bad)
        print(f"\nUnstable nodes (threshold {args.threshold:.3e}): {len(bad_nodes)}")

        values = {name: field(nc, name) for name in (
            "ssh_rhs", "ssh_rhs_old", "d_eta", "hbar", "eta_n", "salt",
            "hnode", "w", "cfl_z",
        )}

        for node in bad_nodes:
            nearby = neighbours(elements, int(node))
            if args.radius > 1:
                frontier = set(nearby)
                for _ in range(args.radius - 1):
                    frontier |= {n for item in frontier for n in neighbours(elements, item)}
                nearby = frontier
            print(
                f"\nnode={node + 1} lon={nodes[node, 1]:.8f} "
                f"lat={nodes[node, 2]:.8f} neighbors={len(nearby)}"
            )
            for name, array in values.items():
                if array is None:
                    continue
                column = array[node]
                if column.ndim == 0:
                    print(f"  {name}={float(column):.8e}")
                else:
                    finite = np.isfinite(column)
                    print(
                        f"  {name}: min={np.nanmin(column):.6e} "
                        f"max={np.nanmax(column):.6e} "
                        f"nonfinite={column.size - finite.sum()}"
                    )

            neighbor_bad = []
            for item in sorted(nearby):
                if item == node:
                    continue
                ntemp = temp[item]
                if np.nanmax(np.abs(ntemp)) > args.threshold:
                    neighbor_bad.append(item + 1)
            print("  unstable one-ring nodes:", neighbor_bad or "none")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
