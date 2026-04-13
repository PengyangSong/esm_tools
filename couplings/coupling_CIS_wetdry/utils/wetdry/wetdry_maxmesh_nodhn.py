#!/usr/bin/env python3
"""
Interpolate PISM fields onto FESOM QGLB/maxmesh nodes and build nodhn for relax cutting.

Input (couple_out ice_wetdry_cut_fields.nc):
  - shelf_base_ut: usurf - thk on the full PISM grid (m)
  - grounded_ice_mask: 1 where PISM mask == 2 (grounded ice), else 0

  nodhn = grounded_nodhn_value (default -9999) where grounded mask > 0.5 after NN,
          else (shelf_base_ut - bedrock).
wetdry_shelf_on_maxmesh.out: interpolated shelf_base_ut, set to 0 on grounded nodes for cavity
mapping.

Requires: numpy, scipy, netCDF4
"""

from __future__ import annotations

import argparse
import os

import numpy as np
from netCDF4 import Dataset
from scipy.spatial import cKDTree

Rearth = 6371.0  # km


def geo2cart(glon, glat):
    glon = np.asarray(glon, dtype=np.float64)
    glat = np.asarray(glat, dtype=np.float64)
    x = Rearth * np.cos(np.radians(glat)) * np.cos(np.radians(glon))
    y = Rearth * np.cos(np.radians(glat)) * np.sin(np.radians(glon))
    z = Rearth * np.sin(np.radians(glat))
    return x.ravel(), y.ravel(), z.ravel()


def read_nod2d(mesh_dir, nod2d_basename="nod2d.out"):
    fpath = os.path.join(mesh_dir, nod2d_basename)
    with open(fpath) as f:
        n = int(f.readline().strip())
        data = np.loadtxt(f, max_rows=n)
    lon = np.asarray(data[:, 1], dtype=np.float64)
    lat = np.asarray(data[:, 2], dtype=np.float64)
    return n, lon, lat


def read_bedrock(global_dir, fname):
    path = os.path.join(global_dir, fname)
    if not os.path.isfile(path):
        raise FileNotFoundError(f"bedrock file not found: {path}")
    b = np.loadtxt(path, dtype=np.float64)
    if b.ndim > 1:
        b = b[:, 0]
    return b


def _find_lon_lat_var(nc, names):
    for nm in names:
        if nm in nc.variables:
            return nm
    return None


def _read_field_for_lonlat_grid(v, time_index=-1):
    """
    Return a 2D (y, x) field for alignment with lon/lat.

    If the variable has a 'time' dimension, take the requested time slice (default
    last). If there is no time dimension (e.g. ice_file_at_ocean.combined.nc after
    ncwa), use the full 2D field. Indexing v[time_index] on a 2D variable would
    wrongly slice the first spatial dimension as if it were time.
    """
    dims = v.dimensions
    if "time" not in dims:
        return np.asarray(v[:], dtype=np.float64)
    if dims[0] == "time":
        return np.asarray(v[time_index], dtype=np.float64)
    ti = list(dims).index("time")
    idx = [slice(None)] * len(dims)
    idx[ti] = time_index
    return np.asarray(v[tuple(idx)], dtype=np.float64)


def _align_lon_lat_data(data, lon1, lat1):
    if data.ndim == 2 and lon1.ndim == 1 and lat1.ndim == 1:
        ny, nx = data.shape
        if lon1.size == nx and lat1.size == ny:
            lon2, lat2 = np.meshgrid(lon1, lat1)
            lon = lon2
            lat = lat2
        else:
            raise ValueError(
                f"expected lon length {nx} and lat length {ny} for field {data.shape}, "
                f"got lon {lon1.size} lat {lat1.size}"
            )
    else:
        lon = lon1
        lat = lat1

    if lon.shape != data.shape or lat.shape != data.shape:
        raise ValueError(
            f"cannot align lon {lon.shape} lat {lat.shape} with data {data.shape}"
        )

    vals = np.asarray(data, dtype=np.float64).ravel()
    lon = np.asarray(lon, dtype=np.float64).ravel()
    lat = np.asarray(lat, dtype=np.float64).ravel()
    return lon, lat, vals


def read_pism_cut_pair(nc_path, shelf_var, grounded_var, time_index=-1):
    with Dataset(nc_path, "r") as nc:
        if shelf_var not in nc.variables:
            raise KeyError(f"variable {shelf_var!r} not in {nc_path}")
        if grounded_var not in nc.variables:
            raise KeyError(f"variable {grounded_var!r} not in {nc_path}")

        lon_nm = _find_lon_lat_var(nc, ("lon", "longitude", "lon_2", "x"))
        lat_nm = _find_lon_lat_var(nc, ("lat", "latitude", "lat_2", "y"))
        if lon_nm is None or lat_nm is None:
            raise KeyError(f"need lon/lat (or longitude/latitude or x/y) in {nc_path}")

        lon1 = np.asarray(nc.variables[lon_nm][:], dtype=np.float64)
        lat1 = np.asarray(nc.variables[lat_nm][:], dtype=np.float64)

        s = _read_field_for_lonlat_grid(nc.variables[shelf_var], time_index)
        g = _read_field_for_lonlat_grid(nc.variables[grounded_var], time_index)

    lon_s, lat_s, shelf_r = _align_lon_lat_data(s, lon1, lat1)
    lon_g, lat_g, ground_r = _align_lon_lat_data(g, lon1, lat1)
    if not (np.allclose(lon_s, lon_g) and np.allclose(lat_s, lat_g)):
        raise ValueError("shelf and grounded fields have different lon/lat layout")

    ok = (
        np.isfinite(shelf_r)
        & np.isfinite(ground_r)
        & np.isfinite(lon_s)
        & np.isfinite(lat_s)
    )
    return lon_s[ok], lat_s[ok], shelf_r[ok], ground_r[ok]


def main():
    p = argparse.ArgumentParser(description="PISM shelf on maxmesh + nodhn for relax")
    p.add_argument("--qglb_global", required=True, help="QGLB global dir (nod2d.out, bedrock)")
    p.add_argument(
        "--pism_cut_nc",
        required=True,
        help="ice_wetdry_cut_fields.nc: shelf_base_ut + grounded_ice_mask",
    )
    p.add_argument("--shelf_ut_var", default="shelf_base_ut")
    p.add_argument("--grounded_mask_var", default="grounded_ice_mask")
    p.add_argument("--bedrock_file", default="bedrock_topography.out", help="file in qglb_global")
    p.add_argument(
        "--nod2d",
        default="nod2d.out",
        help="nod2d file in qglb_global (must match maxmesh; default nod2d.out)",
    )
    p.add_argument("--out_nodhn", required=True, help="output nodhn.out path")
    p.add_argument("--out_shelf_max", required=True, help="output shelf on maxmesh (ascii)")
    p.add_argument(
        "--grounded_nodhn_value",
        type=float,
        default=-9999.0,
        help="nodhn at nodes whose NN PISM cell is grounded ice (mask==2)",
    )
    p.add_argument(
        "--grounded_mask_threshold",
        type=float,
        default=0.5,
        help="treat NN grounded_ice_mask > this as grounded",
    )
    p.add_argument("--max_search_km", type=float, default=800.0)
    args = p.parse_args()

    gdir = args.qglb_global.rstrip("/")
    n_nod, tlon, tlat = read_nod2d(gdir, args.nod2d)
    bedrock = read_bedrock(gdir, args.bedrock_file)
    if bedrock.shape[0] != n_nod:
        raise ValueError(
            f"bedrock length {bedrock.shape[0]} != nod2d count {n_nod} in {gdir}"
        )

    plon, plat, pshelf, pground = read_pism_cut_pair(
        args.pism_cut_nc, args.shelf_ut_var, args.grounded_mask_var
    )

    px, py, pz = geo2cart(plon, plat)
    tree = cKDTree(np.column_stack([px, py, pz]))
    tx, ty, tz = geo2cart(tlon, tlat)
    dists, idx = tree.query(np.column_stack([tx, ty, tz]))

    shelf_max = pshelf[idx].copy()
    grounded_nn = pground[idx].copy()
    shelf_max[dists > args.max_search_km] = 0.0
    grounded_nn[dists > args.max_search_km] = 0.0

    nodhn = shelf_max - bedrock
    grounded_here = grounded_nn > args.grounded_mask_threshold
    nodhn = np.where(grounded_here, args.grounded_nodhn_value, nodhn)

    # cavity mapping: no shelf cavity over grounded ice sheet
    shelf_for_cavity = np.where(grounded_here, 0.0, shelf_max)

    odir = os.path.dirname(os.path.abspath(args.out_nodhn))
    if odir:
        os.makedirs(odir, exist_ok=True)
    np.savetxt(args.out_nodhn, nodhn, fmt="%.8e")
    np.savetxt(args.out_shelf_max, shelf_for_cavity, fmt="%.8e")

    print(
        f"wetdry_maxmesh_nodhn: mode=cut_fields nodes={n_nod}  wrote {args.out_nodhn} "
        f"(grounded_nodhn_value={args.grounded_nodhn_value})"
    )
    print(f"wetdry_maxmesh_nodhn: wrote {args.out_shelf_max} (shelf for cavity, 0 on grounded)")


if __name__ == "__main__":
    main()
