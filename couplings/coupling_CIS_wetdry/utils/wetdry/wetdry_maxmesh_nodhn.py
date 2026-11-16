#!/usr/bin/env python3
"""
Interpolate PISM fields onto FESOM QGLB/maxmesh nodes and build nodhn for relax cutting.

Input (couple_out ice_wetdry_cut_fields.nc):
  - shelf_base_ut: usurf - thk on the full PISM grid (m)
  - grounded_ice_mask: 1 where PISM mask == 2 (grounded ice), else 0
  - floating_ice_mask: 1 where PISM mask == 3 (floating ice shelf), else 0

  nodhn = grounded_nodhn_value (default -9999) where grounded mask > 0.5 after NN,
          else (shelf_base_ut - bedrock).
wetdry_shelf_on_maxmesh.out: interpolated shelf_base_ut, set to 0 on grounded nodes for cavity
mapping.

Continuous fields are transferred with mask-aware bilinear interpolation.  The
PISM latitude axis is sorted before interpolation, so descending latitude
(which can appear vertically inverted in ncview) is handled by coordinates,
not by an unsafe visual flip.  PISM files may store lon/lat as 2-D arrays;
rectilinear 2-D coordinates are reduced to their x/y axes before interpolation.
Masks remain nearest-neighbour/categorical.

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


def _read_grid_variable(v, time_index=-1):
    dims = v.dimensions
    if "time" not in dims:
        return np.asarray(v[:], dtype=np.float64)
    idx = [slice(None)] * len(dims)
    idx[list(dims).index("time")] = time_index
    return np.asarray(v[tuple(idx)], dtype=np.float64)


def _normalise_structured_field(data, lon, lat):
    """Return ascending lon/lat axes and data indexed as [lat, lon]."""
    if lon.ndim != 1 or lat.ndim != 1 or data.ndim != 2:
        raise ValueError("expected 1-D lon/lat and a 2-D PISM field")
    if data.shape == (lat.size, lon.size):
        field = data
    elif data.shape == (lon.size, lat.size):
        field = data.T
    else:
        raise ValueError(
            f"field shape {data.shape} incompatible with lat/lon sizes "
            f"{lat.size}/{lon.size}"
        )
    ilon = np.argsort(lon)
    ilat = np.argsort(lat)
    return lon[ilon], lat[ilat], field[np.ix_(ilat, ilon)]


def _target_lon(lon, source_lon):
    out = np.asarray(lon, dtype=np.float64).copy()
    if np.nanmin(source_lon) >= 0.0 and np.nanmin(out) < 0.0:
        out[out < 0.0] += 360.0
    elif np.nanmax(source_lon) <= 180.0 and np.nanmax(out) > 180.0:
        out[out > 180.0] -= 360.0
    return out


def _grid_lonlat(lon, lat):
    if np.asarray(lon).ndim == 2 and np.asarray(lat).ndim == 2:
        return np.asarray(lon), np.asarray(lat)
    return np.meshgrid(lon, lat)


def bilinear_grid_sample(lon, lat, field, target_lon, target_lat, valid=None):
    """Bilinear sample with valid-corner weight renormalisation."""
    if np.asarray(lon).ndim == 2:
        # PISM polar grids often have genuinely curvilinear lon/lat.  Use a
        # small local, distance-weighted stencil in Cartesian-on-the-sphere
        # coordinates.  This is the curvilinear-grid equivalent of bilinear
        # interpolation, avoids a very large global Delaunay triangulation,
        # and still supports valid-corner renormalisation for the shelf mask.
        source_lon = np.asarray(lon).ravel()
        source_lat = np.asarray(lat).ravel()
        source_field = np.asarray(field).ravel()
        good = np.isfinite(source_field)
        if valid is not None:
            good &= np.asarray(valid).ravel() > 0.5
        if np.count_nonzero(good) < 3:
            return np.full(len(np.asarray(target_lon)), np.nan)
        sx, sy, sz = geo2cart(source_lon[good], source_lat[good])
        tx, ty, tz = geo2cart(_target_lon(target_lon, source_lon), target_lat)
        tree = cKDTree(np.column_stack((sx, sy, sz)))
        k = min(8, np.count_nonzero(good))
        dist, idx = tree.query(np.column_stack((tx, ty, tz)), k=k)
        if k == 1:
            dist = dist[:, None]
            idx = idx[:, None]
        weights = 1.0 / np.maximum(dist, 1.0e-12)
        values = source_field[good][idx]
        exact = dist[:, 0] <= 1.0e-12
        out = np.sum(weights * values, axis=1) / np.sum(weights, axis=1)
        out[exact] = values[exact, 0]
        return out
    x = _target_lon(target_lon, lon)
    y = np.asarray(target_lat, dtype=np.float64)
    ix = np.clip(np.searchsorted(lon, x, side="right") - 1, 0, len(lon) - 2)
    iy = np.clip(np.searchsorted(lat, y, side="right") - 1, 0, len(lat) - 2)
    x0, x1 = lon[ix], lon[ix + 1]
    y0, y1 = lat[iy], lat[iy + 1]
    tx = np.divide(x - x0, x1 - x0, out=np.zeros_like(x), where=(x1 != x0))
    ty = np.divide(y - y0, y1 - y0, out=np.zeros_like(y), where=(y1 != y0))
    weights = np.column_stack(((1 - tx) * (1 - ty), tx * (1 - ty),
                               (1 - tx) * ty, tx * ty))
    values = np.column_stack((field[iy, ix], field[iy, ix + 1],
                              field[iy + 1, ix], field[iy + 1, ix + 1]))
    good = np.isfinite(values)
    if valid is not None:
        good &= np.column_stack((valid[iy, ix], valid[iy, ix + 1],
                                 valid[iy + 1, ix], valid[iy + 1, ix + 1])) > 0.5
    weights = np.where(good, weights, 0.0)
    denom = weights.sum(axis=1)
    out = np.full(x.shape, np.nan, dtype=np.float64)
    have = denom > 0.0
    out[have] = np.sum(weights[have] * values[have], axis=1) / denom[have]
    return out


def nearest_grid_sample(lon, lat, field, target_lon, target_lat, valid=None):
    xx, yy = _grid_lonlat(lon, lat)
    good = np.isfinite(field)
    if valid is not None:
        good &= valid > 0.5
    good = good.ravel()
    if not np.any(good):
        return np.full(len(np.asarray(target_lon)), np.nan)
    if np.asarray(lon).ndim == 2:
        sx, sy, sz = geo2cart(xx.ravel()[good], yy.ravel()[good])
        tx, ty, tz = geo2cart(target_lon, target_lat)
        tree = cKDTree(np.column_stack((sx, sy, sz)))
        _, idx = tree.query(np.column_stack((tx, ty, tz)))
    else:
        tx = _target_lon(target_lon, np.asarray(lon).ravel())
        query = np.column_stack((tx, np.asarray(target_lat, dtype=np.float64)))
        tree = cKDTree(np.column_stack((xx.ravel()[good], yy.ravel()[good])))
        _, idx = tree.query(query)
    return field.ravel()[good][idx]


def read_elem2d(path):
    with open(path) as f:
        n = int(f.readline().strip())
        return np.loadtxt(f, dtype=np.int64, max_rows=n).reshape((-1, 3)) - 1


def smooth_active_field(values, active, elem, passes=1, blend=0.5):
    """Graph smoothing restricted to floating-shelf nodes."""
    if passes <= 0 or blend <= 0.0:
        return values.copy()
    neighbours = [[] for _ in range(len(values))]
    for tri in elem:
        for i, j in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
            neighbours[i].append(j)
            neighbours[j].append(i)
    out = values.copy()
    for _ in range(passes):
        old = out.copy()
        for i in np.where(active)[0]:
            nb = np.asarray([j for j in neighbours[i]
                             if active[j] and np.isfinite(old[j])], dtype=int)
            if nb.size:
                out[i] = (1.0 - blend) * old[i] + blend * np.mean(np.r_[old[i], old[nb]])
    return out


def read_pism_cut_grid(nc_path, shelf_var, grounded_var, floating_var=None, time_index=-1):
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

        shelf = _read_grid_variable(nc.variables[shelf_var], time_index)
        grounded = _read_grid_variable(nc.variables[grounded_var], time_index)
        floating = None
        if floating_var and floating_var in nc.variables:
            floating = _read_grid_variable(nc.variables[floating_var], time_index)

    # PISM commonly stores lon/lat as 2-D (y,x) coordinate variables even when
    # the grid is rectilinear.  Reduce that representation to the corresponding
    # axes, retaining the data orientation explicitly.  Do not transpose or
    # flip based on ncview appearance.
    curvilinear = lon1.ndim == 2 or lat1.ndim == 2
    if curvilinear:
        if lon1.ndim != 2 or lat1.ndim != 2 or lon1.shape != lat1.shape:
            raise ValueError("lon/lat must both be 2-D with identical shape")
        if shelf.ndim != 2 or shelf.shape != lon1.shape:
            raise ValueError("2-D lon/lat shape does not match PISM shelf field")
        lon_axis = lon1[0, :]
        lat_axis = lat1[:, 0]
        lon_rect = np.broadcast_to(lon_axis[None, :], lon1.shape)
        lat_rect = np.broadcast_to(lat_axis[:, None], lat1.shape)
        if not (np.allclose(lon1, lon_rect, equal_nan=True)
                and np.allclose(lat1, lat_rect, equal_nan=True)):
            # Keep genuinely curvilinear coordinates for LinearND interpolation.
            return lon1, lat1, shelf, grounded, floating, "curvilinear"
        lon1, lat1 = np.asarray(lon_axis), np.asarray(lat_axis)
        curvilinear = False

    lat_descending = bool(lat1.size > 1 and lat1[0] > lat1[-1])
    lon, lat, shelf = _normalise_structured_field(shelf, lon1, lat1)
    _, _, grounded = _normalise_structured_field(grounded, lon1, lat1)
    if floating is not None:
        _, _, floating = _normalise_structured_field(floating, lon1, lat1)
    return lon, lat, shelf, grounded, floating, lat_descending


def read_pism_anomaly(nc_path, anomaly_var="topg_anomaly", time_index=-1):
    """Read a PISM-grid bedrock anomaly using the same coordinate layout."""
    with Dataset(nc_path, "r") as nc:
        if anomaly_var not in nc.variables:
            raise KeyError(f"variable {anomaly_var!r} not in {nc_path}")
        lon_nm = _find_lon_lat_var(nc, ("lon", "longitude", "lon_2", "x"))
        lat_nm = _find_lon_lat_var(nc, ("lat", "latitude", "lat_2", "y"))
        if lon_nm is None or lat_nm is None:
            raise KeyError(f"need lon/lat coordinates in {nc_path}")
        lon1 = np.asarray(nc.variables[lon_nm][:], dtype=np.float64)
        lat1 = np.asarray(nc.variables[lat_nm][:], dtype=np.float64)
        anomaly = _read_grid_variable(nc.variables[anomaly_var], time_index)
    curvilinear = lon1.ndim == 2 or lat1.ndim == 2
    if lon1.ndim == 2 or lat1.ndim == 2:
        if lon1.ndim != 2 or lat1.ndim != 2 or lon1.shape != lat1.shape:
            raise ValueError("anomaly lon/lat must both be 2-D with identical shape")
        return lon1, lat1, anomaly
    return _normalise_structured_field(anomaly, lon1, lat1)


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
    p.add_argument("--floating_mask_var", default="floating_ice_mask")
    p.add_argument("--bedrock_file", default="bedrock_topography.out", help="file in qglb_global")
    p.add_argument(
        "--topg_anomaly_nc",
        default=None,
        help="optional PISM-grid topg anomaly; mapped with the selected continuous interpolator",
    )
    p.add_argument(
        "--out_bedrock",
        default=None,
        help="optional FESOM-node bedrock file (QGLB bedrock plus mapped anomaly)",
    )
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
        help="nodhn at nodes whose nearest PISM cell is grounded ice (mask==2)",
    )
    p.add_argument(
        "--grounded_mask_threshold",
        type=float,
        default=0.5,
        help="treat NN grounded_ice_mask > this as grounded",
    )
    p.add_argument("--max_search_km", type=float, default=800.0)
    p.add_argument("--smooth_passes", type=int, default=1,
                   help="floating-shelf graph smoothing passes (default 1; 0 disables)")
    p.add_argument("--smooth_blend", type=float, default=0.5,
                   help="blend toward local floating-shelf neighbor mean (default 0.5)")
    p.add_argument("--interpolation", choices=("bilinear", "nearest"), default="bilinear",
                   help="continuous PISM-field interpolation method (default bilinear)")
    args = p.parse_args()

    gdir = args.qglb_global.rstrip("/")
    n_nod, tlon, tlat = read_nod2d(gdir, args.nod2d)
    bedrock = read_bedrock(gdir, args.bedrock_file)
    if bedrock.shape[0] != n_nod:
        raise ValueError(
            f"bedrock length {bedrock.shape[0]} != nod2d count {n_nod} in {gdir}"
        )

    plon, plat, pshelf, pground, pfloating, lat_descending = read_pism_cut_grid(
        args.pism_cut_nc, args.shelf_ut_var, args.grounded_mask_var,
        args.floating_mask_var
    )

    # Backward compatibility for old cut files without floating_ice_mask.
    if pfloating is None:
        pfloating = ((pground <= args.grounded_mask_threshold)
                     & np.isfinite(pshelf) & (np.abs(pshelf) > 1.0e-8)).astype(float)
        print("wetdry_maxmesh_nodhn: floating_ice_mask absent; inferred from nonzero shelf_base_ut")
    print(
        f"wetdry_maxmesh_nodhn: PISM grid nx={plon.size} ny={plat.size} "
        f"latitude_input={lat_descending if isinstance(lat_descending, str) else ('descending' if lat_descending else 'ascending')}; "
        "interpolation uses coordinate-aware grid handling"
    )

    # Nearest neighbour is retained only for categorical mask classification
    # and for distance-limited fallbacks.  Continuous fields use bilinear
    # interpolation in the PISM coordinate system.
    pgrid_lon, pgrid_lat = _grid_lonlat(plon, plat)
    px, py, pz = geo2cart(pgrid_lon, pgrid_lat)
    tree = cKDTree(np.column_stack([px, py, pz]))
    tx, ty, tz = geo2cart(tlon, tlat)
    dists, idx = tree.query(np.column_stack([tx, ty, tz]))

    grounded_nn = pground.ravel()[idx].copy()
    floating_nn = pfloating.ravel()[idx].copy()
    grounded_nn[dists > args.max_search_km] = 0.0
    floating_nn[dists > args.max_search_km] = 0.0

    if args.interpolation == "bilinear":
        shelf_max = bilinear_grid_sample(
            plon, plat, pshelf, tlon, tlat, valid=pfloating
        )
    else:
        shelf_max = nearest_grid_sample(plon, plat, pshelf, tlon, tlat, valid=pfloating)
    missing = ~np.isfinite(shelf_max)
    if np.any(missing):
        shelf_max[missing] = nearest_grid_sample(
            plon, plat, pshelf, tlon[missing], tlat[missing], valid=pfloating
        )
    shelf_max[~np.isfinite(shelf_max)] = 0.0
    shelf_max[dists > args.max_search_km] = 0.0
    # Do not create a shelf cavity at nodes classified as open ocean/grounded.
    shelf_max[floating_nn <= 0.5] = 0.0

    raw_shelf_max = shelf_max.copy()
    if args.smooth_passes > 0:
        elem_path = os.path.join(gdir, "elem2d.out")
        if os.path.isfile(elem_path):
            elem = read_elem2d(elem_path)
            active = floating_nn > 0.5
            shelf_max = smooth_active_field(
                shelf_max, active, elem, args.smooth_passes, args.smooth_blend
            )
            shelf_max[~active] = 0.0
        else:
            print(f"wetdry_maxmesh_nodhn: no {elem_path}; smoothing skipped")

    bedrock_used = bedrock.copy()
    anomaly_nodes = np.zeros(n_nod, dtype=np.float64)
    if args.topg_anomaly_nc:
        alon, alat, anomaly = read_pism_anomaly(args.topg_anomaly_nc)
        if args.interpolation == "bilinear":
            anomaly_nodes = bilinear_grid_sample(alon, alat, anomaly, tlon, tlat)
        else:
            anomaly_nodes = nearest_grid_sample(alon, alat, anomaly, tlon, tlat)
        missing = ~np.isfinite(anomaly_nodes)
        if np.any(missing):
            anomaly_nodes[missing] = nearest_grid_sample(
                alon, alat, anomaly, tlon[missing], tlat[missing]
            )
        anomaly_nodes = np.where(np.isfinite(anomaly_nodes), anomaly_nodes, 0.0)
        agrid_lon, agrid_lat = _grid_lonlat(alon, alat)
        ax, ay, az = geo2cart(agrid_lon, agrid_lat)
        atree = cKDTree(np.column_stack([ax, ay, az]))
        adists, _ = atree.query(np.column_stack([tx, ty, tz]))
        anomaly_nodes[adists > args.max_search_km] = 0.0
        bedrock_used += anomaly_nodes
        print(
            "wetdry_maxmesh_nodhn: mapped bedrock anomaly "
            f"min={np.nanmin(anomaly_nodes):.6e} "
            f"max={np.nanmax(anomaly_nodes):.6e} "
            f"mean={np.nanmean(anomaly_nodes):.6e}"
        )
    if args.out_bedrock:
        odir = os.path.dirname(os.path.abspath(args.out_bedrock))
        if odir:
            os.makedirs(odir, exist_ok=True)
        np.savetxt(args.out_bedrock, bedrock_used, fmt="%.8e")

    nodhn = shelf_max - bedrock_used
    grounded_here = grounded_nn > args.grounded_mask_threshold
    nodhn = np.where(grounded_here, args.grounded_nodhn_value, nodhn)

    # cavity mapping: no shelf cavity over grounded ice sheet
    shelf_for_cavity = np.where(grounded_here, 0.0, shelf_max)

    odir = os.path.dirname(os.path.abspath(args.out_nodhn))
    if odir:
        os.makedirs(odir, exist_ok=True)
    np.savetxt(args.out_nodhn, nodhn, fmt="%.8e")
    np.savetxt(args.out_shelf_max, shelf_for_cavity, fmt="%.8e")
    shelf_dir = os.path.dirname(os.path.abspath(args.out_shelf_max))
    np.savetxt(os.path.join(shelf_dir, "wetdry_shelf_on_maxmesh_raw.out"), raw_shelf_max, fmt="%.8e")
    np.savetxt(os.path.join(shelf_dir, "wetdry_shelf_on_maxmesh_smoothed.out"), shelf_max, fmt="%.8e")
    with open(os.path.join(shelf_dir, "wetdry_geometry_interp_diagnostics.txt"), "w") as diag:
        diag.write(f"interpolation={args.interpolation}\n")
        diag.write("pism_latitude_input=" + (lat_descending if isinstance(lat_descending, str) else ('descending' if lat_descending else 'ascending')) + "\n")
        diag.write(f"smooth_passes={args.smooth_passes}\n")
        diag.write(f"smooth_blend={args.smooth_blend}\n")
        diag.write(f"floating_nodes={int(np.count_nonzero(floating_nn > 0.5))}\n")
        diag.write(f"raw_shelf_min={np.nanmin(raw_shelf_max):.8e}\n")
        diag.write(f"raw_shelf_max={np.nanmax(raw_shelf_max):.8e}\n")
        diag.write(f"used_shelf_min={np.nanmin(shelf_max):.8e}\n")
        diag.write(f"used_shelf_max={np.nanmax(shelf_max):.8e}\n")

    print(
        f"wetdry_maxmesh_nodhn: mode=cut_fields interpolation={args.interpolation} nodes={n_nod}  wrote {args.out_nodhn} "
        f"(grounded_nodhn_value={args.grounded_nodhn_value})"
    )
    print(f"wetdry_maxmesh_nodhn: wrote {args.out_shelf_max} (shelf for cavity, 0 on grounded)")


if __name__ == "__main__":
    main()
