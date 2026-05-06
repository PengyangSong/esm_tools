# -*- coding: utf-8 -*-
"""
Map a FESOM variable to column points with longitude/latitude metadata.

For cavity/interface variables this writes a true 2D field with dimensions
``(time, horizontal)`` because downstream PISM ``-ocean given`` / ``-ocean th``
forcing reads 2D time-dependent fields and does not use vertical coordinates.

For multi-level variables this keeps the legacy 3D output
``(time, level, horizontal)`` with a ``depth`` coordinate.
"""

import argparse
import sys
import time

import numpy as np
from netCDF4 import Dataset
import pyfesom2 as pf
import xarray as xr

sys.path.append("./")
sys.path.append("./pyfesom")

NAN_REPLACE = 9.96921e36


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--FESOM_PATH", nargs="*", required=True)
    parser.add_argument("--FESOM_VARIABLE", nargs="+")
    parser.add_argument("--FESOM_MESH", nargs="+")
    parser.add_argument("--FESOM_YEARS", nargs="+")
    parser.add_argument("--FESOM_OUTPUT", nargs="*", default="fesom_output.nc4")
    parser.add_argument("--FESOM_MESH_ROTATED", nargs="*", required=True)
    return parser.parse_args()


def detect_surface_layer(dataset, varname, mesh):
    """
    True for native 2D fields (time + horizontal) or a 3D field with a singleton
    vertical/surface dimension.

    Cavity/interface variables are often (time, nod2d) with no depth axis; the
    legacy check (any non-time dim of size 1) misses that and wrongly routes
    them through read_3d_field, producing bogus vertical levels from mesh.zlev.

    Horizontal size is matched against mesh nodal counts from pyfesom2 (typically
    ``mesh.n2d``); output lon/lat are nodal (``mesh.x2`` / ``mesh.y2``).
    """
    var = dataset[varname]
    if var.ndim not in (2, 3):
        return False

    time_dim = "time" if "time" in var.dims else "T"
    non_time_dims = [d for d in var.dims if d != time_dim]
    non_time_sizes = [int(dataset.sizes[d]) for d in non_time_dims]
    if any(size == 1 for size in non_time_sizes):
        return True
    nodal_sizes = {
        int(x)
        for x in (
            getattr(mesh, "n2d", None),
            getattr(mesh, "nod2d", None),
        )
        if x is not None
    }
    if var.ndim == 2 and nodal_sizes:
        if any(int(dataset.sizes[d]) in nodal_sizes for d in non_time_dims):
            return True
    return False


def load_mesh(mesh_path, mesh_rotated):
    print("* Read the mesh MESHPATH=" + " ".join(mesh_path))
    euler_angle = [0, 0, 0] if mesh_rotated else [50, 15, -90]
    print("* Using Euler angles ", euler_angle)
    mesh = pf.load_mesh(*mesh_path, abg=euler_angle, usepickle=False)
    print("*** mesh.zlev = ", mesh.zlev)
    return mesh


def load_input(args):
    filename = (
        str(args.FESOM_PATH[0])
        + str(args.FESOM_VARIABLE[0])
        + ".fesom."
        + str(args.FESOM_YEARS[0])
        + ".nc"
    )
    dataset = xr.open_dataset(filename, decode_times=False)
    varname = args.FESOM_VARIABLE[0]
    if varname not in dataset.variables:
        raise KeyError(f"Variable '{varname}' not found in input dataset")
    return dataset, varname


def read_surface_field(dataset, varname, no_timesteps):
    da = dataset[varname]
    timedim = "time" if "time" in da.dims else "T"
    if da.ndim == 3:
        da = da.squeeze()
    if da.ndim != 2:
        raise ValueError(f"Unsupported surface variable shape for {varname}: {da.shape}")

    spatial_dims = [d for d in da.dims if d != timedim]
    if len(spatial_dims) != 1:
        raise ValueError(f"Expected one spatial dimension for {varname}, got dims={da.dims}")
    da = da.transpose(timedim, spatial_dims[0])
    values = np.asarray(da.values)
    if values.shape[0] != no_timesteps:
        raise ValueError(
            f"Unexpected time extent in {varname}: shape={values.shape}, time length={no_timesteps}"
        )
    values = np.where(np.isnan(values), NAN_REPLACE, values)
    return values


def read_3d_field(args, mesh, no_timesteps):
    # Levels filled here follow ``mesh.zlev[:-1]`` (see enumerate below); vertical
    # metadata in write_3d_output must use the same count (not len(zlev)).
    depth_levels = mesh.zlev[:-1]
    no_zlevels = len(depth_levels)
    output = np.zeros((no_timesteps, no_zlevels, mesh.n2d), dtype=np.float64)
    years = (
        int(args.FESOM_YEARS[0])
        if args.FESOM_YEARS[0] == args.FESOM_YEARS[1]
        else [int(args.FESOM_YEARS[0]), int(args.FESOM_YEARS[1])]
    )

    for itime in np.arange(0, no_timesteps, 1, dtype=np.int32):
        print("*   TIME(" + str(itime) + ")")
        for ilevel, depth in enumerate(depth_levels):
            idepth = int(depth)
            print("*   depth=" + str(idepth) + " (" + str(depth) + ") ++ time(" + str(itime) + ")")
            if isinstance(years, list):
                level_data = pf.get_data(
                    args.FESOM_PATH[0],
                    args.FESOM_VARIABLE[0],
                    years,
                    mesh,
                    depth=idepth,
                    how="mean",
                    use_cftime=True,
                )
            else:
                level_data = pf.get_data(
                    args.FESOM_PATH[0],
                    args.FESOM_VARIABLE[0],
                    years,
                    mesh,
                    depth=idepth,
                    how="mean",
                )
            output[itime, ilevel, :] = np.where(np.isnan(level_data), NAN_REPLACE, level_data)

    return output


def write_common_metadata(output, dataset):
    output.type = "FESOM_ocean_data"
    output.title = "FESOM hydrographic ocean data"
    output.model_domain = "Generic"
    output.insitution = "Alfred Wegener Institute for Polar and Marine Research (AWI)"
    output.address = "Handelshafen, Bremerhaven, Germany"
    output.department = "Climate"
    output.contact = "Christian Rodehacke"
    output.web = "http://www.awi.de"
    output.acknowledgment = "ZUWEISS, BMWF, Germany"
    output.history = "Created " + time.ctime(time.time())

    time_var = output.createVariable("time", "f4", ("time",))
    time_var.long_name = "time"
    time_var.axis = "T"
    time_var.units = ""
    if hasattr(dataset.variables["time"], "calendar"):
        time_var.calendar = dataset.variables["time"].calendar

    horizontal_var = output.createVariable("horizontal", "i4", ("horizontal",))
    horizontal_var.long_name = "horizontal_grid_id"
    horizontal_var.description = "Number of horizontal grid element"
    horizontal_var.axis = "X"

    lon_var = output.createVariable("longitude", "f4", ("horizontal",))
    lat_var = output.createVariable("latitude", "f4", ("horizontal",))
    lon_var.long_name = "longitude"
    lon_var.units = "degrees east"
    lat_var.long_name = "latitude"
    lat_var.units = "degrees north"

    return time_var, horizontal_var, lon_var, lat_var


def write_surface_output(args, dataset, mesh, values):
    output = Dataset(args.FESOM_OUTPUT[0], "w", format="NETCDF3_64BIT_OFFSET")
    output.createDimension("time", None)
    output.createDimension("horizontal", mesh.n2d)

    time_var, horizontal_var, lon_var, lat_var = write_common_metadata(output, dataset)
    field_var = output.createVariable(
        args.FESOM_VARIABLE[0],
        "f4",
        ("time", "horizontal"),
        fill_value=NAN_REPLACE,
    )
    field_var.long_name = getattr(dataset.get(args.FESOM_VARIABLE[0]), "long_name", args.FESOM_VARIABLE[0])
    field_var.units = getattr(dataset.get(args.FESOM_VARIABLE[0]), "units", "")
    field_var.coordinates = "longitude latitude"
    field_var.description = ""

    timedim = "time" if "time" in dataset.dims else "T"
    horizontal_var[:] = np.arange(0, mesh.n2d, dtype=np.int32)
    lon_var[:] = mesh.x2
    lat_var[:] = mesh.y2
    time_var[:] = np.asarray(dataset[timedim].values)[: int(dataset.sizes[timedim])]
    field_var[:] = values
    output.close()


def write_3d_output(args, dataset, mesh, values):
    n_z = int(values.shape[1])
    z_coords = np.asarray(mesh.zlev[:-1], dtype=np.float64)
    if z_coords.shape[0] != n_z:
        raise ValueError(
            f"3D vertical mismatch: data has level count {n_z}, mesh.zlev[:-1] has {z_coords.shape[0]}"
        )

    output = Dataset(args.FESOM_OUTPUT[0], "w", format="NETCDF3_64BIT_OFFSET")
    output.createDimension("time", None)
    output.createDimension("level", n_z)
    output.createDimension("horizontal", mesh.n2d)

    time_var, horizontal_var, lon_var, lat_var = write_common_metadata(output, dataset)
    level_var = output.createVariable("level", "i4", ("level",))
    depth_var = output.createVariable("depth", "f4", ("level",))
    field_var = output.createVariable(
        args.FESOM_VARIABLE[0],
        "f4",
        ("time", "level", "horizontal"),
        fill_value=NAN_REPLACE,
    )

    level_var.long_name = "model_level"
    level_var.description = "Level number"
    level_var.axis = "Z"
    depth_var.long_name = "depth"
    depth_var.units = "m"
    depth_var.description = "depth below sea level"
    depth_var.positive = "down"

    field_var.long_name = getattr(dataset.get(args.FESOM_VARIABLE[0]), "long_name", args.FESOM_VARIABLE[0])
    field_var.units = getattr(dataset.get(args.FESOM_VARIABLE[0]), "units", "")
    field_var.coordinates = "longitude latitude"
    field_var.description = ""

    timedim = "time" if "time" in dataset.dims else "T"
    level_var[:] = np.arange(0, n_z, dtype=np.int32)
    depth_var[:] = z_coords
    horizontal_var[:] = np.arange(0, mesh.n2d, dtype=np.int32)
    lon_var[:] = mesh.x2
    lat_var[:] = mesh.y2
    time_var[:] = np.asarray(dataset[timedim].values)[: int(dataset.sizes[timedim])]
    field_var[:] = values
    output.close()


def main():
    args = parse_arguments()
    print("* Start at " + time.ctime(time.time()))

    mesh = load_mesh(args.FESOM_MESH, bool(args.FESOM_MESH_ROTATED))
    dataset, varname = load_input(args)
    timedim = "time" if "time" in dataset.dims else "T"
    no_timesteps = int(dataset.sizes[timedim])
    print("no_timesteps: ", no_timesteps)

    is_surface_layer = detect_surface_layer(dataset, varname, mesh)
    if is_surface_layer:
        print("*   single-layer surface/interface field detected")
        values = read_surface_field(dataset, varname, no_timesteps)
        write_surface_output(args, dataset, mesh, values)
    else:
        print("*   multi-level field detected")
        values = read_3d_field(args, mesh, no_timesteps)
        write_3d_output(args, dataset, mesh, values)

    dataset.close()
    print("End   at " + time.ctime(time.time()))
    print("    ... Bye")


if __name__ == "__main__":
    main()
