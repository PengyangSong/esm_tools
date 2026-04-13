#!/usr/bin/env python3

import os
from netCDF4 import Dataset
import numpy as np

def split_FESOM2_restart(input_file, output_folder, shared_vars):
    os.makedirs(output_folder, exist_ok=True)
    dsin = Dataset(input_file, 'r')

    # Extract shared variable data and attributes
    shared_data = {}
    for name in shared_vars:
        if name in dsin.variables:
            var = dsin.variables[name]
            shared_data[name] = {
                'data': var[:],
                'datatype': var.datatype,
                'dimensions': var.dimensions,
                'attributes': {attr: var.getncattr(attr) for attr in var.ncattrs()}
            }

    # Process all variables except shared ones
    for var_name in dsin.variables:
        if var_name in shared_vars:
            continue

        var = dsin.variables[var_name]
        data = var[:]
        dims = var.dimensions

        # Reorder dimensions if needed
        if dims in [('time', 'node', 'nz'), ('time', 'node', 'nz_1'),
                    ('time', 'elem', 'nz'), ('time', 'elem', 'nz_1')]:
            data = data.transpose(0, 2, 1)
            dims = (dims[0], dims[2], dims[1])

        output_path = os.path.join(output_folder, f'{var_name}.nc')
        dsout = Dataset(output_path, 'w', format='NETCDF4')

        # Create dimensions
        for dim_name in dims:
            if dim_name not in dsout.dimensions:
                dim = dsin.dimensions[dim_name]
                dsout.createDimension(dim_name, len(dim) if not dim.isunlimited() else None)

        # Create shared dimensions
        for name in shared_vars:
            for dim_name in shared_data[name]['dimensions']:
                if dim_name not in dsout.dimensions:
                    dim = dsin.dimensions[dim_name]
                    dsout.createDimension(dim_name, len(dim) if not dim.isunlimited() else None)

        # Copy shared variables
        for name in shared_vars:
            meta = shared_data[name]
            meta_var = dsout.createVariable(name, meta['datatype'], meta['dimensions'])
            meta_var[:] = meta['data']
            meta_var.setncatts(meta['attributes'])

        # Create main variable
        print(f"{var_name}: {dims} → shape {np.shape(data)}")
        out_var = dsout.createVariable(var_name, var.datatype, dims)
        out_var[:] = data
        out_var.setncatts({attr: var.getncattr(attr) for attr in var.ncattrs()})

        # Copy global attributes
        dsout.setncatts({attr: dsin.getncattr(attr) for attr in dsin.ncattrs()})
        dsout.close()

    dsin.close()

def copy_variable(source_var, target_var, folder):
    source_path = os.path.join(folder, f'{source_var}.nc')
    target_path = os.path.join(folder, f'{target_var}.nc')

    if os.path.exists(target_path):
        print(f"Skipped: {target_var}.nc already exists.")
        return

    with Dataset(source_path, 'r') as dssrc:
        with Dataset(target_path, 'w', format='NETCDF4') as dstgt:
            for dim_name in dssrc.dimensions:
                dim = dssrc.dimensions[dim_name]
                dstgt.createDimension(dim_name, len(dim) if not dim.isunlimited() else None)

            for var_name in dssrc.variables:
                new_name = target_var if var_name == source_var else var_name
                var = dssrc.variables[var_name]
                out_var = dstgt.createVariable(new_name, var.datatype, var.dimensions)
                out_var[:] = var[:]
                out_var.setncatts({attr: var.getncattr(attr) for attr in var.ncattrs()})

            dstgt.setncatts({attr: dssrc.getncattr(attr) for attr in dssrc.ncattrs()})

    print(f"Copied: {source_var}.nc → {target_var}.nc")

def main():
    restart_dir = os.environ.get('last_exp_dir')
    restart_year = os.environ.get('last_restart_year')

    if not restart_dir or not restart_year:
        raise ValueError("Environment variables last_exp_dir and last_restart_year must be set.")

    shared_vars = ['time', 'iter']

    split_FESOM2_restart(
        input_file=os.path.join(restart_dir, f'fesom.{restart_year}.oce.restart.nc'),
        output_folder=os.path.join(restart_dir, f'fesom.{restart_year}.oce.restart'),
        shared_vars=shared_vars
    )

    split_FESOM2_restart(
        input_file=os.path.join(restart_dir, f'fesom.{restart_year}.ice.restart.nc'),
        output_folder=os.path.join(restart_dir, f'fesom.{restart_year}.ice.restart'),
        shared_vars=shared_vars
    )

    oce_folder = os.path.join(restart_dir, f'fesom.{restart_year}.oce.restart')
    copy_variable('salt', 'salt_M1', oce_folder)
    copy_variable('temp', 'temp_M1', oce_folder)

if __name__ == '__main__':
    main()

