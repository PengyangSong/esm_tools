import os
import glob
import numpy as np
from netCDF4 import Dataset


class FESOM2_mesh:
    def __init__(self, dim_n, dim_e, dim_z, zbar, zmid,
                 lon, lat, elem, nlvl, elvl,
                 bot_raw, srf_raw,
                 map_n, map_e,
                 nlvl_cavity, elvl_cavity,
                 is_submesh=False, with_cavity=False,
                 use_partial_cell=True):
        self.dim_n = dim_n
        self.dim_e = dim_e
        self.dim_z = dim_z
        self.zbar = zbar
        self.zmid = zmid
        self.lon = lon
        self.lat = lat
        self.elem = elem
        self.nlvl = nlvl
        self.elvl = elvl
        self.bot_raw = bot_raw
        self.srf_raw = srf_raw
        self.map_n = map_n
        self.map_e = map_e
        self.nlvl_cavity = nlvl_cavity
        self.elvl_cavity = elvl_cavity
        self.is_submesh = is_submesh
        self.with_cavity = with_cavity
        self.use_partial_cell = use_partial_cell
        
        # Fix for python because index starts from zero
        self.nl_min = 0
        self.nl_max = self.dim_z -1

    def __repr__(self):
        return f"<FESOM2_mesh: {self.dim_n} nodes, {self.dim_e} elems, {self.dim_z} layers>"


def ensure_trailing_slash(path_in: str) -> str:
    if not path_in.endswith(os.sep):
        return path_in + os.sep
    return path_in


def fesom2_read_mesh(mesh_path: str) -> FESOM2_mesh:
    mesh_path = ensure_trailing_slash(mesh_path)

    # --- Read node information ---
    node_file = os.path.join(mesh_path, "nod2d.out")
    if not os.path.isfile(node_file):
        raise FileNotFoundError("Missing nod2d.out file")
    else:
        with open(node_file) as f:
            num_nodes = int(f.readline().strip())
            node_data = np.loadtxt(f, max_rows=num_nodes)

    # --- Read element information ---
    elem_file = os.path.join(mesh_path, "elem2d.out")
    if not os.path.isfile(elem_file):
        raise FileNotFoundError("Missing elem2d.out file")
    else:
        with open(elem_file) as f:
            num_elems = int(f.readline().strip())
            elem_data = np.loadtxt(f, dtype=int, max_rows=num_elems)

    # --- Read vertical layer information ---
    aux_file = os.path.join(mesh_path, "aux3d.out")
    if not os.path.isfile(aux_file):
        raise FileNotFoundError("Missing aux3d.out file")
    else:
        with open(aux_file) as f:
            num_layers = int(f.readline().strip())
            zbar = np.loadtxt(f, max_rows=num_layers)
            bot_raw = np.loadtxt(f, max_rows=num_nodes)
            zmid = 0.5 * (zbar[:-1] + zbar[1:])

    # --- Read cavity information ---
    cav_file = os.path.join(mesh_path, "cavity_depth@node.out")
    if os.path.isfile(cav_file):
        srf_raw = np.loadtxt(cav_file)
    else:
        srf_raw = np.zeros_like(bot_raw)

    # --- Read number of layers per node and element ---
    nlvl = np.loadtxt(os.path.join(mesh_path, "nlvls.out")).astype(int)
    elvl = np.loadtxt(os.path.join(mesh_path, "elvls.out")).astype(int)

    # --- Cavity layer information ---
    cavity_nlvl_file = os.path.join(mesh_path, "cavity_nlvls.out")
    cavity_elvl_file = os.path.join(mesh_path, "cavity_elvls.out")
    if os.path.isfile(cavity_nlvl_file) and os.path.isfile(cavity_elvl_file):
        with_cavity = True
        nlvl_cavity = np.loadtxt(cavity_nlvl_file).astype(int)
        elvl_cavity = np.loadtxt(cavity_elvl_file).astype(int)
    else:
        with_cavity = False
        nlvl_cavity = np.ones(num_nodes, dtype=int)
        elvl_cavity = np.ones(num_elems, dtype=int)

    # --- Submesh mapping ---
    map_n_file = os.path.join(mesh_path, "map_nod.out")
    map_e_file = os.path.join(mesh_path, "map_elem.out")
    if os.path.isfile(map_n_file) and os.path.isfile(map_e_file):
        is_submesh = True
        map_n = np.loadtxt(map_n_file).astype(int)
        map_e = np.loadtxt(map_e_file).astype(int)
    else:
        is_submesh = False
        map_n = np.arange(1, num_nodes+1).astype(int)
        map_e = np.arange(1, num_elems+1).astype(int)

    # Modify bottom depth according to FESOM2 src code
    thers_zbar_lev = 5
    thers_depth = zbar[thers_zbar_lev-1]  # Python index starts from zero
    bot_raw = np.where(bot_raw > thers_depth, thers_depth, bot_raw)

    # Fix for python because index starts from zero
    elem_data = elem_data - 1
    nlvl = nlvl - 1
    elvl = elvl - 1
    nlvl_cavity = nlvl_cavity - 1
    elvl_cavity = elvl_cavity - 1
    map_n = map_n - 1
    map_e = map_e - 1

    # --- Construct FESOM2_mesh object ---
    mesh = FESOM2_mesh(
        dim_n=num_nodes,
        dim_e=num_elems,
        dim_z=num_layers-1,
        zbar=zbar,
        zmid=zmid,
        lon=node_data[:, 1],
        lat=node_data[:, 2],
        elem=elem_data,
        nlvl=nlvl,
        elvl=elvl,
        bot_raw=bot_raw,
        srf_raw=srf_raw,
        map_n=map_n,
        map_e=map_e,
        nlvl_cavity=nlvl_cavity,
        elvl_cavity=elvl_cavity,
        is_submesh=is_submesh,
        with_cavity=with_cavity,
    )

    return mesh


def fesom2_init_thickness_ale(mesh: "FESOM2_mesh") -> "FESOM2_mesh":
    """
    From FESOM2: Calculate initial thickness at elements and nodes.
    Updates the FESOM2_mesh object with thickness-related fields.
    """

    zbar = mesh.zbar
    zmid = mesh.zmid

    # --- basic calculation n2e ---
    ebot_raw = np.mean(mesh.bot_raw[mesh.elem], axis=1)
    esrf_raw = np.mean(mesh.srf_raw[mesh.elem], axis=1)

    # --- initialise arrays ---
    ebot = np.zeros(mesh.dim_e)
    esrf = np.zeros(mesh.dim_e)
    bottom_elem_thickness = np.zeros(mesh.dim_e)
    elem_dz = np.zeros((mesh.dim_e, mesh.dim_z))

    nbot = np.zeros(mesh.dim_n)
    nsrf = np.zeros(mesh.dim_n)
    bottom_node_thickness = np.zeros(mesh.dim_n)
    node_dz = np.zeros((mesh.dim_n, mesh.dim_z))

    # --- init_bottom_elem_thickness ---
    if mesh.use_partial_cell:
        for ij in range(mesh.dim_e):
            nle = int(mesh.elvl[ij])
            dd = ebot_raw[ij]
            if dd < zbar[nle]:
                if nle == mesh.nl_max + 1:
                    ebot[ij] = max(dd, zbar[nle] + (zbar[nle] - zmid[nle-1]))
                else:
                    ebot[ij] = max(zmid[nle], dd)
                bottom_elem_thickness[ij] = zbar[nle-1] - ebot[ij]
            else:
                ebot[ij] = min(zmid[nle-1], dd)
                bottom_elem_thickness[ij] = zbar[nle-1] - ebot[ij]
    else:
        for ij in range(mesh.dim_e):
            nle = int(mesh.elvl[ij])
            ebot[ij] = zbar[nle]
            bottom_elem_thickness[ij] = zbar[nle-1] - zbar[nle]

    # --- init_surface_elem_depth ---
    for ij in range(mesh.dim_e):
        ule = int(mesh.elvl_cavity[ij])
        if ule == mesh.nl_min:
            continue
        if mesh.use_partial_cell:
            dd = esrf_raw[ij]
            if dd < zbar[ule]:
                esrf[ij] = max(zmid[ule], dd)
            else:
                esrf[ij] = min(zmid[ule-1], dd)
        else:
            esrf[ij] = zbar[ule]

    # --- define nod_in_elem2D, nlevels_nod2D_min ---
    nod_in_elem2D = np.full((mesh.dim_n, 10), np.nan)
    nod_in_elem2D_num = np.zeros(mesh.dim_n, dtype=int)

    for ij in range(mesh.dim_e):
        for k in range(3):
            n = mesh.elem[ij, k]
            nod_in_elem2D_num[n] += 1
            nod_in_elem2D[n, nod_in_elem2D_num[n] - 1] = ij # Python index

    nlevels_nod2D_min = np.zeros(mesh.dim_n, dtype=int)
    for ij in range(mesh.dim_n):
        elems = nod_in_elem2D[ij, :nod_in_elem2D_num[ij]].astype(int)
        nlevels_nod2D_min[ij] = np.min(mesh.elvl[elems])

    # --- init_bottom_node_thickness ---
    if mesh.use_partial_cell:
        for ij in range(mesh.dim_n):
            nln = int(mesh.nlvl[ij])
            elems = nod_in_elem2D[ij, :nod_in_elem2D_num[ij]].astype(int)
            nbot[ij] = np.min(ebot[elems])
            bottom_node_thickness[ij] = zbar[nln-1] - nbot[ij]
    else:
        for ij in range(mesh.dim_n):
            nln = int(mesh.nlvl[ij])
            nbot[ij] = zbar[nln]
            bottom_node_thickness[ij] = zbar[nln-1] - nbot[ij]

    # --- init_surface_node_depth ---
    for ij in range(mesh.dim_n):
        uln = int(mesh.nlvl_cavity[ij])
        if uln == mesh.nl_min:
            continue
        if mesh.use_partial_cell:
            elems = nod_in_elem2D[ij, :nod_in_elem2D_num[ij]].astype(int)
            nsrf[ij] = np.max(esrf[elems])
        else:
            nsrf[ij] = zbar[uln]

    # --- init_thickness_ale elem and node ---
    for ij in range(mesh.dim_e):
        k1 = int(mesh.elvl_cavity[ij])
        k2 = int(mesh.elvl[ij]) - 1
        elem_dz[ij, k1] = esrf[ij] - zbar[k1+1]
        for k in range(k1+1, k2):
            elem_dz[ij, k] = zbar[k] - zbar[k+1]
        elem_dz[ij, k2] = bottom_elem_thickness[ij]

    for ij in range(mesh.dim_n):
        k1 = int(mesh.nlvl_cavity[ij])
        k2 = int(mesh.nlvl[ij]) - 1
        node_dz[ij, k1] = nsrf[ij] - zbar[k1+1]
        for k in range(k1+1, k2):
            node_dz[ij, k] = zbar[k] - zbar[k+1]
        node_dz[ij, k2] = bottom_node_thickness[ij]

    # --- attach results to mesh ---
    mesh.ebot = ebot
    mesh.esrf = esrf
    mesh.elem_dz = elem_dz

    mesh.nbot = nbot
    mesh.nsrf = nsrf
    mesh.node_dz = node_dz

    mesh.nod_in_elem2D = nod_in_elem2D
    mesh.nod_in_elem2D_num = nod_in_elem2D_num
    mesh.nlevels_nod2D_min = nlevels_nod2D_min

    return mesh


def fesom2_hbar2hnode(hbar: np.ndarray, mesh: "FESOM2_mesh", ale_type: str) -> np.ndarray:
    """
    Treatment for hnode, either initialise or restart.
    Cannot simply fill with zeros.

    Parameters
    ----------
    hbar : np.ndarray
        1D array of size (mesh.dim_n,) with hbar values.
    mesh : FESOM2_mesh
        FESOM2_mesh object containing node_dz, nlvl_cavity, nlevels_nod2D_min, zbar, dim_n.
    ale_type : str
        ALE type: 'linfs', 'zstar', or 'zlevel'.

    Returns
    -------
    hnode : np.ndarray
        Node thickness array (dim_n, dim_z).
    """
    
    # --- Check ale ---
    if ale_type == "linfs":
        hbar = np.zeros_like(hbar)

    # --- Make sure cavity is always linfs ---
    for ij in range(mesh.dim_n):
        if mesh.nlvl_cavity[ij] > mesh.nl_min:
            hbar[ij] = 0.0

    # --- Default: linfs ---
    hnode = mesh.node_dz.copy()

    # --- hbar to hnode ---
    if ale_type == "zstar":
        for ij in range(mesh.dim_n):
            k1 = int(mesh.nlvl_cavity[ij])
            k2 = int(mesh.nlevels_nod2D_min[ij]) - 1
            zstar_coeff = hbar[ij] / (mesh.zbar[k1] - mesh.zbar[k2]) + 1.0
            hnode[ij, k1:k2] *= zstar_coeff

    elif ale_type == "zlevel":
        hnode[:, mesh.nl_min] = hnode[:, mesh.nl_min] + hbar

    hnode = hnode.T # Python netCDF form

    return hnode


def detect_varshape(var_in: np.ndarray, mesh: "FESOM2_mesh") -> dict:
    """
    Identify the dimensional structure of a FESOM2 variable.
    Returns a dict with keys: is_3d, type ('node'/'elem'), is_transpose (bool).
    """
    sz = var_in.shape
    dim_n, dim_e, dim_z = mesh.dim_n, mesh.dim_e, mesh.dim_z

    info = {"is_3d": False, "type": "", "is_transpose": False}

    # Case 1: node-based 3D
    if sz == (dim_n, dim_z) or sz == (dim_n, dim_z+1):
        info.update({"is_3d": True, "type": "node"})
    # Case 2: element-based 3D
    elif sz == (dim_e, dim_z) or sz == (dim_e, dim_z+1):
        info.update({"is_3d": True, "type": "elem"})
    # Case 3: node-based 3D, transpose
    elif sz == (dim_z, dim_n) or sz == (dim_z+1, dim_n):
        info.update({"is_3d": True, "type": "node", "is_transpose": True})
    # Case 4: element-based 3D, transpose
    elif sz == (dim_z, dim_e) or sz == (dim_z+1, dim_e):
        info.update({"is_3d": True, "type": "elem", "is_transpose": True})
    # Case 5: node-based 2D
    elif sz == (dim_n,):
        info.update({"is_3d": False, "type": "node"})
    # Case 6: element-based 2D
    elif sz == (dim_e,):
        info.update({"is_3d": False, "type": "elem"})
    else:
        raise ValueError(f"Unrecognized variable shape: {sz}")

    return info


def fesom2_map_field(mesh1: "FESOM2_mesh", mesh2: "FESOM2_mesh", var_in: np.ndarray, strategy: str = "zero") -> np.ndarray:
    """
    Unified remapping for 2D/3D FESOM2 fields.
    Supports node/element fields, surface or full-depth, with optional extrapolation.
    """

    # Pre check
    if var_in.size == 1:
        return var_in.copy()
    else:
        #var_in = var_in[0] # skip time dimension
        var_in = np.squeeze(var_in[0])

    info = detect_varshape(var_in, mesh1)

    if info["is_transpose"]:
        var_in = var_in.T

    # Determine mapping scenario
    if info["type"] == "node":
        dim1, dim2 = mesh1.dim_n, mesh2.dim_n
        map1, map2 = mesh1.map_n, mesh2.map_n
        lev_srf1, lev_srf2 = mesh1.nlvl_cavity, mesh2.nlvl_cavity
        lev_bot1, lev_bot2 = mesh1.nlvl, mesh2.nlvl
    elif info["type"] == "elem":
        dim1, dim2 = mesh1.dim_e, mesh2.dim_e
        map1, map2 = mesh1.map_e, mesh2.map_e
        lev_srf1, lev_srf2 = mesh1.elvl_cavity, mesh2.elvl_cavity
        lev_bot1, lev_bot2 = mesh1.elvl, mesh2.elvl

    dim_max = max(map1.max(), map2.max()) + 1 # Python index

    # Prepare output shape
    if info["is_3d"]:
        var_out = np.full((dim2, var_in.shape[1]), np.nan)
    else:
        var_out = np.full((dim2,), np.nan)

    # --- Full-depth field remapping (3D) ---
    if info["is_3d"]:
        for k in range(var_in.shape[1]):
            layer_in = var_in[:, k].copy()

            # Mask cavity and bed in source
            cavity_mask = (k < lev_srf1)
            bed_mask = (k >= lev_bot1)
            layer_in[cavity_mask | bed_mask] = np.nan

            # Interpolate to intermediate mesh
            var_max = np.full((dim_max,), np.nan)
            var_max[map1] = layer_in

            # Map to target mesh
            layer_out = var_max[map2]

            # Mask cavity and bed in target
            cavity_mask_out = (k < lev_srf2)
            bed_mask_out = (k >= lev_bot2)
            layer_out[cavity_mask_out | bed_mask_out] = np.nan

            var_out[:, k] = layer_out

        # Apply extrapolation if needed
        if strategy == "extrap" and info["type"] == "node":
            var_out = fesom2_extrap_nod3D(var_out, mesh2)

        if strategy == "zero":
            var_out[np.isnan(var_out)] = 0.0

    # --- Surface field remapping (2D) ---
    else:
        var_tmp = var_in.copy()

        # Mask cavity in source mesh
        cavity_mask = lev_srf1 > mesh1.nl_min
        var_tmp[cavity_mask] = np.nan
        
        # Interpolate to intermediate mesh
        var_max = np.full(dim_max, np.nan)
        var_max[map1] = var_tmp.flatten()

        # Map to target mesh
        var_out = var_max[map2].reshape(-1, 1)

        # 2D var only filled with zero
        if strategy == "zero":
            var_out[np.isnan(var_out)] = 0.0

    if info["is_transpose"]:
        var_out = var_out.T

    return var_out


def fesom2_extrap_nod3D(arr: np.ndarray, mesh: "FESOM2_mesh") -> np.ndarray:
    """
    Extrapolate missing values in 3D node-based fields (temp/salt).
    Horizontal extrapolation followed by vertical extrapolation.

    Parameters
    ----------
    arr : np.ndarray
        Array of shape (dim_n, dim_z) with NaNs for missing values.
    mesh : FESOM2_mesh
        FESOM2_mesh object containing adjacency and level information.

    Returns
    -------
    arr : np.ndarray
        Extrapolated array (dim_n, dim_z).
    """

    dim_n, dim_z = mesh.dim_n, mesh.dim_z

    # --- Horizontal extrapolation ---
    for k in range(dim_z):
        work_array = arr[:, k].copy()
        success = True

        while success:
            success = False
            tmp_array = work_array.copy()

            # loop over nodal vertices
            for ij in range(dim_n):
                if (np.isnan(work_array[ij]) and
                    k < mesh.nlvl[ij] and
                    k >= mesh.nlvl_cavity[ij]):

                    sum_val = 0.0
                    count = 0

                    # loop over adjacent elements
                    for n in range(mesh.nod_in_elem2D_num[ij]):
                        el = int(mesh.nod_in_elem2D[ij, n])
                        if (k >= mesh.elvl[el] or
                            k < mesh.elvl_cavity[el]):
                            continue

                        # loop over vertices of adjacent element
                        for m in range(3):
                            nd = mesh.elem[el, m]
                            if (not np.isnan(work_array[nd]) and
                                k < mesh.nlvl[nd] and
                                k >= mesh.nlvl_cavity[nd]):
                                sum_val += work_array[nd]
                                count += 1

                    # found valid neighbouring node
                    if count > 0:
                        tmp_array[ij] = sum_val / count
                        # optional debug print:
                        print(f"Extrapolation: nod={ij+1}, lev={k+1}")
                        success = True

            work_array = tmp_array

        arr[:, k] = work_array

    # --- Vertical extrapolation ---
    for ij in range(dim_n):
        # fill cavity (bottom-up)
        for k in range(mesh.nl_max-1, -1, mesh.nl_min-1):
            if not np.isnan(arr[ij, k+1]) and np.isnan(arr[ij, k]):
                arr[ij, k] = arr[ij, k+1]

        # fill bottom (top-down)
        for k in range(mesh.nl_min+1, mesh.nl_max+1):
            if not np.isnan(arr[ij, k-1]) and np.isnan(arr[ij, k]):
                arr[ij, k] = arr[ij, k-1]

    return arr


def fesom2_wetdry_restart():
    print("Modify FESOM2 restart file within a wet-dry mesh. Start!\n...")

    # --- Read environment variables ---
    mesh1_dir = os.getenv("last_mesh_dir")
    mesh2_dir = os.getenv("next_mesh_dir")
    fdir_in   = os.getenv("last_restart_dir")
    fdir_out  = os.getenv("next_restart_dir")

    if not all([mesh1_dir, mesh2_dir, fdir_in, fdir_out]):
        raise RuntimeError("Environment variables mesh1_dir, mesh2_dir, fdir_in, fdir_out must be set.")

    # --- Read meshes as FESOM2_mesh objects ---
    mesh1 = fesom2_read_mesh(mesh1_dir)
    mesh2 = fesom2_read_mesh(mesh2_dir)
    
    mesh1 = fesom2_init_thickness_ale(mesh1)
    mesh2 = fesom2_init_thickness_ale(mesh2)
    
    # --- Settings ---
    dyn_zero = False
    ale_type = "zstar"
    
    # --- Get list of NetCDF files ---
    file_list = glob.glob(os.path.join(fdir_in, "*.nc"))
    file_list.sort()
    
    for fname_in in file_list:
        fname = os.path.basename(fname_in)
        fname_out = os.path.join(fdir_out, fname)
    
        print(f"\nProcessing file: {fname}")
    
        # --- Open input NetCDF ---
        nc_in = Dataset(fname_in, "r")
    
        # --- Create output NetCDF ---
        nc_out = Dataset(fname_out, "w", format="NETCDF4")
        nc_out.setncattr("creation method", "modified FESOM2 restart for wet dry mesh strategy")
    
        # --- Define dimensions ---
        for dimname, dim in nc_in.dimensions.items():
            dimlen_in = len(dim)
            print(f" input: dimname: {dimname}\t dimlen: {dimlen_in}")
    
            if dimlen_in == mesh1.dim_n:
                dimlen_out = mesh2.dim_n
            elif dimlen_in == mesh1.dim_e:
                dimlen_out = mesh2.dim_e
            else:
                dimlen_out = dimlen_in
    
            nc_out.createDimension(dimname, dimlen_out)
            print(f" output: dimname: {dimname}\t dimlen: {dimlen_out}")
    
        # --- Define variables ---
        for varname, var in nc_in.variables.items():
            print(f" input var: {varname}\t xtype: {var.dtype}\t dimids: {var.dimensions}")
            nc_out.createVariable(varname, var.dtype, var.dimensions)
            print(f" output var: {varname}\t xtype: {var.dtype}\t dimids: {var.dimensions}")
    
        # --- Transfer and remap variables ---
        for varname, var in nc_in.variables.items():
            print(f" Reading variable: {varname}")
            var_in = var[:]
    
            # For temp and salt, NaN should be solved by extrapolation
            if varname.startswith("temp") or varname.startswith("salt"):
                strategy = "extrap"
            else:
                strategy = "zero"
    
            # hnode to hbar
            if varname == "hnode":
                hbar_file = os.path.join(fdir_in, "hbar.nc")
                with Dataset(hbar_file, "r") as hbar_nc:
                    var_in = hbar_nc.variables["hbar"][:]
    
            # Remapping
            var_out = fesom2_map_field(mesh1, mesh2, var_in, strategy)
    
            # Apply dynamic zeroing if needed
            if dyn_zero:
                if varname.startswith(("u", "v", "w", "ssh", "hbar", "hnode")):
                    var_out = np.zeros_like(var_out)
    
            # hbar to hnode
            if varname == "hnode":
                var_out = fesom2_hbar2hnode(var_out, mesh2, ale_type)
    
            nc_out.variables[varname][:] = var_out
            print(f" Written variable: {varname}")
    
        # --- Close files ---
        nc_in.close()
        nc_out.close()
    
    print("Modify FESOM2 restart file within a wet-dry mesh. End!\n...")

