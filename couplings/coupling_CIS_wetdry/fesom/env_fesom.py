"""
FESOM env for awiesm-pism-vilma CIS when script_dir is couplings/coupling_CIS/fesom.

Wetdry: ICE_TO_FESOM forced on, FESOM_PREP_ICEBERG_DISCHARGE, awicm chunk tags,
RESTART_DIR_fesom for ice→ocean restart handling.

Per-chunk wetdry paths use WETDRY_CHUNK_TAG (YYYYMMDD-YYYYMMDD) so max/submesh dirs are
not overwritten between iterative chunks. Override WETDRY_MAXMESH_DIR / WETDRY_SUBMESH_DIR
in the runscript if needed.
"""

import os


def _wetdry_chunk_tag(general):
    """Path-safe tag from chunk start/end (same calendar convention as CHUNK_*_DATE_*)."""
    cs = general["chunk_start_date"]
    ce = general["chunk_end_date"]

    def _ymd(d):
        return f"{int(d.syear):04d}{int(d.smonth):02d}{int(d.sday):02d}"

    return f"{_ymd(cs)}-{_ymd(ce)}"


def prepare_environment(config):
    general = config["general"]
    couple_dir = general["experiment_couple_dir"].rstrip("/")
    chunk_tag = _wetdry_chunk_tag(general)
    configured_mesh_dir = config["fesom"]["mesh_dir"]
    active_submesh_dir = f"{couple_dir}/fesom_wetdry_submesh_active"
    active_submesh_nod2d = f"{active_submesh_dir}/nod2d.out"
    mesh_state_file = f"{couple_dir}/wetdry_last_mesh_dir.txt"
    mesh_nx_state_file = f"{couple_dir}/wetdry_last_mesh_nx.txt"

    def _mesh_real_if_valid(path):
        if not path:
            return ""
        real = os.path.realpath(path)
        return real if os.path.isfile(os.path.join(real, "nod2d.out")) else ""

    configured_mesh_dir_real = _mesh_real_if_valid(configured_mesh_dir)
    active_mesh_dir_real = _mesh_real_if_valid(active_submesh_dir)

    state_mesh_dir_real = ""
    if os.path.isfile(mesh_state_file):
        try:
            with open(mesh_state_file, "r", encoding="utf-8") as f:
                state_mesh_dir_real = _mesh_real_if_valid(f.read().strip())
        except OSError:
            state_mesh_dir_real = ""

    if not state_mesh_dir_real and configured_mesh_dir_real:
        # Bootstrap state from configured mesh (e.g. chunk 1 base mesh).
        state_mesh_dir_real = configured_mesh_dir_real
        try:
            with open(mesh_state_file, "w", encoding="utf-8") as f:
                f.write(state_mesh_dir_real + "\n")
        except OSError:
            pass

    state_mesh_nx = None
    if os.path.isfile(mesh_nx_state_file):
        try:
            with open(mesh_nx_state_file, "r", encoding="utf-8") as f:
                nx_txt = f.read().strip()
                nx_val = int(nx_txt)
                if nx_val > 0:
                    state_mesh_nx = nx_val
        except (OSError, ValueError):
            state_mesh_nx = None

    # Bootstrap wetdry nx state from configured mesh nod2d if missing.
    if state_mesh_nx is None and configured_mesh_dir_real:
        nod2d_path = os.path.join(configured_mesh_dir_real, "nod2d.out")
        try:
            with open(nod2d_path, "r", encoding="utf-8") as f:
                nx_val = int(f.readline().strip().split()[0])
                if nx_val > 0:
                    state_mesh_nx = nx_val
                    with open(mesh_nx_state_file, "w", encoding="utf-8") as outf:
                        outf.write(f"{nx_val}\n")
        except (OSError, ValueError, IndexError):
            state_mesh_nx = None

    if state_mesh_nx is not None:
        # Ensure chunk prepare/namcouple generation uses the previous chunk's submesh node count.
        config["fesom"]["nx"] = state_mesh_nx
        # OASIS namcouple uses component grid sizes (fesom.grids.feom.nx),
        # so update that explicitly as well.
        if "grids" in config["fesom"] and "feom" in config["fesom"]["grids"]:
            config["fesom"]["grids"]["feom"]["nx"] = state_mesh_nx

    mesh_dir_for_chunk = (
        active_mesh_dir_real or configured_mesh_dir_real or os.path.realpath(configured_mesh_dir)
    )
    # Optional explicit remap source from runscript. If unset, use previous-state mesh.
    configured_last_mesh = config["fesom"].get("wetdry_last_mesh_dir", "")
    if configured_last_mesh:
        last_mesh_for_remap = _mesh_real_if_valid(configured_last_mesh) or os.path.realpath(
            configured_last_mesh
        )
    elif state_mesh_dir_real:
        # Chunk handoff state (absolute physical path) written by coupling.
        last_mesh_for_remap = state_mesh_dir_real
    else:
        last_mesh_for_remap = mesh_dir_for_chunk

    environment_dict = {
        "ICE_TO_FESOM": 1,
        "FESOM_PREP_ICEBERG_DISCHARGE": int(
            config["fesom"].get("use_icebergs", False).__bool__()
        ),
        "FESOM_TO_ICE": int(general["first_run_in_chunk"]),
        # Coupling scripts should work on physical paths (no symlink ambiguity).
        # FESOM namelist MeshPath remains controlled by runscript YAML (fesom.mesh_dir).
        "MESH_DIR_fesom": mesh_dir_for_chunk,
        "MESH_ROTATED_fesom": config["fesom"]["mesh_rotated"],
        "DATA_DIR_fesom": config["fesom"]["experiment_outdata_dir"],
        "COUPLE_DIR": general["experiment_couple_dir"],
        "number_of_years_for_forcing": config["model1"]["chunk_size"],
        "CHUNK_SIZE_pism_standalone": config["model2"]["chunk_size"],
        "CHUNK_START_DATE_fesom": general["chunk_start_date"],
        "CHUNK_END_DATE_fesom": general["chunk_end_date"],
        "CHUNK_START_DATE_awicm": general["chunk_start_date"],
        "CHUNK_END_DATE_awicm": general["chunk_end_date"],
        "PYFESOM_PATH": "/pf/a/a270124/pyfesom2/",
        "EXP_ID": general["command_line_config"]["expid"],
        "iter_coup_regrid_method_ice2oce": "INTERPOLATE",
        "MACHINE": config["computer"]["name"],
        "ICEBERG_DIR": config["fesom"].get("iceberg_dir", ""),
        "RESTART_DIR_fesom": config["fesom"].get("experiment_restart_in_dir", ""),
        "WETDRY_CHUNK_TAG": chunk_tag,
        "WETDRY_MAXMESH_DIR": f"{couple_dir}/fesom_wetdry_maxmesh_{chunk_tag}",
        "WETDRY_SUBMESH_DIR": f"{couple_dir}/fesom_wetdry_submesh_{chunk_tag}",
        "WETDRY_SUBMESH_ACTIVE_DIR": active_submesh_dir,
        "WETDRY_LAST_MESH_DIR": last_mesh_for_remap,
        "WETDRY_LAST_MESH_STATE_FILE": mesh_state_file,
        "WETDRY_LAST_MESH_NX_STATE_FILE": mesh_nx_state_file,
        # submesh_partition: fesom_ini scratch work directory (per chunk)
        "WETDRY_SUBMESH_PARTITION_WORK_DIR": f"{couple_dir}/fesom_wetdry_partition_work_{chunk_tag}",
        "WETDRY_RESTART_REMAP_WORK_DIR": f"{couple_dir}/fesom_wetdry_restart_remap_{chunk_tag}",
        # Default: chunk calendar start year; override if your fesom.<year>.oce.restart uses another year.
        "WETDRY_FESOM_RESTART_YEAR": str(int(general["chunk_start_date"].syear)-1),
    }

    environment_dict["FUNCTION_PATH"] = os.path.dirname(os.path.abspath(__file__))
    print(environment_dict)
    return environment_dict
