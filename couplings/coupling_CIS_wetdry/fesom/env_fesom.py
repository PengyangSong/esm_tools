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

    environment_dict = {
        "ICE_TO_FESOM": 1,
        "FESOM_PREP_ICEBERG_DISCHARGE": int(
            config["fesom"].get("use_icebergs", False).__bool__()
        ),
        "FESOM_TO_ICE": int(general["first_run_in_chunk"]),
        "MESH_DIR_fesom": config["fesom"]["mesh_dir"],
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
        # submesh_partition: fesom_ini scratch work directory (per chunk)
        "WETDRY_SUBMESH_PARTITION_WORK_DIR": f"{couple_dir}/fesom_wetdry_partition_work_{chunk_tag}",
        "WETDRY_RESTART_REMAP_WORK_DIR": f"{couple_dir}/fesom_wetdry_restart_remap_{chunk_tag}",
        # Default: chunk calendar start year; override if your fesom.<year>.oce.restart uses another year.
        "WETDRY_FESOM_RESTART_YEAR": str(int(general["chunk_start_date"].syear)-1),
    }

    environment_dict["FUNCTION_PATH"] = os.path.dirname(os.path.abspath(__file__))
    print(environment_dict)
    return environment_dict
