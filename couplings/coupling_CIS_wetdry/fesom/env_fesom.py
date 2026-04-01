"""
FESOM env for awiesm-pism-vilma CIS when script_dir is couplings/coupling_CIS/fesom.

Wetdry: ICE_TO_FESOM forced on, FESOM_PREP_ICEBERG_DISCHARGE, awicm chunk tags,
RESTART_DIR_fesom for ice→ocean restart handling.
"""

import os


def prepare_environment(config):
    environment_dict = {
        "ICE_TO_FESOM": 1,
        "FESOM_PREP_ICEBERG_DISCHARGE": int(config["fesom"].get("use_icebergs", False)),
        "FESOM_TO_ICE": int(config["general"]["first_run_in_chunk"]),
        "MESH_DIR_fesom": config["fesom"]["mesh_dir"],
        "MESH_ROTATED_fesom": config["fesom"]["mesh_rotated"],
        "DATA_DIR_fesom": config["fesom"]["experiment_outdata_dir"],
        "COUPLE_DIR": config["general"]["experiment_couple_dir"],
        "number_of_years_for_forcing": config["model1"]["chunk_size"],
        "CHUNK_SIZE_pism_standalone": config["model2"]["chunk_size"],
        "CHUNK_START_DATE_fesom": config["general"]["chunk_start_date"],
        "CHUNK_END_DATE_fesom": config["general"]["chunk_end_date"],
        "CHUNK_START_DATE_awicm": config["general"]["chunk_start_date"],
        "CHUNK_END_DATE_awicm": config["general"]["chunk_end_date"],
        "PYFESOM_PATH": "/pf/a/a270124/pyfesom2/",
        "EXP_ID": config["general"]["command_line_config"]["expid"],
        "iter_coup_regrid_method_ice2oce": "INTERPOLATE",
        "MACHINE": config["computer"]["name"],
        "ICEBERG_DIR": config["fesom"].get("iceberg_dir", ""),
        "RESTART_DIR_fesom": config["fesom"].get("experiment_restart_in_dir", ""),
    }

    environment_dict["FUNCTION_PATH"] = os.path.dirname(os.path.abspath(__file__))
    print(environment_dict)
    return environment_dict
