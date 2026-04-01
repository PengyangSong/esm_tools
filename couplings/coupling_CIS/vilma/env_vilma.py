"""
VILMA env for awiesm-pism-vilma CIS when script_dir is couplings/coupling_CIS/vilma.
"""

import os


def prepare_environment(config):
    vilma = config.get("vilma", {})
    general = config["general"]
    environment_dict = {
        "ICE_TO_VILMA": 1,
        "VILMA_TO_ICE": 1,
        "VILMA_GRID_input": vilma.get("grid_input", "n128"),
        "COUPLE_DIR": general["experiment_couple_dir"],
        "solidearth_ice_thickness_file": general["experiment_couple_dir"]
        + "/ice_thickness.nc",
        "ADD_UNCHANGED_ICE": vilma.get("add_unchanged_ice", False),
        "EISLASTFILE_vilma": vilma.get("experiment_input_dir", "")
        + "/"
        + vilma.get("eislastfile", "eislastfile.nc"),
        "RUN_NUMBER_vilma": general["run_number"],
        "RUN_DATE_STAMP": general["run_datestamp"],
        "LAST_RUN_DATE_STAMP": general["last_run_datestamp"],
        "RESTART_DIR_vilma": vilma.get("experiment_restart_out_dir", ""),
        "OUTDATA_DIR_vilma": vilma.get("experiment_outdata_dir", ""),
        "INITIAL_YEAR_vilma": general["initial_date"].syear,
        "NYEAR_vilma_standalone": general.get("nyear", 1000),
        "FINAL_YEAR_vilma": general["final_date"].syear,
        "EISLASTCONF_vilma": vilma.get("experiment_config_dir", "")
        + "/inp/"
        + vilma.get("eislastconf", "loadh.inp"),
    }
    if environment_dict["ADD_UNCHANGED_ICE"] is False:
        environment_dict["ADD_UNCHANGED_ICE"] = 0
    elif environment_dict["ADD_UNCHANGED_ICE"] is True:
        environment_dict["ADD_UNCHANGED_ICE"] = 1

    environment_dict["FUNCTION_PATH"] = os.path.dirname(os.path.abspath(__file__))
    return environment_dict
