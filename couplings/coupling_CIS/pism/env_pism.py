"""
PISM workflow environment for couplings/coupling_CIS/pism (CIS couple_in / couple_out).
"""

import os


def _cfg_int(val, default=0):
    """Coerce esm-tools config values (incl. *WithProvenance) to int."""
    if val is None:
        return int(default)
    try:
        return int(val)
    except (TypeError, ValueError):
        if hasattr(val, "value"):
            return int(val.value)
        raise


def _cfg_float(val, default=0.0):
    """Coerce esm-tools config values (incl. *WithProvenance) to float."""
    if val is None:
        return float(default)
    try:
        return float(val)
    except (TypeError, ValueError):
        if hasattr(val, "value"):
            return float(val.value)
        raise


def prepare_environment(config):
    setup_name = config["general"]["setup_name"]
    pism = config.get("pism", config.get(setup_name, {}))
    general = config["general"]
    default_input_grid = general["experiment_couple_dir"] + "/ice.griddes"

    environment_dict = {
        "PISM_TO_SOLID_EARTH": _cfg_int(pism.get("coupled_to_solidearth"), 1),
        "SOLID_EARTH_TO_PISM": _cfg_int(pism.get("coupled_to_solidearth"), 1),
        "PISM_TO_OCEAN": 0,
        "OCEAN_TO_PISM": _cfg_int(general.get("first_run_in_chunk"), 0),
        "COUPLE_DIR": general["experiment_couple_dir"],
        "RESTART_DIR_pism": pism.get("experiment_restart_in_dir", general["experiment_couple_dir"]),
        "ice_bedrock_change_file": general["experiment_couple_dir"] + "/bedrock_change.nc",
        "POOL_DIR_pism": pism.get("pool_dir", ""),
        "DOMAIN_pism": pism.get("domain", "nhem"),
        "EXE_pism": pism.get("executable", "pismr"),
        "RES_pism": pism.get("resolution", "10km"),
        "RUN_NUMBER_pism": general["run_number"],
        "pism_solidearth_initialize_method": pism.get("solidearth_initialize_method", "regrid"),
        "pism_solidearth_initialize_dummyrun_file": pism.get("solidearth_initialize_dummyrun_file", ""),
        "INPUT_GRID_pism": pism.get("input_grid", default_input_grid),
        "INPUT_FILE_pism": pism.get("cli_input_file_pism"),
        "YR0_INI_pism": general["initial_date"].syear,
        "CURRENT_YEAR_pism": general["current_date"].syear,
        "END_YEAR_pism": general["end_date"].syear,
        "NYEAR_pism_standalone": general.get("nyear", 1),
        "RUN_DATE_STAMP": general["run_datestamp"],
        "LAST_RUN_DATE_STAMP": general["last_run_datestamp"],
        "EXP_ID": general["command_line_config"]["expid"],
        "OUTPUT_DIR_pism": pism.get("experiment_outdata_dir", general.get("experiment_couple_dir", "")),
        "DATA_DIR_pism": pism.get("experiment_outdata_dir", general.get("experiment_couple_dir", "")),
        "MACHINE": config["computer"]["name"],
        "iterative_coupling_atmosphere_pism_ablation_method": pism.get("ablation_method", "PDD"),
        "DEBM_EXE": pism.get("debm_path", ""),
        "MY_OBLIQUITY": pism.get("debm_obl", "23.441"),
        "DEBM_BETA": pism.get("debm_beta", 999),
        "iterative_coupling_atmosphere_pism_regrid_method": pism.get("regrid_method", "DOWNSCALE"),
        "REDUCE_TEMP": _cfg_int(pism.get("reduce_temp"), 0),
        "REDUCE_TEMP_BY": pism.get("reduce_temp_by", 1),
        "USE_YMONMEAN": pism.get("use_ymonmean", 0),
        "MULTI_YEAR_MEAN_SMB": pism.get("multi_year_mean_smb", 1),
        "TEMP2_BIAS_FILE": pism.get("temp2_bias_file", ""),
        "DOWNSCALING_LAPSE_RATE": pism.get("lapse_rate", -0.005),
        "DOWNSCALE_PRECIP": pism.get("downscale_precip", 1),
    }
    if pism.get("outdata_targets"):
        environment_dict["latest_ex_file_pism"] = pism["outdata_targets"].get("ex_file", "")
    version = pism.get("version", "1.2")
    environment_dict["VERSION_pism"] = version.replace("github", "").replace("index", "").replace("snowflake", "")[:3]
    environment_dict["EX_INT"] = pism.get("ex_interval", "monthly")
    environment_dict["YR0_pism"] = general["start_date"].syear
    environment_dict["M0_pism"] = general["start_date"].smonth
    environment_dict["D0_pism"] = general["start_date"].sday
    environment_dict["END_MONTH_pism"] = general["end_date"].smonth
    environment_dict["END_DAY_pism"] = general["end_date"].sday
    environment_dict["CHUNK_START_DATE_pism"] = general["chunk_start_date"]
    environment_dict["CHUNK_END_DATE_pism"] = general["chunk_end_date"]
    environment_dict["CHUNK_START_YEAR_pism"] = general["chunk_start_date"].syear
    environment_dict["CHUNK_END_YEAR_pism"] = general["chunk_end_date"].syear
    environment_dict["SPINUP_FILE_pism"] = pism.get("spinup_file", "")
    if config["general"]["run_number"] > 1 and pism.get("restart_in_targets"):
        environment_dict["latest_restart_file_pism"] = pism["restart_in_targets"].get("restart", "")
    if "model2" in config:
        environment_dict["CHUNK_SIZE_pism_standalone"] = config["model2"].get("chunk_size", 1)

    environment_dict["PISM_TO_ATMOSPHERE"] = _cfg_int(general.get("last_run_in_chunk"), 0)
    environment_dict["ATMOSPHERE_TO_PISM"] = _cfg_int(general.get("first_run_in_chunk"), 0)
    chunk_sy = str(general["chunk_start_date"].syear)
    nyear_pism = _cfg_int(general.get("nyear"), 1)
    end_y = int(chunk_sy) + nyear_pism - 1
    expid = general["command_line_config"]["expid"]
    exp_input_dir = (pism.get("experiment_input_dir") or "").rstrip("/")
    if exp_input_dir:
        environment_dict["first_year_in_chunk_input"] = (
            f"{exp_input_dir}/{expid}_pismr_input_{chunk_sy}0101-{end_y}1231.nc"
        )
    else:
        environment_dict["first_year_in_chunk_input"] = ""
    restart_out = pism.get("restart_out_targets") or {}
    environment_dict["last_year_in_chunk_restart"] = restart_out.get("restart", "")

    environment_dict["NYEAR"] = _cfg_int(general.get("nyear"), 1)
    environment_dict["MIN_MON_SELECT"] = _cfg_int(pism.get("select_min_glacial_depth"), 1)
    environment_dict["CRITICAL_THK_FOR_MASK_pism"] = _cfg_float(pism.get("thk_threshold"), 5.0)
    environment_dict["fesom_use_iceberg"] = int(
        bool(pism.get("iceberg_coupling", False))
    )

    environment_dict["FUNCTION_PATH"] = os.path.dirname(os.path.abspath(__file__))
    return environment_dict
