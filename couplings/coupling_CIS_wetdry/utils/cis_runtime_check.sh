#!/usr/bin/env bash
# Validate the native runtime used by coupling_CIS_wetdry.
#
# Usage:
#   cis_runtime_check.sh [runtime-root] [qglb-dir] [priority-nodes]
#
# Individual paths may also be supplied through the WETDRY_* environment
# variables.  The script is intentionally independent of esm_tools so it can
# be run once after installing a user's sandbox.
set -uo pipefail

runtime_root="${1:-${WETDRY_RUNTIME_ROOT:-}}"
qglb_dir="${2:-${WETDRY_QGLB_DIR:-}}"
priority_nodes="${3:-${WETDRY_PRIORITY_NODES:-}}"

if [[ -n "${runtime_root}" ]]; then
    runtime_root="${runtime_root%/}"
fi

meshtools_bin="${WETDRY_MESHTOOLS_BIN:-}"
[[ -z "${meshtools_bin}" && -n "${runtime_root}" ]] && meshtools_bin="${runtime_root}/meshtools_triangle/bin"
fesom_ini="${WETDRY_FESOM_INI:-}"
[[ -z "${fesom_ini}" && -n "${runtime_root}" ]] && fesom_ini="${runtime_root}/bin/fesom_ini.x"
rscript="${RSCRIPT:-${WETDRY_RSCRIPT:-Rscript}}"
r_library="${R_LIBS_USER:-${WETDRY_R_LIBRARY:-}}"
[[ -z "${r_library}" && -n "${runtime_root}" ]] && r_library="${runtime_root}/libs"
python_exe="${WETDRY_PYTHON:-python3}"
pyfesom_path="${PYFESOM_PATH:-${WETDRY_PYFESOM_PATH:-}}"

errors=0
warnings=0
pass_msg() { printf 'OK   %s\n' "$*"; }
warn_msg() { printf 'WARN %s\n' "$*"; warnings=$((warnings + 1)); }
fail_msg() { printf 'FAIL %s\n' "$*"; errors=$((errors + 1)); }

check_command() {
    local name="$1"
    local command_name="${2:-$1}"
    if command -v "$command_name" >/dev/null 2>&1; then
        pass_msg "$name: $(command -v "$command_name")"
    else
        fail_msg "$name: command not found ($command_name)"
    fi
}

check_executable() {
    local name="$1" path="$2"
    if [[ -n "$path" && -x "$path" ]]; then
        pass_msg "$name: $path"
        if command -v ldd >/dev/null 2>&1 && file "$path" 2>/dev/null | grep -q ELF; then
            local ldd_env=()
            if [[ -n "${WETDRY_NATIVE_LD_LIBRARY_PATH:-}" ]]; then
                ldd_env=(env "LD_LIBRARY_PATH=${WETDRY_NATIVE_LD_LIBRARY_PATH}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}")
            fi
            if "${ldd_env[@]}" ldd "$path" 2>/dev/null | grep -q 'not found'; then
                fail_msg "$name: unresolved shared library (ldd $path)"
            else
                pass_msg "$name: shared libraries resolved"
            fi
        fi
    elif [[ -n "$path" ]]; then
        fail_msg "$name: executable missing or not executable ($path)"
    else
        fail_msg "$name: no path configured"
    fi
}

printf 'CIS wet/dry runtime preflight\n'
printf '  runtime root: %s\n' "${runtime_root:-<not set>}"
printf '  mesh bin:     %s\n' "${meshtools_bin:-<not set>}"
printf '  FESOM ini:    %s\n' "${fesom_ini:-<not set>}"
printf '  Rscript:      %s\n' "$rscript"
printf '  R library:    %s\n' "${r_library:-<default R library>}"
printf '  Python:       %s\n' "$python_exe"
printf '\nCompiled tools\n'
check_executable relax150614.x "${meshtools_bin:+${meshtools_bin}/relax150614.x}"
check_executable rg.x "${meshtools_bin:+${meshtools_bin}/rg.x}"
check_executable fesom_ini.x "$fesom_ini"

if [[ -n "$meshtools_bin" ]]; then
    export PATH="${meshtools_bin}:${PATH}"
fi
if [[ -n "$r_library" ]]; then
    export R_LIBS_USER="$r_library"
fi
if [[ -n "$pyfesom_path" ]]; then
    export PYTHONPATH="${pyfesom_path}${PYTHONPATH:+:${PYTHONPATH}}"
fi

printf '\nCommand-line tools\n'
check_command CDO cdo
check_command NCO ncks
check_command ncdump ncdump
check_command ncrename ncrename

printf '\nPython and R packages\n'
if "$python_exe" -c 'import netCDF4, numpy, xarray, pyfesom2' >/dev/null 2>&1; then
    pass_msg "Python imports: netCDF4 numpy xarray pyfesom2"
else
    fail_msg "Python imports: netCDF4 numpy xarray pyfesom2 (using $python_exe)"
fi
if [[ -n "$r_library" ]]; then
    r_check="R_LIBS_USER=$(printf '%q' "$r_library") $rscript"
else
    r_check="$rscript"
fi
if env ${r_library:+R_LIBS_USER="$r_library"} "$rscript" -e 'library(spheRlab); library(ncdf4)' >/dev/null 2>&1; then
    pass_msg "R imports: spheRlab ncdf4"
else
    fail_msg "R imports: spheRlab ncdf4 (using $rscript; R library ${r_library:-default})"
fi

printf '\nExperiment inputs\n'
if [[ -n "$qglb_dir" ]]; then
    if [[ -f "$qglb_dir/global/nod2d.out" ]]; then
        pass_msg "QGLB mesh: $qglb_dir/global/nod2d.out"
    else
        fail_msg "QGLB mesh missing: $qglb_dir/global/nod2d.out"
    fi
else
    warn_msg 'QGLB directory not supplied (pass it as argument 2 or WETDRY_QGLB_DIR)'
fi
if [[ -n "$priority_nodes" ]]; then
    if [[ -f "$priority_nodes" ]]; then
        pass_msg "priority nodes: $priority_nodes"
    else
        fail_msg "priority-node file missing: $priority_nodes"
    fi
else
    warn_msg 'priority-node file not supplied (pass it as argument 3 or WETDRY_PRIORITY_NODES)'
fi

printf '\nResult: %d failure(s), %d warning(s)\n' "$errors" "$warnings"
if [[ "$errors" -gt 0 ]]; then
    exit 1
fi
exit 0
