#!/usr/bin/env bash
# Build a private native runtime for coupling_CIS_wetdry.
#
# This is intentionally a setup-time script.  It does not compile anything
# during a model coupling event.  Run it from any directory; all installations
# are below a separate runtime root unless overridden.  This avoids cloning
# dependencies into the esm_tools source tree when the script is run from
# coupling_CIS_wetdry/utils.
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=${CIS_RUNTIME_ROOT:-${SANDBOX_ROOT:-${PWD}/cis-runtime}}
ROOT_DIR=$(mkdir -p "${ROOT_DIR}" && cd -- "${ROOT_DIR}" && pwd)

# Keep a complete installation transcript while still showing progress on
# screen.  Set CIS_SANDBOX_LOG to choose another location.
SANDBOX_LOG=${CIS_SANDBOX_LOG:-${ROOT_DIR}/sandbox.log}
mkdir -p "$(dirname "${SANDBOX_LOG}")"
exec > >(tee -a "${SANDBOX_LOG}") 2>&1
printf '\n===== CIS sandbox start %s =====\n' "$(date '+%Y-%m-%d %H:%M:%S %Z')"

MESHTOOLS_REPO=${MESHTOOLS_REPO:-https://gitlab.awi.de/psong/meshtools_triangle.git}
MESHTOOLS_REF=${MESHTOOLS_REF:-feature/core3-paleo-workflow}
SPHERLAB_REPO=${SPHERLAB_REPO:-https://github.com/FESOM/spheRlab.git}
SPHERLAB_REF=${SPHERLAB_REF:-master}
FESOM2_REPO=${FESOM2_REPO:-https://github.com/FESOM/fesom2.git}
FESOM2_REF=${FESOM2_REF:-2.6.7}
PYFESOM2_VERSION=${PYFESOM2_VERSION:-}
CONDA_ENV_NAME=${CONDA_ENV_NAME:-meshtools-core3}
CONDA_ENV_PATH=${CONDA_ENV_PATH:-${ROOT_DIR}/conda-env}
R_LIBRARY=${R_LIBRARY:-${ROOT_DIR}/libs}

module load git
module load gcc/12.1.0
module load cdo
module load conda
module load netcdf-fortran/4.5.4-gcc12.1.0
module load netcdf-c/4.8.1-gcc12.1.0
module load r

nc_config=${NETCDF_FORTRAN_CONFIG:-$(command -v nf-config || true)}
if [[ -z "${nc_config}" ]]; then
    echo "ERROR: nf-config not found after loading NetCDF-Fortran" >&2
    exit 1
fi
nf_prefix=$(${nc_config} --prefix)
nc_c_config=${NETCDF_C_CONFIG:-$(command -v nc-config || true)}
if [[ -z "${nc_c_config}" ]]; then
    echo "ERROR: nc-config not found after loading NetCDF-C" >&2
    exit 1
fi
nc_c_prefix=$(${nc_c_config} --prefix)
export LD_LIBRARY_PATH="${nf_prefix}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export R_LIBS_USER="${R_LIBRARY}"

clone_if_missing() {
    local repository=$1 destination=$2 ref=$3
    if [[ -d "${destination}/.git" ]]; then
        printf 'Using existing repository: %s\n' "${destination}"
    elif [[ -e "${destination}" ]]; then
        printf 'ERROR: %s exists but is not a Git repository.\n' "${destination}" >&2
        return 1
    else
        git clone --branch "${ref}" "${repository}" "${destination}"
    fi
}

clone_if_missing "${MESHTOOLS_REPO}" "${ROOT_DIR}/meshtools_triangle" "${MESHTOOLS_REF}"

if [[ -f "${CONDA_ENV_PATH}/conda-meta/history" ]]; then
    printf 'Using existing Conda environment prefix: %s\n' "${CONDA_ENV_PATH}"
else
    conda env create --prefix "${CONDA_ENV_PATH}" \
        -f "${ROOT_DIR}/meshtools_triangle/environment.yml"
fi

# Some Conda activation hooks inspect optional variables without guarding
# against nounset mode.
set +u
conda activate "${CONDA_ENV_PATH}"
set -u

# CDO alone does not guarantee that ncks is installed.  NetCDF-C and
# NetCDF-Fortran are provided by the site modules above; `netcdf-c` is not a
# Conda-forge package name.  ncdump/ncrename therefore come from the loaded
# NetCDF-C module while ncks comes from Conda.
conda install --yes --channel conda-forge --prefix "${CONDA_ENV_PATH}" nco

python_exe="${CONDA_ENV_PATH}/bin/python"
export PYTHONNOUSERSITE=1
export PATH="${CONDA_ENV_PATH}/bin:${PATH}"
if [[ -n "${PYFESOM2_VERSION}" ]]; then
    "${python_exe}" -m pip install --no-user "pyfesom2==${PYFESOM2_VERSION}" 'setuptools<81' 'numpy<2' 'pandas<2.0'
else
    "${python_exe}" -m pip install --no-user pyfesom2 'setuptools<81' 'numpy<2' 'pandas<2.0'
fi

(
    cd "${ROOT_DIR}/meshtools_triangle"
    make allclean
    make -j4
    make check
    make test
    make test-python
)

clone_if_missing "${SPHERLAB_REPO}" "${ROOT_DIR}/spheRlab" "${SPHERLAB_REF}"
mkdir -p "${R_LIBRARY}"
Rscript -e 'install.packages("ncdf4", repos="https://cloud.r-project.org", lib=Sys.getenv("R_LIBS_USER"))'
R CMD INSTALL --library="${R_LIBRARY}" "${ROOT_DIR}/spheRlab"

clone_if_missing "${FESOM2_REPO}" "${ROOT_DIR}/fesom2" "${FESOM2_REF}"
(
    cd "${ROOT_DIR}/fesom2/mesh_part"
    # configure.sh is not idempotent on this checkout (it uses plain ln for
    # fesom_ini.x). Reconfigure only when the compiled target is absent.
    if [[ ! -x build/fesom_ini ]]; then
        bash -l configure.sh
    else
        printf 'Using existing FESOM mesh_part build.\n'
    fi
)

mkdir -p "${ROOT_DIR}/bin"
ln -sfn ../fesom2/mesh_part/build/fesom_ini "${ROOT_DIR}/bin/fesom_ini.x"
ln -sfn ../meshtools_triangle/bin/relax150614.x "${ROOT_DIR}/bin/relax150614.x"
ln -sfn ../meshtools_triangle/bin/relax150614.x "${ROOT_DIR}/bin/relax.x"
ln -sfn ../meshtools_triangle/bin/rg.x "${ROOT_DIR}/bin/rg.x"

# Write a sourceable runtime description for esm_tools jobs.
runtime_env="${ROOT_DIR}/runtime.env"
# Resolve the active Conda prefix rather than relying on `conda activate` being
# available in scheduler-launched shells.  All users can source this file and
# use the same absolute executables from this shared runtime root.
conda_prefix=$(${python_exe} -c 'import sys; print(sys.prefix)')
rscript_exe=$(command -v Rscript)
{
    printf 'export WETDRY_RUNTIME_ROOT=%q\n' "${ROOT_DIR}"
    printf 'export WETDRY_MESHTOOLS_BIN=%q\n' "${ROOT_DIR}/meshtools_triangle/bin"
    printf 'export WETDRY_FESOM_INI=%q\n' "${ROOT_DIR}/bin/fesom_ini.x"
    printf 'export R_LIBS_USER=%q\n' "${R_LIBRARY}"
    printf 'export CONDA_PREFIX=%q\n' "${conda_prefix}"
    printf 'export WETDRY_PYTHON=%q\n' "${python_exe}"
    printf 'export RSCRIPT=%q\n' "${rscript_exe}"
    # Keep site NCO/CDO first in model jobs; CIS Python uses WETDRY_PYTHON
    # explicitly and does not require the Conda bin directory on PATH.
    printf 'export PATH="%s:${PATH}"\n' "${ROOT_DIR}/bin:${ROOT_DIR}/meshtools_triangle/bin"
    # Keep private libraries available to coupling wrappers without changing the
    # model process' global loader path.  Albedo module libraries remain
    # authoritative for ECHAM/FESOM.
    # Native mesh tools are compiled against the setup-time NetCDF modules;
    # keep their libraries separate from the Python/Conda stack.
    printf 'export WETDRY_NATIVE_LD_LIBRARY_PATH=%q\n' "${CONDA_ENV_PATH}/lib:${nf_prefix}/lib:${nc_c_prefix}/lib"
    printf 'export WETDRY_PRIVATE_LD_LIBRARY_PATH=%q\n' "${CONDA_ENV_PATH}/lib"
} > "${runtime_env}"

printf 'Runtime file written: %s\n' "${ROOT_DIR}/runtime.env"
printf 'Set WETDRY_QGLB_DIR, WETDRY_PRIORITY_NODES and PYFESOM_PATH before checking.\n'
# Validate the runtime just written, rather than stale WETDRY_* variables that
# may have been exported by an older sandbox in the calling shell.
# shellcheck disable=SC1090
source "${runtime_env}"
"${SCRIPT_DIR}/cis_runtime_check.sh" "${ROOT_DIR}" || {
    echo "Runtime installation completed, but preflight reported missing experiment inputs or dependencies." >&2
    exit 1
}
printf 'CIS wet/dry sandbox setup completed successfully.\n'
printf 'Sandbox log: %s\n' "${SANDBOX_LOG}"
