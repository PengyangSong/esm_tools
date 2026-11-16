# CIS wet/dry coupling

## Purpose

`coupling_CIS_wetdry` implements offline iterative coupling between AWI-ESM2
(ECHAM6, FESOM2, JSBACH and OASIS3-MCT), PISM, and optionally VILMA. Its
distinguishing feature is rebuilding the FESOM wet/dry mesh when PISM changes
the ice-sheet and ice-shelf geometry.

```text
ECHAM  -- atmosphere forcing --------> PISM
ECHAM  <-- orography, mask, runoff ---- PISM

FESOM  -- cavity temperature/salinity
          or prescribed basal melt ---> PISM
FESOM  <-- geometry and freshwater ----- PISM
                    |
                    `-- new wet/dry mesh and remapped FESOM restart

PISM   -- ice load -------------------> VILMA
PISM   <-- bedrock change -------------- VILMA
```

The models run sequentially in independently scheduled chunks and exchange
NetCDF files through `general.experiment_couple_dir`.

## Package layout

```text
coupling_CIS_wetdry/
|-- echam/       ECHAM couple-in/out wrappers and field processing
|-- fesom/       FESOM coupling, dynamic mesh update and restart remapping
|-- general/     shared shell, CDO, file-list and PISM helper functions
|-- pism/        atmosphere/ocean/solid-Earth coupling for PISM
|-- utils/       conversion and preprocessing programs
|   `-- wetdry/  wet/dry mesh and restart utilities
`-- vilma/       PISM-VILMA coupling and AWI-ESM placeholders
```

Each model directory normally contains:

- `coupling_cis2<model>.functions`: top-level `couple_in` wrapper;
- `coupling_<model>2cis.functions`: top-level `couple_out` wrapper;
- directional coupling implementations;
- `env_<model>.py`: conversion from esm_tools configuration to the shell
  environment used by the coupling functions.

## Iterative scheduling

The top-level iterative controller declares the participating setups and their
chunk lengths. For example:

```yaml
general:
    iterative_coupling: true

model1:
    setup_name: awiesm
    runscript: awiesm.yaml
    chunk_size: 5
    chunk_unit: years

model2:
    setup_name: pism
    runscript: pism.yaml
    chunk_size: 10
    chunk_unit: years
```

Individual model runscripts schedule coupling jobs around the ordinary esm_tools
workflow:

```yaml
workflow:
    next_run_triggered_by: couple_out
    subjobs:
        couple_in:
            run_before: prepcompute
            run_only: first_run_in_chunk
        couple_out:
            run_after: tidy
            run_only: last_run_in_chunk
```

The important runtime state includes `first_run_in_chunk`,
`last_run_in_chunk`, `chunk_number`, `chunk_start_date`, and `chunk_end_date`.

## PISM coupling

### PISM couple-in

The entry point is `pism/coupling_cis2pism.functions`:

```text
cis2pism
|-- solidearth2pism
`-- awiesm2pism
    |-- atmosphere2pism
    `-- ocean2pism
```

The principal outputs consumed by PISM are:

```text
latest_atmo_forcing_file.nc
latest_ocean_forcing_file.nc
```

Solid-Earth coupling is controlled by `pism.coupled_to_solidearth`. Set it to
`0` for the two-model AWI-ESM-PISM configuration.

### Atmosphere to PISM

`pism/coupling_atmosphere2pism.functions` converts ECHAM output into PISM
surface forcing. Depending on the selected method it processes temperature,
precipitation, radiation, humidity, cloud cover, wind, emissivity and
transmissivity.

The workflow supports elevation-dependent downscaling, interpolation and
conservative remapping. Surface mass balance can be prepared through PDD or
DEBM and related energy-balance pathways.

Representative configuration is:

```yaml
pism:
    regrid_method: DOWNSCALE
    lapse_rate: -0.005
    downscale_precip: 1
    ablation_method: DEBM
```

### FESOM to PISM

The first stage is executed during FESOM `couple_out` by
`fesom/coupling_fesom2ice.functions`. It expects annual cavity-interface files:

```text
cavity_Tsurf.fesom.<year>.nc
cavity_Ssurf.fesom.<year>.nc
cavity_Bmelt.fesom.<year>.nc
```

It verifies the requested years, converts the unstructured fields to a
horizontal representation with longitude and latitude, optionally applies the
`GLAC_flag` mask, and creates:

```text
ocean_file_for_ice.nc
ocean.griddes
ocean_names_for_ice.dat
```

During PISM `couple_in`, `pism/coupling_ocean2pism.functions` temporally
processes and remaps these fields to the PISM grid. Two canonical modes are
supported.

#### Temperature-salinity mode

```yaml
pism:
    fesom_to_pism_ocean: temperature_salinity
```

Aliases include `ts`, `tempsalt`, `th`, and `ocean_th`.

```text
FESOM Tsurf --> PISM theta_ocean
FESOM Ssurf --> PISM salinity_ocean
PISM option: -ocean th
```

#### Prescribed basal-melt mode

```yaml
pism:
    fesom_to_pism_ocean: basal_melt
```

Aliases include `basal`, `basal_melt_rate`, `shelfbmassflux`, `given`, and
`direct`.

```text
FESOM Tsurf      --> PISM shelfbtemp
FESOM basal melt --> PISM shelfbmassflux
PISM option: -ocean given
```

The coupling writes the necessary esm_pism option overrides so that PISM uses
the ocean model matching the generated forcing file. Removed combined or
transition modes fail explicitly rather than silently selecting a pathway.

Supported spatial processing includes `INTERPOLATE`, `REMAP`, `EXTRAPOLATE`,
and `NONE`. The usual CDO interpolation type is distance weighted (`dis`).

### PISM couple-out

The entry point is `pism/coupling_pism2cis.functions`:

```text
pism2cis
|-- pism2awiesm
|   |-- pism2atmosphere
|   `-- pism2ocean
`-- pism2solidearth
```

`pism2atmosphere` generates an updated domain mask, ice-thickness and
orography anomalies, and glacial runoff/discharge for ECHAM and JSBACH.

`pism2ocean` can prepare basal shelf freshwater, iceberg freshwater or heat,
iceberg-model input, basal temperature-gradient fields and geometry masks. It
also creates the wet/dry geometry product:

```text
ice_file_at_ocean.combined.nc
```

Freshwater conversions use seawater density. Optional residual freshwater
accounting is controlled with:

```yaml
pism:
    residual_fw_to_ocean: false
    residual_fw_negative_tolerance: 1.0e-14
```

## FESOM wet/dry coupling

### FESOM couple-in

The entry point is `fesom/coupling_cis2fesom.functions`:

```text
cis2fesom
|-- ice2fesom
`-- solidearth2fesom
```

The solid-Earth-to-FESOM path is currently a placeholder. The wet/dry workflow
is implemented in `fesom/coupling_ice2fesom.functions`.

### Wet/dry sequence

```text
ice2fesom
|-- freeze the previous restart-mesh identity
|-- generate a new wet/dry submesh
|   |-- cut maximum mesh and submesh
|   |-- generate node and element mappings
|   |-- derive node and cavity depths
|   |-- enforce minimum cavity layers
|   |-- repair neighboring-column vertical overlap
|   |-- validate all mesh edges
|   `-- initialize and partition the mesh
|-- remap the previous FESOM restart to the new mesh
|-- publish the new active mesh
|-- persist its physical path and horizontal node count
|-- optionally prepare iceberg discharge
`-- remap PISM forcing to the FESOM grid
```

Enable this behavior explicitly:

```yaml
fesom:
    iterative_wetdry_mesh: true
```

### Mesh identity and state

The workflow distinguishes four mesh roles:

1. `fesom.mesh_dir`: configured initial mesh and bootstrap fallback;
2. `WETDRY_LAST_MESH_DIR`: mesh on which the input FESOM restart was written;
3. `WETDRY_SUBMESH_DIR`: newly generated mesh for the current handoff;
4. `WETDRY_SUBMESH_ACTIVE_DIR`: link to the most recently published mesh.

Default chunk-specific paths are based on a `YYYYMMDD-YYYYMMDD` chunk tag:

```text
${COUPLE_DIR}/fesom_wetdry_submesh_<chunk-tag>
${COUPLE_DIR}/fesom_wetdry_restart_remap_<chunk-tag>
${COUPLE_DIR}/fesom_wetdry_partition_work_<chunk-tag>
```

The active mesh and persistent state are represented by:

```text
${COUPLE_DIR}/fesom_wetdry_submesh_active
${COUPLE_DIR}/wetdry_last_mesh_dir.txt
${COUPLE_DIR}/wetdry_last_mesh_nx.txt
```

The active-mesh link is changed only after restart remapping has used the old
mesh. This publication order prevents the new mesh from accidentally being
used as the source description of the old restart.

### Mesh stages

#### 1. Submesh cutting

`iterative_coupling_ice_fesom_wetdry_submesh_cutting` prepares the maximum
mesh, incorporates PISM geometry and runs mesh cutting and reordering.

Required commands are `relax150614.x` and `rg.x`. Their directory may be
provided through `WETDRY_MESHTOOLS_BIN`.

During cutting, `wetdry_maxmesh_nodhn.py` transfers PISM shelf geometry to the
QGLB nodes. Continuous fields use mask-aware bilinear interpolation; the
PISM latitude coordinate is sorted numerically, so descending latitude (the
orientation commonly seen as `InvP` in ncview) is handled correctly. The
grounded and floating masks remain categorical. For genuinely curvilinear
2-D PISM lon/lat coordinates, the same continuous transfer uses a local
Cartesian-on-the-sphere stencil. A limited floating-shelf-only graph smoother
is applied by default to reduce one-cell shelf-base steps.

Controls are:

```text
WETDRY_GEOMETRY_INTERPOLATION   bilinear (default) or nearest
WETDRY_GEOMETRY_SMOOTH_PASSES   default 1; set 0 to disable
WETDRY_GEOMETRY_SMOOTH_BLEND    default 0.5
WETDRY_FLOATING_MASK_VAR        default floating_ice_mask
```

The intermediate fields are written under the maximum-mesh directory as
`wetdry_shelf_on_maxmesh_raw.out`, `wetdry_shelf_on_maxmesh_smoothed.out`,
and `wetdry_geometry_interp_diagnostics.txt`.

#### 2. Submesh mapping

`iterative_coupling_ice_fesom_wetdry_submesh_mapping` runs
`utils/wetdry/submesh_mapping.py` and produces:

```text
map_nod.out
map_elem.out
```

#### 3. Depth generation

`iterative_coupling_ice_fesom_wetdry_submesh_depth_from_map` uses
`utils/wetdry/wetdry_submesh_depth_outputs.py` to combine maximum-mesh
bedrock, mapped PISM shelf geometry and FESOM vertical levels. Representative
outputs are `depth@node.out`, `cavity_depth@node.out`, and `aux3d.out`.

#### 4. Cavity validation

`iterative_coupling_ice_fesom_wetdry_submesh_cavity_min_layers` uses
`utils/wetdry/wetdry_submesh_cavity_min_layers.py` to ensure that new ocean
cavities contain enough vertical layers.

Relevant controls include:

```text
WETDRY_CAVITY_MIN_LAYER_SPAN
WETDRY_CAVITY_BED_INDEX_OFFSET
WETDRY_CAVITY_DEPTH_ADJUST_EPS
WETDRY_CAVITY_MIN_LAYERS_THIN_ONLY
WETDRY_CAVITY_DEPTH_POSITIVE_DOWN
WETDRY_CAVITY_ZMID_CONVENTION
```

After the per-node minimum-layer adjustment, the workflow applies the
FESOM--ISSM-style edge repair in
`utils/wetdry/wetdry_repair_edge_vertical_overlap.py`. For every edge with no
shared active vertical layer, it deepens the shallower ocean-bed column to the
next layer boundary reached by its deeper neighbor. The shelf-base field is
not changed. The repair is bounded by `WETDRY_EDGE_REPAIR_MAX_DEEPENING`
(default 500 m), followed by a complete edge check. Diagnostics are written
to `edge_overlap_repair_diagnostics.txt`; the pre-repair depth is retained as
`depth@node.out.before_edge_repair`.

```text
WETDRY_EDGE_REPAIR_MAX_DEEPENING   default 500 m
WETDRY_EDGE_REPAIR_MAX_PASSES      default 10
WETDRY_SKIP_EDGE_VERTICAL_OVERLAP_REPAIR=1   disable (not recommended)
```

#### 5. Initialization and partitioning

`iterative_coupling_ice_fesom_wetdry_submesh_partition` runs `fesom_ini.x`,
creates mesh constants, assigns ocean and glacier flags, and builds the
partition information for the next FESOM run.

The executable can be supplied through `FESOM_INI` or `WETDRY_FESOM_INI`.
Supporting resources are under `utils/wetdry/submesh_partition/`.

### FESOM restart remapping

`iterative_coupling_ice_fesom_modify_fesom_restart` invokes the utilities
`utils/wetdry/fesom2_restart.py` and
`utils/wetdry/split_FESOM2_restart.py` to transfer the previous FESOM state
from `WETDRY_LAST_MESH_DIR` to `WETDRY_SUBMESH_DIR`.

The restart year defaults to the year before the current chunk start and may be
overridden with `WETDRY_FESOM_RESTART_YEAR`.

The central consistency requirement is:

> `WETDRY_LAST_MESH_DIR` must identify the physical mesh on which the input
> FESOM restart was actually written.

### OASIS grid-size synchronization

A wet/dry update can change the number of FESOM surface nodes. The new count is
stored in `wetdry_last_mesh_nx.txt`. Before OASIS configuration is generated,
esm_tools applies the value to both:

```yaml
fesom:
    nx: <node count>
    grids:
        feom:
            nx: <node count>
```

At wet/dry chunk boundaries, old OASIS restart inputs are excluded so OASIS
regenerates interpolation state for the changed FESOM grid. This special
handling is scoped to runs with both `general.iterative_coupling: true` and
`fesom.iterative_wetdry_mesh: true`.

## ECHAM coupling

The ECHAM couple-out entry point is `echam2cis`:

```text
echam2cis
|-- echam2ice
`-- echam2solidearth
```

`echam2ice` selects and concatenates the atmospheric variables and years needed
by PISM and may produce multi-year means. `echam2solidearth` is currently a
placeholder.

The ECHAM couple-in entry point is `cis2echam`:

```text
cis2echam
|-- ice2echam
`-- solidearth2echam
```

`ice2echam` can update ECHAM orography and glacial masks, JSBACH initialization
and restart state, vegetation ratios, glacial runoff and the corresponding
namelist or OASIS inputs. `solidearth2echam` is currently a placeholder. The
first climate chunk normally skips couple-in because there is no previous PISM
update.

## VILMA coupling

`vilma/coupling_ice2vilma.functions` remaps PISM ice thickness and masks to the
VILMA grid and prepares load-history input. `coupling_vilma2ice.functions`
prepares bedrock change for application to PISM's `topg` restart field.

The AWI-ESM-to-VILMA and VILMA-to-AWI-ESM functions are currently explicit
placeholders. Future implementations may handle relative sea level, land-sea
masks, bathymetry and coastal geometry.

## Representative workflow configuration

### Per-user runtime (native, no container)

The coupling can use a user-local runtime such as the `sandbox.sh` bootstrap
installation. Point the runscript at that installation instead of relying on
developer-specific defaults:

```yaml
fesom:
    iterative_wetdry_mesh: true
    wetdry_runtime_root: /path/to/user/work/test
    wetdry_qglb_dir: /path/to/CORE3_QGLB
    wetdry_priority_nodes: /path/to/priority_nodes.PISM.txt
    # Optional overrides when they are not below wetdry_runtime_root:
    # wetdry_rscript: /path/to/Rscript
    # wetdry_r_library: /path/to/R/library
    # wetdry_python: /path/to/python
    # pyfesom_path: /path/to/pyfesom2
```

With `wetdry_runtime_root` set, the environment preparation derives the mesh
binary directory as `<root>/meshtools_triangle/bin`, FESOM initialization as
`<root>/bin/fesom_ini.x`, and the R library as `<root>/libs`. The QGLB mesh,
priority-node file, R executable, and Python/`pyfesom2` installation remain
explicit because they are experiment or site data rather than generic tools.

The runtime must provide `relax150614.x`, `rg.x`, `fesom_ini.x`, R with
`spheRlab` and `ncdf4`, Python with `pyfesom2`, and CDO/NCO/NetCDF utilities.
Run the sandbox's package checks after installation and ensure the same
environment is loaded in the esm_tools batch jobs.

Before launching a coupling job, validate the installation with:

```bash
source /path/to/cis_runtime.env
couplings/coupling_CIS_wetdry/utils/cis_runtime_check.sh \
    "$WETDRY_RUNTIME_ROOT" "$WETDRY_QGLB_DIR" "$WETDRY_PRIORITY_NODES"
```

The template `utils/cis_runtime.env.example` lists the runtime variables. The
checker is deliberately independent of the model workflow and should pass
before the first `couple_in` job is submitted.

```yaml
fesom:
    iterative_wetdry_mesh: true
    choose_general.chunk_number:
        1:
            mesh_dir: "${general.start_mesh}"
        "*":
            mesh_dir: "${general.update_mesh}"
    number_of_years_for_forcing: 2
    workflow:
        next_run_triggered_by: couple_out
        subjobs:
            couple_in:
                run_before: prepcompute
                script: coupling_cis2fesom.functions
                script_dir: ${general.script_dir}/coupling_CIS_wetdry/fesom
                call_function: cis2fesom
                env_preparation: env_fesom.py
                run_only: first_run_in_chunk
                skip_chunk_number: 1
            couple_out:
                run_after: tidy
                script: coupling_fesom2cis.functions
                script_dir: ${general.script_dir}/coupling_CIS_wetdry/fesom
                call_function: fesom2cis
                env_preparation: env_fesom.py
                run_only: last_run_in_chunk
```

```yaml
pism:
    fesom_to_pism_ocean: temperature_salinity
    residual_fw_to_ocean: false
    coupled_to_solidearth: 0
    couplers:
        surface:
            given:
                files:
                    surface_given_file: latest_atmo_forcing_file.nc
        ocean:
            th:
                files:
                    ocean_th_file: latest_ocean_forcing_file.nc
```

## Dependencies

The workflow uses Korn shell, Python, CDO, NCO, NetCDF utilities, R and several
compiled tools. Wet/dry-specific runtime commands include:

```text
relax150614.x
rg.x
fesom_ini.x
Rscript
python or python3
cdo
ncks
ncrename
ncdump
```

Some current defaults are site- or user-specific absolute paths. Portable
deployments should provide these tools through configuration, environment
modules, or a container and should record their versions in experiment
provenance.

## Important coupling products

```text
latest_atmo_forcing_file.nc
latest_ocean_forcing_file.nc
ocean_file_for_ice.nc
ocean.griddes
ocean_names_for_ice.dat
ice_file_at_ocean.combined.nc
ice.griddes
ice_names_for_ocean.dat
wetdry_last_mesh_dir.txt
wetdry_last_mesh_nx.txt
fesom_wetdry_submesh_active
fesom_wetdry_submesh_<chunk-tag>/
fesom_wetdry_restart_remap_<chunk-tag>/
```

## Failure conditions

The coupling stops or cannot complete when, for example:

- required FESOM cavity or PISM output files are missing;
- a forcing mode is unsupported;
- required mesh executables or `fesom_ini.x` are unavailable;
- mapping files or `nod2d.out` are missing;
- NetCDF variables or dimensions do not match expectations;
- the input restart is paired with the wrong source mesh;
- restart remapping or mesh partitioning fails.

## Current limitations

- Some executable defaults are tied to the Albedo environment and personal
  paths.
- General CIS code is duplicated between `coupling_CIS` and
  `coupling_CIS_wetdry`.
- Mesh generation, restart remapping and esm_tools orchestration are tightly
  coupled in `coupling_ice2fesom.functions`.
- Several compiled or generated artifacts are stored in the coupling tree.
- AWI-ESM/VILMA geometry feedback remains a placeholder.
- The branch does not contain a tracked end-to-end automated wet/dry test.

Despite these limitations, the implemented handoff is complete in its central
purpose: PISM changes the geometry, a compatible FESOM mesh and restart are
generated, OASIS dimensions are refreshed, and the following climate chunk can
run on the updated ocean domain.
