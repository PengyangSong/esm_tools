# AWIESM – PISM – VILMA three-model iterative coupling

Runscripts for coupling AWIESM2 (climate) with PISM (ice sheet) and VILMA (solid earth).

## CIS wrapper convention

All coupling wrappers use the **CIS** (Climate-Icesheet-Solidearth) naming:
- `couple_in` calls `cis2<model>`, `couple_out` calls `<model>2cis`
- Each wrapper sources and calls the individual coupling functions sequentially

## File structure

```
awiesm_pism_vilma.yaml   <- top-level iterative coupling controller
├── awiesm.yaml          <- model1: climate (ECHAM + FESOM via OASIS)
├── pism.yaml            <- model2: ice sheet
└── vilma.yaml           <- model3: solid earth (GIA)
```

## Execution order per chunk

1. **AWIESM** (model1): ECHAM + FESOM run for `chunk_size` years.
   - `echam2cis` → `echam2ice` + `echam2solidearth` (placeholder)
   - `fesom2cis` → `fesom2ice` + `fesom2solidearth` (placeholder)
2. **PISM** (model2): ice sheet run for its `chunk_size`.
   - `cis2pism` → `solidearth2pism` (bedrock update from VILMA) + `awiesm2pism` (climate forcing)
   - `pism2cis` → `pism2awiesm` (orography, freshwater for ECHAM) + `pism2solidearth` (ice thickness for VILMA)
3. **VILMA** (model3): solid earth run for its `chunk_size`.
   - `cis2vilma` → `ice2vilma` (ice thickness from PISM) + `awiesm2vilma` (placeholder)
   - `vilma2cis` → `vilma2ice` (bedrock change for next PISM chunk) + `vilma2awiesm` (placeholder)

## Key settings to adjust

All three runscripts must share the **same** `general.base_dir` (and `expid`) so `experiment_couple_dir` is consistent.

- **awiesm_pism_vilma.yaml**: `base_dir`, `account`, chunk sizes for each model.
- **awiesm.yaml**: `version`, `model_dir`, `mesh_dir`, `ini_parent_*`, dates.
- **pism.yaml**: `model_dir`, `pool_dir`, `spinup_file`, `domain`, `resolution`, PISM physics, `coupled_to_solidearth`.
- **vilma.yaml**: `model_dir`, `scenario`, dates.
