# AWIESM – PISM two-model iterative coupling (without VILMA)

Runscripts for coupling AWIESM2 (climate) with PISM (ice sheet), without solid earth feedback.
This template uses the same CIS coupling wrappers as `awiesm-pism-vilma` but with `coupled_to_solidearth: 0`,
so `solidearth2pism` and `pism2solidearth` are skipped at runtime.

## File structure

```
awiesm_pism.yaml   <- top-level iterative coupling controller
├── awiesm.yaml    <- model1: climate (ECHAM + FESOM via OASIS)
└── pism.yaml      <- model2: ice sheet (coupled_to_solidearth: 0)
```

## Differences from awiesm-pism-vilma

- No `model3` (VILMA) in the controller
- `pism.coupled_to_solidearth: 0` disables solid earth coupling functions
- No `solidearth_initialize_method` needed (no bedrock update)

## Key settings to adjust

Both runscripts must share the **same** `general.base_dir` (and `expid`).

- **awiesm_pism.yaml**: `base_dir`, `account`, chunk sizes.
- **awiesm.yaml**: `version`, `model_dir`, `mesh_dir`, `ini_parent_*`, dates.
- **pism.yaml**: `model_dir`, `pool_dir`, `spinup_file`, `domain`, `resolution`, PISM physics.
