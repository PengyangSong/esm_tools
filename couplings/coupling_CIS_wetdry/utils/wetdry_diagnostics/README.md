# Wet/dry FESOM diagnostics

This directory contains offline diagnostics for a FESOM wet/dry mesh transition.
It does not modify the coupling functions or the FESOM executable.

The first tool, `inspect_blowup.py`, reads a FESOM blow-up NetCDF file and the
mesh used for that run. It reports the nodes with non-finite or extreme tracer
values and prints the barotropic fields (`ssh_rhs`, `ssh_rhs_old`, `d_eta`),
vertical CFL, layer thickness, and neighboring nodes.

Example:

```bash
PYTHONNOUSERSITE=1 "$WETDRY_PYTHON" inspect_blowup.py \
  /path/to/fesom.2026.oce.blowup.nc \
  /path/to/fesom_wetdry_submesh_20260101-20301231 \
  --radius 1
```

The diagnostic uses only data already written by FESOM and is intended to
identify the first unstable columns before any model-code instrumentation.

To inspect approximate edge transports around the unstable cluster:

```bash
PYTHONNOUSERSITE=1 "$WETDRY_PYTHON" inspect_edge_transport.py \
  /path/to/fesom.2026.oce.blowup.nc \
  /path/to/fesom_wetdry_submesh_20260101-20301231 \
  --node 27205 27206 27213 27214 27215
```

The edge metric is a local tangent-plane approximation; use it to identify
problematic edges and cavity-level discontinuities, not as a bit-for-bit
reproduction of FESOM's spherical metric.

To check element/node cavity consistency:

```bash
PYTHONNOUSERSITE=1 "$WETDRY_PYTHON" check_cavity_consistency.py \
  /path/to/fesom.2026.oce.blowup.nc \
  /path/to/fesom_wetdry_submesh_20260101-20301231 \
  --node 27205 27206 27213 27214 27215
```
