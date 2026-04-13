#!/usr/bin/env Rscript
# NetCDF mesh description for FESOM (from gen_wetdry_v2).
# TODO: remove hard-coded lib.loc; align with gen_wetdry_v2/Step4_fesom_ini/gen_mesh_cst.R / site config.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: gen_mesh_cst.R <mesh_dir>")
}
meshpath <- args[1]
ofile <- paste0(meshpath, "/mesh_cst.nc", sep = "")

library(spheRlab, lib.loc = "/albedo/home/psong/R/x86_64-conda-linux-gnu-library/4.2")
library(ncdf4, lib.loc = "/albedo/home/psong/R/x86_64-conda-linux-gnu-library/4.2")

rotated <- FALSE
if (rotated) {
  grid <- sl.grid.readFESOM(
    griddir = meshpath, rot = rotated, rot.invert = rotated,
    rot.abg = c(50, 15, -90)
  )
} else {
  grid <- sl.grid.readFESOM(
    griddir = meshpath, rot = rotated, rot.invert = rotated,
    rot.abg = c(0, 0, 0), threeD = FALSE
  )
}

sl.grid.writeCDO(grid, ofile = ofile, netcdf = TRUE, depth = FALSE)
