#!/usr/bin/env Rscript
# NetCDF mesh description for FESOM (from gen_wetdry_v2).
# The runtime selects the R library through R_LIBS_USER/R_LIBS_SITE.  Do not
# embed a developer-specific library path here.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: gen_mesh_cst.R <mesh_dir>")
}
meshpath <- args[1]
ofile <- paste0(meshpath, "/mesh_cst.nc", sep = "")

library(spheRlab)
library(ncdf4)

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
