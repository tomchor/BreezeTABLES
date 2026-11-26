# BreezeTABLES

Breeze.jl simulations for the TABLES large-eddy simulation intercomparison, organized by A. Rauchöcker, S. Serafin and I. Stiperski.

## What to follow

- `Rauchoecker_QJRMS_TABLES.pdf` is the reference for the model setup and the output. It is the compiled version of `ModelIntercomparison_ModelSetupAndOutput.tex` in the organizers' Overleaf project. If the two ever differ, the Overleaf source is the most recent version. The PDF in this repository has the same text as the version compiled from the Overleaf source on 26 November 2026.
- The initial profiles are the six sounding files in `Rauchoecker_QJRMS_TABLES/ModelSetUps/TABLES_soundings/`.
- The manuscript draft is not a reference for the setup. Its text described the inversion wrongly, and the organizers have since corrected it.

## Confirmed by the organizers

A. Rauchöcker confirmed the following by email on 26 November 2026.

- The flat and 100 m ridge cases start with a 4 K inversion on top of the mixed layer, as the protocol says. In the sounding files it lies between 1000 m and 1050 m.
- The surface heat flux of 0.12 K m/s is per horizontal area, not per surface area. The total heat input is then the same in all cases. `tables.jl` already works this way, because Oceananigans applies a bottom flux through the horizontal area of the cell, also on terrain-following grids.
- We may interpolate to the common grid with our own code (issue #4). The routine in `Rauchoecker_QJRMS_TABLES/interpolation_routine/` is the one the organizers used for CM1 output on terrain-following coordinates. It shows the conventions to match, and using it is optional.

## Model table of the manuscript

The organizers asked for a row for Breeze in the model table of the manuscript (`tab:model_list` in `ModelIntercomparison_Manuscript.tex` of the Overleaf project), describing the current setup. The row has to be updated when the wall model (issue #2) or the closure (issue #3) changes. The current setup is:

| Model name | Eqns. | Space discretization | Terrain implementation | SGS | Wall model |
|---|---|---|---|---|---|
| Breeze | Compress. | finite volume | terrain following | none (implicit, fifth-order WENO) | neutral log law |
