# Files from the TABLES organizers

This folder holds files from the project folder that the organizers of the TABLES intercomparison (A. Rauchöcker, S. Serafin and I. Stiperski) share with the participants. They are unmodified copies, as received on 24 November 2026. The protocol itself is `Rauchoecker_QJRMS_TABLES.pdf` in the root of this repository.

## Contents

- `ModelSetUps/TABLES_soundings/`: the initial profiles of the six cases (Section 2 of the protocol) and a figure of them.
- `interpolation_routine/`: the organizers' Python routine that renames the output variables and interpolates them to the common grid (Section 5 of the protocol), with its documentation.

The rest of the project folder (manuscript draft, figures, meeting notes, LaTeX files) is not tracked. The `.gitignore` of this folder excludes it.

## Sounding files

Lines that start with `|` are comments. The first line of data holds the surface pressure (hPa), potential temperature (K) and water vapor mixing ratio (g/kg). Each line after that holds one level: height (m), potential temperature (K), water vapor mixing ratio (g/kg), u (m/s) and v (m/s). Levels are 10 m apart and start at 10 m. They end at 10000 m in `u10_flat` and `u10_ridge100`, and at 6000 m in the other files.

| File | Potential temperature | Wind |
|---|---|---|
| `input_sounding_teamx_u00_flat` | 300 K up to 1000 m, then 0.08 K/m up to 1050 m (304 K), then 0.003 K/m | none |
| `input_sounding_teamx_u10_flat` | same | u = 10 m/s |
| `input_sounding_teamx_u00_ridge100` | same | none |
| `input_sounding_teamx_u10_ridge100` | same | u = 10 m/s |
| `input_sounding_teamx_u00_ridge1000` | 300 K at the surface, 0.003 K/m | none |
| `input_sounding_teamx_u10_ridge1000` | same | none up to 1000 m, linear to u = 10 m/s at 1400 m, uniform above |

The surface values are 1000 hPa and 300 K in all files. Water vapor and v are zero everywhere. The files of the flat cases are identical to those of the 100 m ridge cases.

The analytic profiles in `tables.jl` reproduce the six files at every level.

## Conventions set by the interpolation routine

- The target levels are `np.arange(5 + 1e-6, 3000, 10)`: 300 levels from 5 m to 2995 m. For the ridge cases they are heights above the foot of the ridge, since the routine interpolates from the true height of each cell (`zhval` in CM1).
- The interpolation is linear in the vertical and is done column by column. There is no horizontal interpolation.
- Target levels outside the range of the model levels of a column become NaN, the default of `scipy.interpolate.griddata`. This includes the points below the terrain.
- The routine is applied to the finished statistics. A comment in the script recommends moving staggered variables to cell centers beforehand.
- In the variable mapping for CM1 (`config_interpolation.py`), the snapshot holds the velocities interpolated to cell centers (`uinterp`, `vinterp`, `winterp`), the potential temperature (`th`) and the pressure (`prs`, which is the full pressure in CM1).
- The script handles the statistics of the ridge cases as (z, x) fields and those of the flat cases as profiles in z. The protocol asks for the dimensions (x, z) in all cases.
- The README of the routine gives the vertical range of the common grid as [5, 4995] m. The script and the protocol use 3000 m.
