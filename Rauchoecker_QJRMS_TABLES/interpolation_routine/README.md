# TABLES data+postprocessing

Dear colleagues,

This repository provides a summary of the Research questions we want to work on, necessary variables to answer them and information on the data format and common grid of the output files.

In addition, it includes a function and script for the following tasks for already postprocessed data:
* renaming the output variables according to the naming convention below
* interpolating the model output to a common cartesian grid in the vertical

The script requires numpy, scipy and xarray to be installed. The environment file attached installs these dependencies (and all others I need for postprocessing during my work) by running

```python
conda install -f xreinstats.yml
```

## Usage of interpolation file

### Rename variables
The file "config_interpolation.py" defines the name of the desired variables in the model-specific output files and maps them to the variable name according to the naming convention with a dictionary; one each for the averaged statistics, snapshot variables and spectra. 
Import the desired dictionary to the file "interpolation_tables.py" to use it to rename variables, interpolate them to the common grid or both.
If your variables are already named accordingly, you can comment the relevant line out.

### Interpolating
Adapt the following variables to fit your data:

* model: the model used
* name: name of at least one researcher of the group (in case more groups work with the same model)
* case: "flat_u00", "flat_u10", "ridge100_u00", "ridge100_u10", "ridge1000_u00" or "ridge1000_u10"
* content: "statistics" or "snapshot"

These are needed to name the output file according to the agreed format. In addition, the following are needed:

* filename
* path
* outpath

The script will automatically interpolate to the new vertical grid. 
Most models define their coordinate at the mass level (i.e. the center of the grid cell / volume), which for the ~10 m grid spacing we agreed upon would place them 5 m, 15 m, 25 m, etc above the ground.  
Due to different grid necessities and terrain solutions in the models (terrain following coordinates vs. immersed boundaries, stretching), the vertical coordinate will differ for most models and an interpolation is inevitable.
To limit the effect of that interpolation for codes that use an isotropic grid spacing close to the surface, the vertical coordinate for the interpolation will therefore be z = [5, 4995] m.
Because tests revealed that the lowest model sometimes was below the lowest level of the target grid due to rounding errors, a small offset of 10$^{-6}$ was added to the coordinate levels of the target grid to compensate for that effect.

The grid spacing in the horizontal directions should be 20 m for all models and we want to limit interpolations as much as possible. 
Also the horizontal coordinate is often defined at the mass level (i.e. 10 m, 30 m, 50 m, etc horizontally for 20 m grid spacing), but we were not sure whether all models involved would do that.
Hence, we decided to deduct the minimum value of the horizontal coordinates from each coordinate value to have all data at the same location irrespective of the definition of individual models (at 0 m, 20 m, 40 m, etc):
$$
x_{out} = x_{model} - min(x_{model})
$$ 
$$
y_{out} = y_{model} - min(y_{model})
$$

## Research questions

Here is a summary of the questions we came up with:

1. How does the variability in profiles of mean and turbulent quantities between the different models compare to the variability within the turbulent flow itself?
    * Or in other words, how well are LES results comparable across different codes in relation to the internal variability of the flow itself?
    * Here the variance from the ensemble average from all models can be used as a measure of inter-model variability or model uncertainty.
    * The variance of each variable around the average for each individual model is the natural variability inherent to the turbulent flow.
    * The focus will be on the vertical structure at characteristic locations along the terrain: valley floor, ridge top, center of upstream and downstream slope.
    * Here we will also explore the Reynolds stress anisotropy (Stiperski et al. 2019, Stiperski et al. 2021, Stiperski and Calaf 2023) as well as scaled variances (surface layer scaling, mixed layer scaling) to test the variability in scaling relations. 
2. How does the CBL structure over mountains differ from the flat terrain CBL, and do all models simulate the same structure?
    * It is known from lower resolution LES simulations (e.g. Wagner et al. 2015, Leukauf et al. 2015) that the mountain CBL does not feature a pronounced entrainment layer with a negative sensible heat flux. Is this feature equally present in all the models? 
    * Additionally, for tall mountains, multiple stacked thermally driven circulation cells (consisting of upslope flow detaching from the sides with subsidence over the center of the valley) are known to occur in low wind conditions. How well are these cells observed and is the location of separation as well as the strength of the upslope flows and subsidence over the center of the valley equally reproduced in all the models?
3. Are there significant differences between turbulence statistics simulated by fully compressible models and incompressible models / models that use the Boussinesq approximation for simulations of a high ridge?
    *	he Boussinesq approximation is valid for small vertical displacements, but might not be accurate enough to represent large vertical displacements in deep convective boundary layers over tall mountains. Here we will explore the performance of all the models over the tall ridge. 
4. For the same forcing, are there differences in coherent structures observed over flat and mountainous terrain? If so, are these differences comparable between the models?
    * It is well established that the coherent structures over flat terrain undergo a transition from rolls to cells as the wind forcing diminishes and the stratification becomes increasingly more convective \citep{Salesky2017}. On the other hand, over mountainous terrain the debate is still ongoing as to the type of coherent structures one can expect, given that even under weak wind forcing, the internally induced thermally-driven circulations dominate the boundary layer flows (Schmidt and Schumann 1989, Babić and De Wekker 2019).

Additional research questions for potential future papers might include (although would require additional output):

5.	Comparison of TKE budget / heat budget terms (at least the resolved ones) across models and different locations on the transect
6.	turbulent spectra (calculated in the homogeneous along-valley/ along-mountain direction) for u, v, w, $\theta$ at different locations across the ridge


## Necessary Variables

To answer the first 5 research questions stated above, the following  variables are essential:

* Mean of u, v, w, θ and pressure (question 1, 2, 3)
    - naming convention: “u_mean”, “v_mean”, “w_mean”, “t_mean”, “p_mean”
* Resolved components of the Reynolds stress tensor (question 1, 2, 3)
    - naming convention: “upup_res”, “upvp_res”, “upwp_res”, “vpvp_res”, “vpwp_res”, “wpwp_res”
* Subgrid components of the Reynolds Stress tensor (question 1, 2, 3)
    - naming convention: “upup_sfs”, “upvp_sfs”, “upwp_sfs”, “vpvp_sfs”, “vpwp_sfs”, “wpwp_sfs”
* $u_*$ and $\theta_*$ at the surface
    - naming convention: “ust”, “tst”
* One 3D snapshot from the last simulation time step for u, v, w, t, p (question 4)
    - naming convention: "u", "v", "w", "t", "p"

In addition, these variables would be desirable to answer the additional research questions:

* θ’θ’,u’θ’, v’θ’, w’θ’, u’p’, v’p’, w’p’, dissipation TKE, dissipation heat, density
    - naming convention: "tptp", "uptp", "vptp", "wptp", "uppp", "vppp", wppp", "eps", "eps_t", "rho"
* u’u’u’, u’u’v’, u’u’w’, v’v’u’, v’v’v’, v’v’w’, w’w’u’, w’w’v’, w’w’w’, u’u’θ’, v’v’θ’, w’w’θ’
    - naming convention: "upupup", "upupvp", "upupwp", "vpvpup", "vpvpvp", "vpvpwp", "wpwpup", "wpwpvp", "wpwpwp", "upuptp", "vpvptp", "wpwptp"
* turbulent spectra for u, v, w, $\theta$

As discussed earlier, turbulent statistics can be calculated either by using recursive averages or by calculating them a posteriori from 1-min instantaneous output of the last 30 minutes of simulation runtime. When calculating statistics on the instantaneous output, it is preferable to have enough grid points in the spanwise direction (→ 1024 grid points spanwise) to guarantee a large enough sample compared to the recursively calculated averages.


## Averaging
As discussed earlier, turbulent statistics can be calculated in two ways:
1. off-line calculation a posteriori: from 1-min instantaneous output of the last 30 minutes of simulation runtime:
    $$
        <\overline{a'b'}>_{ty} = <\overline{(a - <\overline{a}>_{ty})*(b - <\overline{b}>_{ty})}>_{ty}
    $$
    $$
        <\overline{a'b'}>_{ty} (x,z) = \frac{1}{N_t*N_y} \sum_{t=1}^{N_t} \sum_{y=1}^{N_y} [a(x,y,z,t)-<\overline{a}>_{ty}(x,z)]*[b(x,y,z,t)-<\overline{b}>_{ty}(x,z)]
    $$
2. recursive averages: online calculation of averages using the relation
    $$
        <\overline{a'b'}>_{ty} = <\overline{ab}>_{ty} - <\overline{a}>_{ty}*<\overline{b}>_{ty}
    $$
The angled brackets indicate that the averages should be calculated in both time (t) and spanwise (y) direction.

When calculating statistics a posteriori on the instantaneous output (Version 1), it is important to have enough grid points in the spanwise direction (therefore the simulation has to have 1024 grid points in the spanwise direction) to guarantee a large enough sample compared to the recursively calculated averages.

When Version 2 is used (on-line recursive calculation) then less grid points are necessary in the spanwise direction since all time steps are used to create aggregated averages in time and y direction. 


## Desired Data format

The data should be uploaded to a common repository.


Please provide the data in the following format:
* netcdf files with each field as a separate variable, with dimensions (x, z) except for 3D snapshots with dimensions (x, y, z). For each case we ask for the following: 
    - one file containing 1st and 2nd order statistics (plus potentially components for balances → eps, eps_t; or those fr budgets in separate file as well...)
    - one file containing the snapshots for u, v, w, θ and p
    - one file containing the spectra for u, v, w, θ and p
* data interpolated to a common cartesian grid
    - vertical grid: z = [5, 4995] m, Δz = 10 m
    - horizontal grid: x = [0, 10240] m, Δx = 20 m
    - a routine to convert data to a common grid is provided in interpolation_tables.py (function “intdata_2d” for statistics, function “intdata_3d” for snapshots)
* Naming convention: “tables_case_content_model_name.nc”
    - case = [flat_u00, flat_u10, ridge100_u00, ridge100_u10, ridge1000_u00, ridge1000_u10]
    - content = [statistics, spectra, snapshots]
    - e.g. “tables_ridge100_u10_snapshots_cm1_rauchoecker.nc”
    - e.g. “tables_flat_u00_statistics_wrf_serafin.nc”


## References

Babić, N. and S.F.J. De Wekker, 2019: Characteristics of roll and cellular convection in a deep and wide semiarid valley: A large-eddy simulation study, Atmospheric Research, 223, 74–87, doi:10.1016/j.atmosres.2019.03.009.

Leukauf, D., A. Gohm, N.W. Rotach and J.S. Wagner, 2015: The Impact of the Temperature Inversion Breakup on the Exchange of Heat and Mass in an Idealized Valley: Sensitivity to the Radiative Forcing, Journal of Applied Meteorology and Climatology, 54(11), 2199–2216, doi:10.1175/JAMC-D-15-0091.1.

Salesky, S.T., M. Chamecki and E. Bou-Zeid, 2017: On the Nature of the Transition Between Roll and Cellular Organization in the Convective Boundary Layer, Boundary-Layer Meteorology, 163 (1), 41–68, doi:10.1007/s10546-016-0220-3.

Schmidt, H. and U. Schumann 1989: Coherent structure of the convective boundary layer derived from large-eddy simulations, Journal of Fluid Mechanics, 200, 511–562, doi:10.1017/S0022112089000753.

Stiperski, I., M. Calaf and M.W. Rotach, 2019: Scaling, Anisotropy, and Complexity in Near-Surface Atmospheric Turbulence. Journal of Geophysical Research: Atmosphere, 124 (3), 1428–1448, doi:10.1029/2018JD029383.

Stiperski, I., M. Chamecki and M. Calaf 2021: Anisotropy of Unstably Stratified Near-Surface Turbulence. Boundary-Layer Meteorology, 180(3), 363–384, doi:10.1007/s10546-021-00634-0.

Stiperski, I. and M. Calaf, 2023: Generalizing Monin-Obukhov Similarity Theory (1954) for Complex Atmospheric Turbulence. Physical Review Letters, 130 (12), 124001, doi:10.1103/PhysRevLett.130.124001.

Wagner, J.S., A. Gohm and M.W. Rotach 2015: The impact of valley geometry on daytime thermally driven flows and vertical transport processes, Quarterly Journal of the Royal Meteorological Society, 141(690), 1780–1794, doi:10.1002/qj.2481.
