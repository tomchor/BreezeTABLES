import os
import numpy as np
import xarray as xr
from scipy.interpolate import griddata

from config_interpolation import (vars_stats, varlist_stats, 
                                  varlist_snapshots, vars_snapshots,
                                  z_true, zold, zfold, x, y, dcoords)

#%% Set output strings and paths
# model and name - for output names
model = "cm1"
name = "rauchoecker"
case = "flat_u00"
content = "statistics"

# filename and path to directory of data to be interpolated
if content=="statistics":
    path = "../out_postproc/"
    filename = "tables_{}_{}_{}.nc".format(case, model, name)
else:
    path = "../../cm1_mod/run/"
    filename = "cm1out_tables_{}.nc".format(case)
outname = "tables_{}_{}_{}_{}.nc".format(case, content, model, name)

print('filename: ', filename)
print('outname: ', outname)

# path to directory where interpolated data should be saved
outpath = "../out_postproc/tables/"
infile = os.path.join(path,filename)
outfile = os.path.join(outpath,outname)

def intdata(ds, dsn, dz=0., zcoord="zh", zfcoord="zf", zinterp="zhval", staggered=False):
    """
    Interpolate data to a common vertical grid.
    This function is needed for data snapshots.

    Parameters:
    ----------
    ds : xr.Dataset
        Simulation output on model-specific vertical grid.
    dsn : xr.Dataset
        Dataset with interpolated vertical grid (empty at beginning).
    dz : float
        height above surface for lowest model level.
        Default is 0.
    zcoord : str
        name of z coordinate in ds.
        Default is 'zh' (in CM1).
    zfcoord : str
        name of z coordinate in ds at flux levels (if applicable).
        Default is 'zf' (in CM1).
    zinterp : str
        name of variable containing true height in ds.
        Default is 'zhval' (in CM1).
    staggered : bool
        If true, assume staggered grid (-> zf needed).
        If false, set zfcoord to "".

    Returns:
    ----------
    dsn : xr.Dataset
        Dataset with interpolated vertical grid (now filled).
    TODO: Make zf optional - code fails if not there!
    """

    if not staggered:
        zfcoord = ""
    
    print("zcoord = " + zcoord)
    print("zfcoord = " + zfcoord)
    print("zinterp = " + zinterp)

    print(dsn.coords)

    for var in list(ds.variables):
        print("interpolating "+var)
        print(ds[var].coords)
        if (((zcoord in ds[var].coords) & (var!=zcoord) & (var!=zfcoord)) |
            ((zfcoord in ds[var].coords) & (var!=zfcoord) & (var!=zcoord))):
            if (("x" in ds[var].coords) & ("y" in ds[var].coords)):
                # for data snapshots (3D)
                data = np.zeros((len(dsn.y), len(dsn.z), len(dsn.x)))
                for ix, x in enumerate(ds.x):
                    _, zpoints = np.meshgrid(dsn.y, dsn.z)
                    if zcoord==zinterp:
                        zc = zcoord
                        intheight = ds[zcoord]
                    else:
                        if zcoord in ds[var].coords:
                            zc = zcoord
                            intheight = ds[zinterp].isel(x=ix, y=-1).data
                        elif zfcoord in ds[var].coords:
                            zc = zfcoord
                            intheight = (ds[zinterp].isel(x=ix, y=-1)-dz).data
                    data[:,:,ix] = xr.apply_ufunc(griddata,
                                                  intheight,
                                                  ds[var].sel(x=x).T,
                                                  zpoints.T,
                                                  input_core_dims=[[zc], [zc], ["z"]],
                                                  output_core_dims=[["z"]],
                                                  exclude_dims=set((zc,)),
                                                  vectorize=True,
                                                  dask="parallelized")
                    dsn[var] = xr.DataArray(data=data, coords=dsn.coords)
            elif (("x" in ds[var].coords) & ("y" not in ds[var].coords)):
                # for data statistics of ridge cases (2D)
                data = np.zeros((len(dsn.z), len(dsn.x)))
                for ix, x in enumerate(ds.x):
                    if zcoord==zinterp:
                        intheight = ds[zcoord]
                    else:
                        if zcoord in ds[var].coords:
                            intheight = ds[zinterp].isel(x=ix).data
                        elif zfcoord in ds[var].coords:
                            intheight = (ds[zinterp].isel(x=ix)-dz).data
                    data[:,ix] = griddata(intheight, ds[var].sel(x=x), dsn.z, method="linear")
                dsn[var] = xr.DataArray(data=data, coords=dsn.coords)
            elif (("x" not in ds[var].coords) & ("y" not in ds[var].coords)):
                # for data statistics of flat cases (1D)
                data = np.zeros((len(dsn.z)))
                if zcoord==zinterp:
                    intheight = ds[zcoord]
                else:
                    if zcoord in ds[var].coords:
                        intheight = ds[zinterp].data
                    elif zfcoord in ds[var].coords:
                        intheight = (ds[zinterp]-dz).data
                data = griddata(intheight, ds[var], dsn.z, method="linear")
                dsn[var] = xr.DataArray(data=data, coords=dsn.coords)
        else:
            if ((var==zfcoord) | (var==zcoord)):
                pass
            else:
                dsn[var] = ds[var]
    return dsn 



if __name__=="__main__":
    # choose case: snapshots or statistics
    if content=="snapshots":
        varlist = varlist_snapshots
        varsd = vars_snapshots
    elif content=="statistics":
        varlist = varlist_stats
        varsd = vars_stats
    else:
        raise NameError("Unknown output.")

    # select base coordinate for vertical grid
    if "ridge" in filename:
        # in CM1, the variable "zhval" gives absolute height.
        # It is only part of the output for simulations with terrain.
        zinterp = z_true
        if "statistics" in outname:
            varsd.update({x: "x"})
        varlist.append(z_true)
    else:
        # use vertical coordinate for flat simulations instead.
        zinterp = zold
    
    # new vertical grid
    znew = np.arange(5.+1e-6, 3000., 10.)

    print("interpolate...")

    print("open dataset")
    ds = xr.open_dataset(infile)
    
    # select last time step if multiple time steps in dataset 
    if "time" in ds.coords:
        if len(ds.time.shape)>0:
            ds = ds.isel(time=-1)
        ds = ds.drop("time")

    for dc in dcoords:
        if dc in list(ds.variables):
            # if coordinates are given in km instead of m, change that here
            if ds[dc].max()<=1000.:
                ds[dc] = ds[dc]*1000.
                ds[dc].attrs["units"] = "m"
            print(ds[dc].max())
            # if coordinates are given in km instead of m, change that here
            if ((ds[dc].min()<0) & ((dc==x) | (dc==y))):
                ds[dc] = ds[dc]-ds[dc].min()

            print(ds[dc].min(), ds[dc].max())

    # select variables that are needed, discard the rest
    ds = ds[varlist]

    # rename variables selected earlier
    print("rename variables")
    ds = ds.rename(varsd)

    # last sanity check: statistics should not contain y coordinate
    if content=="statistics":
        for var in list(ds.variables):
            if ("y" in ds[var].coords) & (var != "y"):
                print("Calculate spanwise average of " + var)
                ds[var] = ds[var].mean(dim="y")
        ds = ds.drop_vars("y")


    print(ds.coords)

    if zfold in ds.coords:
        # accounts for grid staggering if staggered data are not already destaggered
        # need to replace zf in this line with appropriate name in your dataset
        ds = ds.sel(zf=ds[zfold][:-1])
        dz = (ds[zold].data-ds[zfold].data)
        staggered = True
        # the next line obviously leads to some problems above the ridge when grid is stretched
        # therefore better to destagger already beforehand
    else:
        dz = 0.
        staggered = False
    print("interpolate variables")
    if  (("x" in ds.coords) & ("y" not in ds.coords)):
        # for ridge data
        coords = dict(z=(["z"], znew),
                      x=(["x"], ds["x"].data))
        dsn = xr.Dataset(coords=coords)
        dsn = intdata(ds, dsn, dz=dz, zcoord=zold, zfcoord=zfold, zinterp=zinterp, staggered=staggered)
    elif (("x" not in ds.coords) & ("y" not in ds.coords)):
        # for flat data
        coords = dict(z=(["z"], znew))
        dsn = xr.Dataset(coords=coords)
        dsn = intdata(ds, dsn, dz=dz, zcoord=zold, zfcoord=zfold, zinterp=zinterp, staggered=staggered)
    elif (("x" in ds.coords) & ("y" in ds.coords)):
        # for snapshots
        coords = dict(y=(["y"], ds["y"].data),
                      z=(["z"], znew),
                      x=(["x"], ds["x"].data))
        dsn = xr.Dataset(coords=coords)
        dsn = intdata(ds, dsn, dz=dz, zcoord=zold, zfcoord=zfold, zinterp=zinterp, staggered=staggered)
    else:
        print("Unknown dimensions {}...".format(list(ds.coords)))


    print("saving renamed dataset {}...".format(outname))
    dsn.to_netcdf(outfile)
