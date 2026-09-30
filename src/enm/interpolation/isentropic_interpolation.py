"""
Interpolate model-level fields onto target isentropic (theta) levels.
"""
import ctypes
import numpy as np
import xarray as xr
from .. import _native

##############################################################
# Wrapper function that gets theta levels and obtains 
# variables on these levels

## Log of checks
# 28/09/26
# -  checked linear interpolation onto theta levels from model levels
# - checked horizontal interpolation of winds
# - checked computation of zeta (gives +/- values in norther/southern hemispheres)
# - interpolation of pressure from model levels to theta


def apply_isentropic_interpolation(ds_int,ds_surf,
                                   exp_type=None,
                                    a = 6371299.,Omega=7.292e-5,g=9.80665): # physical constants

    # get theta levels onto which to interpolate model-level data
    thlev = _get_theta_levels(exp_type=exp_type)
    
    # Compute isentropic vorticity
    ## linearly interpolate zonal wind and potential temperature onto latitude midpoints (slow, rewrite in C++ for speed)
    ds_int['U'] = ds_int['u'] * np.cos(np.radians(ds_int.latitude)) / a
    target_vars = ['U', 'theta']
    lat_mask = ds_int.latitude > -90
    lat_indices = np.where(lat_mask)[0]
    U_th_lat_mids = ds_int[target_vars].map(
        lambda var: .5*(var.roll(latitude=-1, roll_coords=False) + var).isel(latitude=lat_indices)
    )

    ## linearly interpolate meridional wind and potential temperature onto longitude midpoints (slow, rewrite in C++ for speed)
    ds_int['V'] = ds_int['v'] * np.cos(np.radians(ds_int.latitude)) / a
    target_vars = ['V', 'theta']
    V_th_lon_mids = ds_int[target_vars].map(
        lambda var: .5*(var.roll(longitude=-1, roll_coords=False) + var)
    )

    ## linearly interpolate U at latitude midpoints onto theta levels
    U_interp = _interpolate_to_theta_levels(U_th_lat_mids['theta'].values, {'U':U_th_lat_mids['U'].values}, thlev, fill_value=np.nan)
    
    ## linearly interpolate V at longitude midpoints onto theta levels
    V_interp = _interpolate_to_theta_levels(V_th_lon_mids['theta'].values, {'V':V_th_lon_mids['V'].values}, thlev, fill_value=np.nan)
        
    ## compute isentropic vorticity by differening winds
    lats = ds_int.latitude.values
    latmids = np.where(lats<90,(lats + np.roll(lats,shift=1))/2,np.nan)
    latmids = latmids[~np.isnan(latmids)]
    mu = np.sin(np.deg2rad(latmids)) # compute mu = sin(lat) at latitude midpoints
    dUdmu = np.permute_dims(np.diff(np.permute_dims(U_interp['U'],(0,1,3,2)),axis=-1)/np.diff(mu),(0,1,3,2))
    
    lons = ds_int.longitude.values
    dVdl = np.diff(V_interp['V'],axis=-1,append=V_interp['V'][:,:,:,0:1])/np.diff(lons,append=lons[0])
    dVdl = dVdl[...,1:-1,:] # ignore poles since dUdmu is only computed away from poles
    
    mu = np.sin(np.deg2rad(lats))[1:-1] # recompute mu=sin(lat) at gridpoints excluding poles
    zeta = 2*Omega*mu - np.permute_dims(dUdmu,(0,1,3,2)) + (1/(1-mu**2))*np.permute_dims(dVdl,(0,1,3,2))
    zeta = np.permute_dims(zeta,(0,1,3,2)) # (time,theta,lat,lon)
    
    # Compute isentropic density
    ## interpolate pressure onto theta half-levels
    thdiff = np.diff(thlev)
    thmids = thlev[:-1] + thdiff/2
    thmids = np.hstack([thlev[0]-thdiff[0]/2,thmids,thlev[-1]+thdiff[-1]/2])
    p_interp = _interpolate_to_theta_levels(ds_int['theta'].values, {'p':ds_int['p'].values}, thmids, fill_value=np.nan)
    
    ## difference pressure to obtain isentropic density on theta levels
    ## note that total column mass (assuming hydrostatic balance) is conserved by construction of simple finite difference scheme
    r = -(1/g)*np.diff(np.permute_dims(p_interp['p'],(0,2,3,1)),axis=-1)/np.diff(thmids)
    r = np.permute_dims(r,(0,3,1,2))
    r = r[...,1:-1,:] # cut off poles
    
    # Compute isentropic PV
    ## take ratio of vorticity to density and apply lait scaling
    q = zeta/r
    
    # Apply modifications to derived variables
    ## Cap relative vorticity at 2*Omega in the midlatitudes at lower theta levels
    zeta, thlevh, sel, th_sel, lat_sel = _cap_relative_vorticity(zeta,thlev,
                                                                 nmidth=None,nmidlat=None,Omega=Omega,lats=lats)

    ## Modify PV and isentropic density in the midlatitudes at lower theta levels
    qmod, rmod = _modify_r_q(zeta=zeta[sel],q=q[sel],r=r[sel],
                             thlev=thlev[th_sel],thlevh=thlevh[np.append(th_sel,th_sel[-1]+1)],
                             mband=3,lait2pv=np.ones_like(q[sel]))

    q[sel] = qmod
    r[sel] = rmod
    
    # Define output
    vars_dict = {'q' : q, 'zeta' : zeta, 'r' : r}

    return vars_dict

##############################################################
# Funtions for getting theta levels onto which to interpolate

def _get_theta_levels(exp_type=None):
    if exp_type == 'IGCM':
        thlev = _get_igcm_theta_levels()    
    elif exp_type == 'ERA5':
        thlev = _get_era5_theta_levels()
    else:
        thlev = _get_era5_theta_levels()
    return thlev

def _get_igcm_theta_levels(th0=244.,thtop=400.,dth=2.):
    '''Defines evenly spaced theta levels from th0 to thtop spaced by dth'''
    thlev = np.arange(th0,thtop+dth,dth)
    return thlev


def _get_era5_theta_levels(th0=218.,thlow=320.,dthunder=1.5,
                          thmid=450.,thref=380,expansion_fact=0.982,
                          dz=1.,hden=6.5,kappa=2/7,nthlev_strat=43):
    '''
    Defines theta levels for ERA5 experiments
    
    - Overworld levels: evenly spaced in pseudo-height, then converted to
      theta via theta = thref * exp(z / zscal).
    - Underworld levels: start uniform (spacing dthunder), then the spacing
      is smoothly expanded ("bexpand") as theta approaches thmid, so it
      blends into the first Overworld gap without a jump.
    
    Input:
        th0 - minimum theta level
        thlow - theta level upto which levels are evenly spaced
        dthunder - constant spacing for 'underworld' levels
        thmid - determines where to start exponentially-spaced levels
        thref - reference theta determines splicing between underworld and overworld
        expansion_fact - defines expansion of theta levels up to thmid
        hden - density scale height in km (RT/g for T=222K)
        kappa - ratio R/c_p
        nthlev_strat - number of (exponentially spaced) theta levels in the stratosphere
    Output:
        thlev - array of theta levels
    '''
 
    # --- Overworld: evenly spaced in pseudo-height, converted to theta ---
    nthlev = int(nthlev_strat)
    zscal = hden / kappa        # exponential height scale for theta
    z0 = zscal * np.log(thmid / thref)
    zarr = np.arange(nthlev) * dz + dz + z0
    thlev = thref * np.exp(zarr / zscal)
 
    # --- Underworld: uniform, then expanding spacing up to thmid ---
    nunder = round((thlow - th0) / dthunder)
    thunder = list(th0 + np.arange(nunder) * dthunder)
 
    dthmid = thlev[1] - thlev[0]
    bexpand =  expansion_fact * (dthmid - dthunder) / (thmid - thlow)
 
    th = thlow
    while th < thmid:
        thunder.append(th)
        dtcur = dthunder + bexpand * (th - thlow)
        th = th + dtcur
    thunder = np.array(thunder)
 
    # --- Concatenate underworld and overworld ---
    thlev = np.concatenate([thunder, thlev])
    
    return thlev

########################################################################
# Functions for interpolation onto theta levels
#
# Uses the compiled native backend (``enm._native``) when available -- a
# single OpenMP-parallel merge-based pass per column, shared across all
# variables -- and falls back to a pure-NumPy implementation otherwise (e.g.
# when g++/OpenMP wasn't available at install time; see setup.py).
#
#

def _interpolate_to_theta_levels(th, variables, thlev, fill_value=np.nan, use_numpy=False):
    """
    Linearly interpolate variables from model levels onto target theta
    levels.

    Parameters
    ----------
    th : ndarray, shape (ntime, nlev, nlat, nlon)
        Potential temperature on model levels. Must be monotonic along the
        level axis (increasing or decreasing -- direction is
        auto-detected, and assumed consistent across the whole array).
    variables : dict[str, ndarray]
        Other fields to interpolate, each with the same shape as `th`.
        `th` itself is included in the output automatically.
    thlev : ndarray, shape (ntarget,)
        Target theta levels (K), monotonically increasing.
    fill_value : float, default NaN
        Value used where a target level falls outside a column's theta
        range (i.e. would require extrapolation).

    Returns
    -------
    dict[str, ndarray]
        Each field (plus 'th') interpolated onto shape
        (ntime, ntarget, nlat, nlon), float32.
    """
    ntime, nlev, nlat, nlon = th.shape
    ncols = nlat * nlon
    ntarget = thlev.shape[0]

    thlev = np.ascontiguousarray(thlev, dtype=np.float32)
    names = ["th"] + list(variables.keys())
    all_fields = [th] + [variables[name] for name in variables]

    out = {name: np.empty((ntime, ntarget, nlat, nlon), dtype=np.float32)
           for name in names}

    native_lib = _native.load()

    for t in range(ntime):
        th_t = np.ascontiguousarray(th[t].reshape(nlev, ncols), dtype=np.float32)

        # Level-index direction is a global modeling convention (e.g. model
        # level 0 = top of atmosphere vs. = surface), so one check per
        # timestep is enough -- no need to check column by column.
        increasing = th_t[0].mean() < th_t[-1].mean()
        order = slice(None) if increasing else slice(None, None, -1)
        th_sorted = np.ascontiguousarray(th_t[order])

        fields_sorted = np.ascontiguousarray(
            np.stack([np.ascontiguousarray(f[t].reshape(nlev, ncols))[order]
                      for f in all_fields])
        )  # (nvars, nlev, ncols)

        if native_lib is not None and not use_numpy:
            result = _interp_native(
                native_lib, th_sorted, fields_sorted, thlev, fill_value
            )
        else:
            result = _interp_numpy(
                th_sorted, fields_sorted, thlev, fill_value
            )

        for i, name in enumerate(names):
            out[name][t] = result[i].reshape(ntarget, nlat, nlon)

    return out



def _interp_native(lib, th_sorted, fields_sorted, thlev, fill_value):
    """
    Linearly interpolate using the native kernel.

    ctypes hands the C++ side a raw pointer with no dtype attached -- if
    the arrays weren't already contiguous float32 (float64 is numpy's
    default dtype!), the C++ side would silently reinterpret those bytes
    as float32 and produce garbage, with no error at all. So this
    defensively (re)converts every array here, rather than trusting the
    caller to have done it -- unlike _interp_numpy, which tolerates any
    numeric dtype naturally, this is the one place in the module where
    that guarantee has to be made explicit.
    """
    th_sorted = np.ascontiguousarray(th_sorted, dtype=np.float32)
    fields_sorted = np.ascontiguousarray(fields_sorted, dtype=np.float32)
    thlev = np.ascontiguousarray(thlev, dtype=np.float32)

    nlev, ncols = th_sorted.shape
    nvars = fields_sorted.shape[0]
    ntarget = thlev.shape[0]

    out = np.empty((nvars, ntarget, ncols), dtype=np.float32)
    f32p = ctypes.POINTER(ctypes.c_float)

    lib.enm_interp_theta_levels(
        th_sorted.ctypes.data_as(f32p),
        fields_sorted.ctypes.data_as(f32p),
        ctypes.c_long(nlev),
        ctypes.c_long(ncols),
        ctypes.c_long(nvars),
        thlev.ctypes.data_as(f32p),
        ctypes.c_long(ntarget),
        out.ctypes.data_as(f32p),
        ctypes.c_float(fill_value),
    )
    return out

def _interp_numpy(th_sorted, fields_sorted, thlev, fill_value):
    """Pure-NumPy fallback; same semantics as the native kernel, used when
    the native extension hasn't been built."""
    nlev, ncols = th_sorted.shape
    nvars = fields_sorted.shape[0]
    ntarget = thlev.shape[0]

    idx0 = np.empty((ntarget, ncols), dtype=np.int32)
    for k in range(ntarget):
        idx0[k] = np.count_nonzero(th_sorted <= thlev[k], axis=0) - 1
    valid = (idx0 >= 0) & (idx0 < nlev - 1)
    idx0 = np.clip(idx0, 0, nlev - 2)

    th0 = np.take_along_axis(th_sorted[None, :, :], idx0[:, None, :], axis=1)[:, 0, :]
    th1 = np.take_along_axis(th_sorted[None, :, :], (idx0 + 1)[:, None, :], axis=1)[:, 0, :]
    weight = (thlev[:, None] - th0) / (th1 - th0)

    out = np.empty((nvars, ntarget, ncols), dtype=np.float32)
    for v in range(nvars):
        f0 = np.take_along_axis(fields_sorted[v][None, :, :], idx0[:, None, :], axis=1)[:, 0, :]
        f1 = np.take_along_axis(fields_sorted[v][None, :, :], (idx0 + 1)[:, None, :], axis=1)[:, 0, :]
        interp = f0 + weight * (f1 - f0)
        interp[~valid] = fill_value
        out[v] = interp
    return out

########################################################################
# Functions for modifying derived variables (isentropic vorticity, PV and density)

def _cap_relative_vorticity(zeta,thlev,nmidth=None,nmidlat=None,Omega=None,lats=None):
    '''Cap relative vorticity by 2*Omega on isentropic levels upto nmidth
    and at latitude levels upto nmidlat (counting from both poles).
    Parameters
    ----------
    zeta : ndarray
        vorticity with shape (time, theta, lat, lon).
    thlev : ndarray
        Theta-level coordinates.
    nmidth : int, optional
        Maximum theta index to which the adjustment is applied. Defaults to the
        index of the median theta level.
    nmidlat : int, optional
        Maximum absolute latitude index to which the adjustment is applied,
        measured symmetrically about the equator. Defaults to nlat // 4.
    '''
    zeta = np.array(zeta, copy=True)
    _, nth, nlat, _ = zeta.shape

    if nmidth is None:
        nmidth = np.where((thlev>=np.median(thlev)))[0][0]
    else:
        nmidth = min(int(nmidth), nth - 1)

    if nmidlat is None:
        nmidlat = nlat // 4
    else:
        nmidlat = min(int(nmidlat), nlat // 4)

    th_sel = np.arange(0, nmidth + 1)
    lat_sel = np.concatenate([
        np.arange(0, nmidlat + 1),
        np.arange(nlat - nmidlat, nlat),
    ])
    sel = (slice(None), *np.ix_(th_sel, lat_sel), slice(None))

    # compute half isentropic levels
    dth = np.diff(thlev)
    thlevh = np.empty(len(thlev) + 1, dtype=float)
    thlevh[0] = thlev[0] - 0.5 * dth[0]
    thlevh[1:-1] = thlev[:-1] + 0.5 * dth
    thlevh[-1] = thlev[-1] + 0.5 * dth[-1]

    # cap absolute relative vorticity in selected region by 2*Omega
    fcor = 2.0 * Omega * np.sin(np.deg2rad(lats))  # (nlat,)
    fcor_sel = fcor[lat_sel][None, None, :, None]  # shape (1,1,nlat_sel,1)
    fthres = 2.0 * Omega

    zsub = zeta[sel]
    rel = zsub - fcor_sel
    mask_cap = np.isfinite(rel) & (np.abs(rel) > fthres)
    if np.any(mask_cap):
        zsub[mask_cap] = (fcor_sel + fthres * np.sign(rel))[mask_cap]
        zeta[sel] = zsub
        
    # check that maximum absolute relative vorticity does not exceed imposed limit
    rel_vort_ex = np.any(
                        (np.where(np.isnan(zeta[sel]),0,
                                np.abs(zeta[sel]- fcor_sel)) #  zero NaN values in zeta
                                -fthres)/fthres>1e-14)
    
    if rel_vort_ex:
            print('WARNING: 3-d isentropic relative vorticity not properly capped')
        
    return zeta, thlevh, sel, th_sel, lat_sel

def _modify_r_q(zeta, q, r, thlev, thlevh, mband=3, lait2pv=None):
    """
    Modify PV and isentropic density.
    In each column, isentropic density is redistributed to be constant upto first
    local minimum in PV if PV is initially decreasing.

    Parameters
    ----------
    zeta, q, r : ndarray of shape (time, theta, lat, lon)
        Vorticity, PV and isentropic density fields.
    thlev : ndarray, shape (theta,)
        Theta levels.
    thlevh : ndarray, shape (theta+1,)
        Theta half-levels.
    mband : int
        Number of levels to include in the modified lower-density band.
    lait2pv : ndarray, optional
        Lait scaling factor, same shape as q. If None, treated as unity.

    Returns
    -------
    qmod, rmod, minel
        Modified PV, modified density, and the first minimum level index.
    """
    zeta = np.asarray(zeta)
    q = np.asarray(q)
    r = np.asarray(r)

    if lait2pv is None:
        lait2pv = np.ones_like(q)
    else:
        lait2pv = np.asarray(lait2pv, dtype=float)

    ntime, nth, nlat, nlon = q.shape
    minel = np.full((ntime, nlat, nlon), np.nan, dtype=float)

    if nth < 2:  # the loop version never modifies anything in this case
        return q.copy(), r.copy(), minel

    dth = np.diff(thlevh)[:nth]
    dth4 = dth[None, :, None, None]
    k = np.arange(nth)[None, :, None, None]

    # First finite level in each column
    valid = np.isfinite(q)
    anyvalid = valid.any(axis=1)                      # (t, i, j)
    mbot = np.argmax(valid, axis=1)                   # (t, i, j)

    # Need at least two levels above mbot, and q[mbot] > q[mbot+1]
    mb = np.minimum(mbot, nth - 2)[:, None]
    q0 = np.take_along_axis(q, mb, axis=1)[:, 0]
    q1 = np.take_along_axis(q, mb + 1, axis=1)[:, 0]
    with np.errstate(invalid="ignore"):
        eligible = anyvalid & (mbot < nth - 1) & (q0 > q1)

    # First level at/after mbot where q increases; otherwise the top level
    with np.errstate(invalid="ignore"):
        rise = (np.diff(q, axis=1) > 0) & (np.arange(nth - 1)[None, :, None, None]
                                           >= mbot[:, None])
    has_rise = rise.any(axis=1)
    mmin = np.where(has_rise, np.argmax(rise, axis=1), nth - 1)

    minel[eligible] = mmin[eligible]

    do_mod = eligible & (mmin > mbot)

    # Band: [mbot, min(len(thlev), min(mmin + mband, nth-1) + 1))
    stop = np.minimum(len(thlev), np.minimum(mmin + mband, nth - 1) + 1)
    band = (do_mod[:, None]
            & (k >= mbot[:, None])
            & (k < stop[:, None]))

    # Mass-conserving constant density in the band
    mass = np.where(band, r * dth4, 0.0).sum(axis=1)
    thick = np.where(band, dth4, 0.0).sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        rmean = mass / thick

        # round to the input dtype first, as the loop does on assignment
        rmod = np.where(band, rmean[:, None], r).astype(r.dtype, copy=False)
        qmod = np.where(band, zeta / (lait2pv * rmod), q).astype(q.dtype, copy=False)

    return qmod, rmod