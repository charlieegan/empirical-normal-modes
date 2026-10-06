"""
Interpolate model-level fields onto target isentropic (theta) levels and derive isentropic
density (r), vorticity (zeta) and Ertel PV (q).
"""
import ctypes
import numpy as np
import xarray as xr
from .. import _native

def apply_isentropic_interpolation(ds_int,ds_surf,
                                        exp_type=None,
                                        a = 6371299.,Omega=7.292e-5,g=9.80665,kappa=2./7.,
                                        rdgas=287.,p00=1.e5,zs=0., # extra physical constants, only
                                                                   # needed for the Montgomery potential
                                        nbound=1,mband=2,mumod=0.76,thref=380.,nlat_active=None):
    '''
    Interpolate model-level data onto isentropic levels and compute isentropic vorticity (zeta),
    isentropic density (r) and modified (Lait) PV (q). Translated directly from the IDL code
    backwa_pvinv_nbs.pro using Claude Sonnet 5.

    Parameters
    ----------
    ds_int : as apply_isentropic_interpolation, but must also carry 'pv' (Ertel PV, SI units -- the
             IDL reads it in PVU and divides by 1e6; legacy_tradv.read_tradv's own 'pv' is still in
             PVU, so the caller must do that conversion before passing it in here).
    ds_surf : must carry 'lnsp' (as produced by legacy_tradv.to_current_repo_format) -- unlike
              apply_isentropic_interpolation, this function actually uses it, for surface pressure.
    rdgas : dry air gas constant (IDL: rdgas, 287.).
    p00 : reference pressure for potential temperature (IDL: p00, 1e5 Pa).
    zs : surface geopotential (IDL: zs; g*height, m^2 s^-2), broadcastable to (time,lat,lon). Defaults
         to 0 (flat/no orography) -- legacy_tradv test data has none; substitute real surface
         geopotential once available (e.g. from ERA5) for a faithful comparison.
    See apply_isentropic_interpolation for the remaining parameters.

    Output q, zeta, r and the boundary fields (thetalb, plb, ptop, mbot) are all full nlat in the
    latitude dimension, matching the IDL's array shapes exactly: q/zeta have NaN in the southernmost
    row, since the IDL's own vorticity loop never computes it either (see the note above dUdmu).

    Known limitation shared with the IDL: a column whose lower-boundary theta sits at or above the
    second-highest theta level (mbot >= nth-1) will index out of bounds in the bottom-density formula
    below -- the IDL has exactly the same unguarded assumption (mbot+1 implicitly assumed < nthlev).
    '''
    cp = rdgas/kappa
    expo = (1.-kappa)/kappa
    gasfac = p00/(rdgas*g)

    # get theta levels onto which to interpolate model-level data
    thlev = _get_theta_levels(exp_type=exp_type)
    nth = len(thlev)
    thlevh = _get_half_levels(thlev)
    lait2pv = _get_lait2pv(thlev,thref=thref)
    lats = ds_int.latitude.values

    # Lower boundary of the isentropic domain. Must use the unmodified model-level data.
    theta_ml = ds_int['theta'].values
    p_ml = ds_int['p'].values
    nlev = theta_ml.shape[1]

    # Top of the isentropic domain reliably covered by the data (IDL: topminth, nthlim).
    nthlim, topminth = _get_top_theta_limit(theta_ml,thlevh)

    thetalb, plb = _get_lower_boundary(theta_ml,p_ml,nbound=nbound)
    mbot = _get_lowest_theta_index(thlev,thetalb)
    medianths, mmed = _get_mmed(thetalb,lats,thlev,mumod=mumod)

    # Compute isentropic vorticity -- identical to apply_isentropic_interpolation.
    ds_int['U'] = ds_int['u'] * np.cos(np.radians(ds_int.latitude)) / a
    target_vars = ['U', 'theta']
    # Midpoint between the last and first row (wrapping around a pole) is meaningless and must be
    # dropped; this is a wraparound exclusion by position, not a check for an exact +/-90 gridpoint
    # (a Gaussian grid, e.g. from legacy_tradv, never has one -- only a regular lat/lon grid does).
    lat_indices = np.arange(len(lats) - 1)
    U_th_lat_mids = ds_int[target_vars].map(
        lambda var: .5*(var.roll(latitude=-1, roll_coords=False) + var).isel(latitude=lat_indices)
    )
    ds_int['V'] = ds_int['v'] * np.cos(np.radians(ds_int.latitude)) / a
    target_vars = ['V', 'theta']
    V_th_lon_mids = ds_int[target_vars].map(
        lambda var: .5*(var.roll(longitude=-1, roll_coords=False) + var)
    )
    U_interp = _interpolate_stable(U_th_lat_mids['theta'].values, {'U':U_th_lat_mids['U'].values}, thlev, fill_value=np.nan)
    V_interp = _interpolate_stable(V_th_lon_mids['theta'].values, {'V':V_th_lon_mids['V'].values}, thlev, fill_value=np.nan)

    ## compute isentropic vorticity by differencing winds.
    ## IDL stores u/v on a latitude half-grid with one extra row at each pole (muh(0)=1 at the north
    ## pole, muh(nlat)=-1/0 at the south -- see backwa_pvinv_nbs.pro), which never gets assigned a
    ## value, so those rows are implicitly U=V=0. That lets it difference across the row next to the
    ## north pole using this zero boundary, but its main loop only ever runs j=0..nlat-2
    ## (`for j=0,nlatend-1`, nlatend=nlat-1), so the southernmost row is simply never computed. This
    ## asymmetry (pole boundary condition in the north, row just dropped in the south) is reproduced
    ## here as-is for comparability -- it's a property of the legacy code, not of whether the real
    ## data reaches the poles. mu=1 (not lats[0]) is the fixed boundary value, matching IDL's
    ## convention regardless of where the first real data row actually is (e.g. ~89.46 deg on a
    ## Gaussian grid, vs exactly 90 on a regular lat/lon grid).
    mu_full = np.sin(np.deg2rad(lats))  # mu(j) = sin(lat(j)), all nlat rows
    mu_mids = 0.5*(mu_full[:-1] + mu_full[1:])  # muh(j) for j=1..nlat-1: nlat-1 midpoints
    mu_mids = np.concatenate([[1.], mu_mids])  # prepend muh(0)=1, the fixed north-pole boundary
    U_pole = np.zeros_like(U_interp['U'][:,:,:1,:])  # explicit U=0 boundary, for every theta level
    U_mids_padded = np.concatenate([U_pole, U_interp['U']], axis=2)  # (time,theta,nlat,lon)
    dUdmu = np.diff(U_mids_padded,axis=2) / np.diff(mu_mids)[None,None,:,None]  # (time,theta,nlat-1,lon)

    lons = np.deg2rad(ds_int.longitude.values)
    dVdl = np.diff(V_interp['V'], axis=-1, prepend=V_interp['V'][:,:,:,-1:]) / np.diff(lons, prepend=lons[-1:]-2*np.pi)
    dVdl = dVdl[...,:-1,:] # match dUdmu's row range (0..nlat-2); the last row is never computed (see above)

    mu = np.sin(np.deg2rad(lats))[:-1][None,None,:,None] # sin(lat) at the computed rows, 0..nlat-2
    utermarr = - dUdmu
    vtermarr = (1/(1-mu**2))*dVdl
    zeta = 2*Omega*mu - dUdmu + (1/(1-mu**2))*dVdl # (time,theta,nlat-1,lon)
    nan_row = np.full_like(zeta[:,:,:1,:], np.nan)
    utermarr = np.concatenate([utermarr, nan_row], axis=2)
    vtermarr = np.concatenate([vtermarr, nan_row], axis=2)
    zeta = np.concatenate([zeta, nan_row], axis=2) # (time,theta,nlat,lon); last row never computed

    # Isentropic density via the Montgomery potential (IDL: zarr -> marr/mont -> denth).
    ps = np.exp(ds_surf['lnsp'].values[:,0]) # (time,lat,lon)
    t_ml = theta_ml * (p_ml/p00)**kappa # temperature on model levels (IDL: tarr, first-loop usage)

    ## Hydrostatic integration of geopotential from the surface upward (IDL: zarr), one model half-
    ## layer at a time, exactly as the IDL does it.
    z_ml = np.empty_like(p_ml)
    z_ml[:,-1] = zs + rdgas*t_ml[:,-1]*np.log(ps/p_ml[:,-1])
    for l in range(nlev-2,-1,-1):
        z_ml[:,l] = z_ml[:,l+1] + rdgas*0.5*(t_ml[:,l+1]+t_ml[:,l])*np.log(p_ml[:,l+1]/p_ml[:,l])
    zlb = z_ml[:,-nbound] # geopotential at the lower-boundary model level (IDL: zlb, lbound=1 case)

    ## Pressure (-> temperature) on full theta levels, and pressure+geopotential on the single top
    ## half-level boundary thtop=thlevh[-1] (IDL: tth/tarr and ptop/ztop respectively). The
    ## fallback-to-top-model-level behaviour for thtop mirrors apply_isentropic_interpolation's ptop.
    full_interp = _interpolate_stable(theta_ml, {'p':p_ml}, thlev, fill_value=np.nan)
    p_full = full_interp['p']
    p_full = np.where(np.isfinite(p_full),p_full,p_ml[:,0:1])
    T_full = thlev[None,:,None,None] * (p_full/p00)**kappa # IDL: tth / tarr (second-loop usage)

    half_interp = _interpolate_stable(theta_ml, {'p':p_ml,'z':z_ml}, thlevh, fill_value=np.nan)
    ptop = half_interp['p'][:,-1,:,:]
    ptop = np.where(np.isfinite(ptop),ptop,p_ml[:,0])
    ztop = half_interp['z'][:,-1,:,:]
    ztop = np.where(np.isfinite(ztop),ztop,z_ml[:,0])
    ttop = thlevh[-1] * (ptop/p00)**kappa # IDL: ttop

    ## Integrate M downward from the top (IDL: marr). Computed for every level/column uniformly;
    ## only levels at/above each column's own mbot are actually used below.
    mthtop = cp*T_full[:,-1] + ztop # IDL: mthtop
    dtharr = np.concatenate([np.diff(thlev), [thlevh[-1]-thlev[-1]]]) # full-level spacing (IDL: dtharr)
    dthothh = dtharr / thlevh[1:] # spacing normalised by the intervening half level (IDL: dthothh)
    M = np.empty_like(T_full)
    M[:,-1] = mthtop
    for m in range(nth-2,-1,-1):
        M[:,m] = M[:,m+1] - cp*0.5*(T_full[:,m+1]+T_full[:,m])*dthothh[m]

    ## Differentiate M to get density: a quadratic fit through (thlb,thlev[mbot],thlev[mbot+1]) at the
    ## lower boundary, a standard 3-point finite difference for interior levels, and a one-sided
    ## formula anchored on the top boundary -- the IDL's three separate regimes (lines 2132-2166).
    thm_i,tho_i,thp_i = (x[None,:,None,None] for x in (thlev[:-2],thlev[1:-1],thlev[2:]))
    Mm,Mo,Mp = M[:,:-2], M[:,1:-1], M[:,2:]
    dmdth_int = (Mp-Mm)/(thp_i-thm_i)
    d2mdth2_int = ((Mp-Mo)/(thp_i-tho_i) - (Mo-Mm)/(tho_i-thm_i)) * 2./(thp_i-thm_i)
    r_int = -gasfac*d2mdth2_int*(dmdth_int/cp)**expo # (time,nth-2,lat,lon), for full levels 1..nth-2

    backward_deriv = (M[:,-1]-M[:,-2])/(thlev[-1]-thlev[-2])
    dmdth_top = cp*T_full[:,-1]/thlev[-1]
    d2mdth2_top = (cp*ttop/thlevh[-1] - backward_deriv) / (thlevh[-1]-thlevh[-2])
    r_top = -gasfac*d2mdth2_top*(dmdth_top/cp)**expo # (time,lat,lon), top full level only

    thm_b,thp_b = thlev[mbot], thlev[mbot+1] # (time,lat,lon); see the docstring's known-limitation note
    mm_b = np.take_along_axis(M,mbot[:,None],axis=1)[:,0]
    mp_b = np.take_along_axis(M,(mbot+1)[:,None],axis=1)[:,0]
    dth_b = thp_b - thm_b
    under = dth_b*(thp_b*thm_b - thetalb*thetalb)
    dmdth_bot = (mp_b*(thm_b**2-thetalb**2) + mm_b*(thetalb**2+thp_b**2-2*thp_b*thm_b)
                 - zlb*dth_b**2) / under
    d2mdth2_bot = (mp_b*thm_b - mm_b*thp_b + zlb*dth_b)*2./under
    r_bot = -gasfac*d2mdth2_bot*(dmdth_bot/cp)**expo # (time,lat,lon), mbot level only

    r = np.full(T_full.shape,np.nan)
    r[:,1:-1] = r_int
    r[:,-1] = r_top
    np.put_along_axis(r,mbot[:,None],r_bot[:,None],axis=1) # overwrite the per-column mbot level
    # r keeps the full nlat -- no latitude-differencing dependency, so (unlike zeta) there's no row
    # it can't compute; matches denth's full-row range, and thetalb/plb/ptop/mbot below are already
    # full nlat too (no latitude-differencing dependency, computed directly per-column).

    # Distribute mass uniformly on the lowest band of theta levels, mbot..mbot+mband, so that the mass of
    # each column between the lower boundary pressure and ptop is respected exactly (IDL: lbound=1).
    all_cols = np.ones(mbot.shape,dtype=bool)
    r = _set_bottom_density(r,mbot,mbot+mband,thlevh,thetalb,plb,ptop,all_cols,g=g)

    # Levels below the lower boundary are not part of the domain
    below = np.arange(nth)[None,:,None,None] < mbot[:,None]
    r = np.where(below,np.nan,r).astype(r.dtype,copy=False)
    zeta = np.where(below,np.nan,zeta)
    zeta_unmod = zeta.copy()

    # Cap relative vorticity and flatten density up to the first PV minimum in the polar lower
    # troposphere, then compute modified (Lait) PV
    zeta, r, q, modified = _modify_r_q(zeta,r,lait2pv,thlevh,thetalb,plb,ptop,mbot,mmed,
                                       lats,Omega=Omega,g=g,mband=mband,nlat_active=nlat_active)

    # Define output
    vars_dict = {'U_lat_mids' : U_interp['U'], 'V_lon_mids' : V_interp['V'], 'p' : p_full, # interpolated variables
                 'utermarr' : utermarr, 'vtermarr' : vtermarr, 'zeta_unmod' : zeta_unmod, # contributions to vorticity
                 'mont' : M, 'zlb' : zlb, 'ztop' : ztop, # new vs. apply_isentropic_interpolation
                 'q' : q, 'zeta' : zeta, 'r' : r, # derived variables
                 'thetalb' : thetalb, 'plb' : plb, 'ptop' : ptop, 'mbot' : mbot, # boundary diagnostics
                 'medianths' : medianths, 'mmed' : mmed, 'lait2pv' : lait2pv, 'thlev' : thlev, 'thlevh' : thlevh, # levels
                 'polar_modified' : modified, 'nthlim' : nthlim, 'topminth' : topminth}

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

def _get_half_levels(thlev):
    '''Theta half-levels (IDL thlevh): midpoints between levels, extended by half a spacing at each end.
    Length is len(thlev)+1; thlevh[-1] is the top boundary thtop.'''
    dth = np.diff(thlev)
    thlevh = np.empty(len(thlev) + 1, dtype=float)
    thlevh[0] = thlev[0] - 0.5 * dth[0]
    thlevh[1:-1] = thlev[:-1] + 0.5 * dth
    thlevh[-1] = thlev[-1] + 0.5 * dth[-1]
    return thlevh

def _get_lait2pv(thlev,thref=380.):
    '''Factor converting Ertel PV to modified (Lait) PV: (theta/thref)^4.5 for theta >= thref, else 1'''
    lait2pv = np.ones(len(thlev))
    above = thlev >= thref
    lait2pv[above] = (thlev[above] / thref) ** 4.5
    return lait2pv

########################################################################
# Function for the top of the isentropic domain covered by the data

def _get_top_theta_limit(theta,thlevh):
    '''Index of the highest theta level reliably covered by the data at the top of the model domain
    (IDL: nthlim, topminth).

    topminth is the minimum potential temperature on the top model level over the whole horizontal
    domain: the highest theta guaranteed to be covered by every column at the top of the data. nthlim
    is the index of the last theta half-level at or below topminth, i.e. the count of theta levels fully
    covered by the data (IDL finds the first half-level above topminth and subtracts one, falling back
    to the top theta level, nthlev, if none are above -- both cases reduce to the same formula here,
    since topminth is always well above thlevh[0] in practice). This is a diagnostic boundary consumed
    by later wave-activity stages (which loop over levels 0:nthlim-1); it does not affect q/zeta/r here.

    theta : (time, model_level, lat, lon), level 0 = top. thlevh : (ntheta+1,) half levels.
    Returns nthlim : (time,) int, topminth : (time,)
    '''
    topminth = theta[:,0].min(axis=(-2,-1))
    nthlim = np.searchsorted(thlevh,topminth,side='right') - 1
    return nthlim, topminth

########################################################################
# Functions for the lower boundary of the isentropic domain

def _get_lower_boundary(theta,p,nbound=1):
    '''Theta and pressure of the lower boundary of the isentropic domain (IDL, lbound=1).
    The boundary is model level nlev-nbound. If the lowest model level is warmer than that level
    (statically unstable layer between them), the boundary theta is raised to the lowest-level value,
    i.e. the unstable layer is moved below the boundary. With nbound=1 this never changes anything.
    
    theta, p : (time, model_level, lat, lon), level 0 = top, last level = lowest.
    Returns thetalb, plb : (time, lat, lon)
    '''
    th_surf = theta[:,-1]
    th_lb = theta[:,-nbound]
    thetalb = np.where(th_surf > th_lb, th_surf, th_lb)
    plb = p[:,-nbound]
    return thetalb, plb

def _get_lowest_theta_index(thlev,thetalb):
    '''Index of the lowest theta level at or above the lower-boundary theta (IDL: mbot=where(thlev ge thsurf)(0)).
    Equals len(thlev) where the boundary theta exceeds all theta levels.'''
    return np.searchsorted(thlev,thetalb,side='left')

def _get_mmed(thetalb,lats,thlev,mumod=0.76):
    '''Define the theta range of the polar lower troposphere (IDL: medianths, mmed).
    medianths : zonal-mean lower-boundary theta on the first latitude row (from the north) with
                sin(lat) < mumod.
    mmed : index of the first theta level above medianths. 
    thetalb : (time, lat, lon) on the full latitude grid; lats : full latitude vector (north to south).
    '''
    mu = np.sin(np.deg2rad(lats))
    below = np.where(mu < mumod)[0]
    ntime = thetalb.shape[0]
    if below.size == 0:
        return np.full(ntime,np.nan), np.zeros(ntime,dtype=int)
    medianths = thetalb[:,below[0],:].mean(axis=-1)
    mmed = np.searchsorted(thlev,medianths,side='right')
    return medianths, mmed

def _set_bottom_density(r,mbot,minband,thlevh,thetalb,plb,ptop,cols,g=9.80665):
    '''Overwrite isentropic density on levels mbot..minband (inclusive) with a constant such that the mass
    of the column between plb and ptop is respected exactly.
    IDL: massbot=(plb-ptop)*p00/ga-denint ; denth(mbot:minband)=massbot/(thlevh(minband+1)-thlb)
    where denint integrates density over levels above minband up to the top boundary.
    
    r : (time, theta, lat, lon); mbot, minband, thetalb, plb, ptop : (time, lat, lon) 
    cols : (time, lat, lon) bool, columns to modify.
    Non-finite density above the band (levels above the data) is treated as zero mass.
    '''
    nth = r.shape[1]
    dthh = np.diff(thlevh)[None,:,None,None]
    k = np.arange(nth)[None,:,None,None]
    ok = cols & (minband <= nth-2)
    above = k > minband[:,None]
    term = np.where(above & np.isfinite(r), r*dthh, 0.0)
    denint = term.sum(axis=1,dtype=np.float64)
    massbot = (plb.astype(np.float64) - ptop.astype(np.float64))/g - denint
    top_edge = thlevh[np.clip(minband+1,0,nth)]
    rbot = massbot/(top_edge - thetalb)
    band = ok[:,None] & (k >= mbot[:,None]) & (k <= minband[:,None])
    return np.where(band,rbot[:,None],r).astype(r.dtype,copy=False)

########################################################################
# Functions for interpolation onto theta levels
#
# Requires the compiled native backend (``enm._native``) -- a single
# OpenMP-parallel merge-based pass per column, shared across all variables.
# There is no NumPy fallback; see isentropic_interp.cpp for the kernel.
#

def _interpolate_stable(th, variables, thlev, fill_value=np.nan):
    """
    Interpolate to theta levels in columns that may be statically unstable (theta not monotonic
    with height), e.g. an inversion of any depth near the surface.

    The native kernel (enm._native) brackets each target theta level by scanning from the surface
    upward and taking the LAST level with theta below the target -- equivalent to the IDL's
    `where(theta lt thm)` taking the first hit scanning down from the top -- so it does not require
    the column to be monotonic (see isentropic_interp.cpp). No pre-processing of theta is needed here;
    this wrapper only adds the lower-boundary masking below, which the interpolation kernel itself
    doesn't know about (theta below the column's own lowest-level value is a valid bracket, but theta
    below the *lower-boundary* value, th_bot, is outside the isentropic domain -- see mbot in
    apply_isentropic_interpolation).

    th : (time, model_level, lat, lon), level 0 = top, last level = lowest.
    """
    th_bot = th[:,-1]
    out = _interpolate_to_theta_levels(th,variables,thlev,fill_value=fill_value)
    below = np.asarray(thlev)[None,:,None,None] < th_bot[:,None]
    for name in out:
        out[name][np.broadcast_to(below,out[name].shape)] = fill_value
    return out

def _interpolate_to_theta_levels(th, variables, thlev, fill_value=np.nan):
    """
    Linearly interpolate variables from model levels onto target theta
    levels, via the native kernel (enm._native.enm_interp_theta_levels).

    Parameters
    ----------
    th : ndarray, shape (ntime, nlev, nlat, nlon)
        Potential temperature on model levels. Columns need not be
        monotonic -- see isentropic_interp.cpp.
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
    if native_lib is None:
        raise RuntimeError(
            "enm._native extension (libenm_native.so) is not built; the isentropic interpolation "
            "requires it (no NumPy fallback). Rebuild via setup.py (needs g++/OpenMP)."
        )

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

        result = _interp_native(
            native_lib, th_sorted, fields_sorted, thlev, fill_value
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
    caller to have done it.
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

########################################################################
# Functions for modifying derived variables (isentropic vorticity, PV and density)
#
# IDL logic (applied per column, global data only, where mbot+1 < mmed):
#   1. for theta levels mbot..mmed cap |zeta - f| at 2*Omega 
#   2. search upward from mbot for the first local minimum of (Lait) PV, level minel
#   3. if minel > mbot (PV decreasing away from the boundary), overwrite density on mbot..minel+mband with 
#      a constant that conserves column mass, and recompute PV = zeta/(lait2pv*r)

def _cap_relative_vorticity(zeta,mbot,mmed,fcor,fthres):
    '''Cap relative vorticity |zeta - fcor| at fthres on theta levels mbot..mmed (inclusive).
    zeta : (time, theta, lat, lon);  mbot, mmed : broadcastable to (time, lat, lon);
    fcor : Coriolis parameter 2*Omega*sin(lat), broadcastable to zeta; fthres : cap (IDL: 2*Omega).
    Only pass columns that are in the modification region (see _get_polar_columns).'''
    nth = zeta.shape[1]
    k = np.arange(nth)[None,:,None,None]
    region = (k >= mbot[:,None]) & (k <= mmed[:,None])
    rel = zeta - fcor
    with np.errstate(invalid='ignore'):
        cap = region & (np.abs(rel) > fthres)
    return np.where(cap,fcor + fthres*np.sign(rel),zeta).astype(zeta.dtype,copy=False)

def _first_minimum_index(q,mbot):
    '''Index of the first local minimum of q searching upward (increasing theta index) from mbot.
    IDL: minel=m-1 where m is the first level >= mbot+1 with q(m) > q(m-1) (a NaN also stops the search).
    If q never rises, the top level is returned. Equals mbot when q(mbot+1) > q(mbot).'''
    nth = q.shape[1]
    with np.errstate(invalid='ignore'):
        rise = ~(q[:,1:] <= q[:,:-1])          # rise[:, i] refers to level m = i+1
    m = np.arange(1,nth)[None,:,None,None]
    rise &= m >= (mbot[:,None] + 1)
    has_rise = rise.any(axis=1)
    first = np.argmax(rise,axis=1) + 1
    return np.where(has_rise,first - 1,nth - 1)

def _get_polar_columns(mbot,mmed,lats,nlat_active=None,global_data=None):
    '''Columns in which the IDL applies the polar lower-troposphere modification: global data only, and 
    mbot+1 < mmed. lats must be the latitudes of the rows of mbot. nlat_active restricts to the first rows.'''
    if global_data is None:
        global_data = bool(np.min(lats) < 0)
    cols = (mbot + 1 < mmed[:,None,None]) & global_data
    if nlat_active is not None:
        cols &= (np.arange(len(lats)) < nlat_active)[None,:,None]
    return cols

def _modify_r_q(zeta,r,lait2pv,thlevh,thetalb,plb,ptop,mbot,mmed,lats,
                Omega=7.292e-5,g=9.80665,mband=2,nlat_active=None,global_data=None):
    """
    Modify vorticity and isentropic density in the polar lower troposphere, and compute PV.
    
    In columns where the IDL modifies the state (see above), vorticity is capped and density is made
    constant (mass-conserving) from the lower boundary up to mband levels above the first PV minimum.
    Everywhere, q = zeta/(lait2pv*r) (modified/Lait PV, SI units).
    
    Parameters
    ----------
    zeta, r : (time, theta, lat, lon); NaN below the lower boundary. 
    lait2pv : (theta,) ; thlevh : (theta+1,) ; 
    thetalb, plb, ptop, mbot : (time, lat, lon) ; mmed : (time,) ; lats : (lat,) latitudes of the rows
    
    Returns
    -------
    zeta, r, q : modified fields, and `modified` : (time, lat, lon) bool, columns where density was changed.
    """
    zeta = np.array(zeta,copy=True)
    r = np.array(r,copy=True)
    lait = np.asarray(lait2pv)[None,:,None,None]
    fthres = 2.0*Omega
    active = _get_polar_columns(mbot,mmed,lats,nlat_active=nlat_active,global_data=global_data)
    modified = np.zeros(active.shape,dtype=bool)
    
    t_i, j_i, i_i = np.nonzero(active)
    if t_i.size > 0:
        # Gather the active columns into pseudo-4D arrays (1, theta, ncol, 1) and reuse the 4-D routines
        z = zeta[t_i,:,j_i,i_i].T[None,:,:,None]
        rr = r[t_i,:,j_i,i_i].T[None,:,:,None]
        col = lambda x: x[t_i,j_i,i_i][None,:,None]
        mb = col(mbot)
        fcor = (2.0*Omega*np.sin(np.deg2rad(lats)))[j_i][None,None,:,None]
        z = _cap_relative_vorticity(z,mb,mmed[t_i][None,:,None],fcor,fthres)
        with np.errstate(divide='ignore',invalid='ignore'):
            qc = z/(lait*rr)
        minel = _first_minimum_index(qc,mb)
        mod = minel > mb
        rr = _set_bottom_density(rr,mb,minel+mband,thlevh,col(thetalb),col(plb),col(ptop),mod,g=g)
        zeta[t_i,:,j_i,i_i] = z[0,:,:,0].T
        r[t_i,:,j_i,i_i] = rr[0,:,:,0].T
        modified[t_i,j_i,i_i] = mod[0,:,0]
    
    with np.errstate(divide='ignore',invalid='ignore'):
        q = zeta/(lait*r)
    return zeta, r, q, modified