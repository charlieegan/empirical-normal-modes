'''
Code to compute mass and circulation integrals within PV contours on isentropic levels.
This, along with relevant metadata, is used as input for the background state calculation.
'''

import numpy as np
import xarray as xr


def get_bgs_input(ds_der,thref=380.):
    '''Function to get input for background state computation from 3-d derived variables and surface variables'''
    
    pvmaxth = _get_max_PV_on_theta(ds_der['q'])
    ptop = _get_zonal_average_p(ds_der['p'])
    thlev = ds_der.theta
    pvlev = _get_PV_levels(thlev,thref=thref)
    mass_int, circ_int = _compute_mass_and_circulation_integrals(ds_der['r'],ds_der['q'],pvlev,thlev)
    
    ds_bgs_input = xr.Dataset()
    
    return ds_bgs_input

def _compute_mass_and_circulation_integrals(r,q,pvlev,thlev):
    '''Function computes mass and circulation integrals within PV contours on isentropic levels
    Input:
        r - (time,theta,lat,lon) array; 3-d isentropic density
        q - (time,theta,lat,lon) array; 3-d Ertel PV
        pvlev - vector of PV levels
        thlev - vector of theta levels
    Returns:
        mass_int - mass integrals
        circ_int - circulation integrals
    '''
    ntime, _, nlat, nlon = r.shape
    npv = len(pvlev)
    nth = len(thlev)
    
    mass_int = np.zeros([ntime,npv,nth])
    circ_int = np.zeros([ntime,npv,nth])
    
    return mass_int, circ_int

def _get_max_PV_on_theta(q):
    '''Function to get maximum value of PV on each theta level. This is used during post-processing of the background state
    in order to define an upper bound when interpolating the computed background state PV on each isentropic level.
    
    Input:
        q - (time,theta,lat,lon) array; 3-d Ertel PV
    Returns:
        maximum PV on each theta level at each time
    '''
    return np.max(q,axis=(2,3))

def _get_PV_levels(thlev,thref=380.,nlev_st = 56,nlev_tpp = 13,nlev_trop = 5):
    '''Routine for getting PV levels relevant to the specified experiment. Mass and circulation
    integrals will be computed within each corresponding PV contour'''
    # Define evenly spaced PV levels in the stratosphere
    pvlev_st = 5.0 + 1.0 * np.arange(nlev_st, dtype=float)

    # Define logarithmically spaces PV levels in the tropopause
    lnqmin = -1.0
    lnqmax = np.log(pvlev_st[0])
    dlnq = (lnqmax - lnqmin) / float(nlev_tpp)
    lnq = lnqmin + dlnq * np.arange(nlev_tpp, dtype=float)
    pvlev_tpp = np.exp(lnq)

    # Define evely spaced PV levels in the troposphere
    pvlev_trop = (pvlev_tpp[0] / float(nlev_trop)) * np.arange(nlev_trop, dtype=float)

    # Concatenate PV-level arrays
    pvlev = np.concatenate([pvlev_trop, pvlev_tpp, pvlev_st])

    # Define lait-to-pv scaling
    lait2pv = np.ones(len(thlev))
    i_overwrite = np.where(thlev >= thref)
    lait2pv[i_overwrite] = (thlev[i_overwrite] / thref) ** 4.5

    return pvlev, lait2pv

def _get_zonal_average_p(p):
    '''Function to compute zonal average pressure on top theta level at specified latitudes.
    This is used in the background state computation to deal with missing data above the top
    theta level.'''
    ################ Move into interpolation function and do zonal average there to not carry around unnecessary data
    return 

