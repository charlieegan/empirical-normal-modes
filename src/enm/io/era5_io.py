import cdsapi
import eccodes
import xarray as xr
import numpy as np
import os

def download_ERA5_model_level_data_and_get_pressure_and_theta(dtime=None,output_path=None,
                                                              kappa=2./7.,p00=101325.,
                                                              coeff_output_filename='era5_l137_coefficients.npz'):
    
    _download_ERA5_model_level_data(dtime=dtime,output_path=output_path)
    _save_ERA5_model_level_coefficients(output_path=output_path,output_filename=coeff_output_filename)
    ds_interior, ds_surface = _get_model_level_pressure_and_potential_temp(dtime=dtime,output_path=output_path,kappa=kappa,p00=p00)
    
    return ds_interior, ds_surface

def _download_ERA5_model_level_data(dtime=None,output_path=None):
    '''Function loads ERA5 model level data for u, v, T on {date} at {time} and saves as netCDF.
    Separately, surface variables p and z are also loaded and saved as netCDF.
    dtime - list specifying date and time, e.g. ['2024-01-01','12:00']
    output_path - string; where to save output netCDF files
    '''
    
    # define where to save the data (may want to change this to use a temporary file that is immediately deleted)
    if output_path is None:
        output_path = '/storage/research/s2senm/cq934523/ERA5_model_level_data/'
    if dtime is None:
        date = '2024-01-01'
        time = '12:00'
        output_filename_interior = 'ERA5_single_step_test_interior.nc'
        output_filename_surface = 'ERA5_single_step_test_surface.nc'
    else:
        date = dtime[0]
        time = dtime[1]
        dtime_str = f'{date.replace('-','')}{time.replace(':','')}'
        output_filename_interior = f'ERA5_interior_{dtime_str}.nc'
        output_filename_surface = f'ERA5_surface_{dtime_str}.nc'
    
    int_ex = os.path.exists(output_path+output_filename_interior)
    surf_ex = os.path.exists(output_path+output_filename_surface)
    if int_ex and surf_ex:
        print(f"Files '{output_filename_interior}' and {output_filename_surface} already exist. Skipping download.")
        return
    
    c = cdsapi.Client()

    # Generate the slash-separated string for all 137 levels
    LEVEL_LIST = '/'.join([str(i) for i in range(1, 138)])
    
    # Base parameters shared by requests for both interior and surface variables
    base_params = {
        'class': 'ea',
        'expver': '1',
        'stream': 'oper',
        'type': 'an',
        'levtype': 'ml',
        'date': date,
        'time': time,
        'grid': '0.25/0.25',
        'format': 'netcdf',
        'area': '90/-180/-90/180'
    }

    # Download the interior variables t, u and v on all model levels
    print(f"Submitting MARS request for t, u and v on all model levels on {date} at {time} UTC")
    params_3d = {**base_params, 'levelist': LEVEL_LIST, 'param': '130/131/132'}
    c.retrieve('reanalysis-era5-complete', params_3d, output_path+output_filename_interior)

    # Download surface variables z, lnsp on model level 1
    print("Downloading 2D surface variables (z, lnsp)...")
    params_2d = {**base_params, 'levelist': '1', 'param': '129/152'}
    c.retrieve('reanalysis-era5-complete', params_2d, output_path+output_filename_surface)
        
def _save_ERA5_model_level_coefficients(output_path=None,output_filename='era5_l137_coefficients.npz'):
    """
    Checks if the coefficients file exists. If missing, it downloads a minimal
    GRIB file, extracts the L137 A and B coefficients via ecCodes, and saves them.
    """
    if output_path is None:
        output_path = '/storage/research/s2senm/cq934523/ERA5_model_level_data/'
    output_filename = output_path+output_filename
    
    # Check if the file already exists
    if os.path.exists(output_filename):
        print(f"File '{output_filename}' already exists. Skipping download.")
        return

    print(f"'{output_filename}' not found. Initialising download script...")
    c = cdsapi.Client()
    temp_grib = 'temp_metadata_anchor.grib'
    
    try:
        # Download minimal data
        c.retrieve('reanalysis-era5-complete', {
            'class': 'ea',
            'expver': '1',
            'stream': 'oper',
            'type': 'an',
            'levtype': 'ml',
            'levelist': '1',
            'param': '130',
            'date': '2024-01-01',
            'time': '12:00',
            'grid': '0.25/0.25',
            'data_format': 'grib',
            'area': [90, -180, 89.75, -179.75] 
        }, temp_grib)
        
        # Extract the underlying pv array directly using low-level ecCodes
        print("Parsing binary grid coordinate structures with ecCodes...")
        with open(temp_grib, 'rb') as f:
            gid = eccodes.codes_grib_new_from_file(f)
            if gid is None:
                raise IOError("Failed to parse the downloaded GRIB file.")
            
            # Extract the raw structural vertical coordinate array
            pv_array = eccodes.codes_get_array(gid, 'pv')
            eccodes.codes_release(gid)

        # Split the pv array into A and B segments (L137 structure = 138 entries each)
        num_boundaries = 138
        a_coeff = np.array(pv_array[:num_boundaries])
        b_coeff = np.array(pv_array[num_boundaries:])
        
        # 5. Save output
        np.savez(output_filename, a=a_coeff, b=b_coeff)
        print(f"Coefficients saved to '{output_filename}'")
        
    except Exception as e:
        print(f"An error occurred during retrieval/processing: {e}")
        raise e
        
    finally:
        # Remove temporary files
        if os.path.exists(temp_grib):
            os.remove(temp_grib)
        
def _get_model_level_pressure_and_potential_temp(dtime=None,output_path=None,kappa=2./7.,p00=101325.):
    '''Function recovers pressure and potential temperature on model levels using surface pressure,
    model level coefficients A and B, and model level temperature. Potential temperature is defined as
    
    Theta = T(p_r/p)^kappa
    
    where p_r is reference pressure, and kappa=R/c_p with R the gas constant and c_p specific heat capacity at constant pressure.
    '''
    if output_path is None:
        output_path = '/storage/research/s2senm/cq934523/ERA5_model_level_data/'
    if dtime is None:
        date = '2024-01-01'
        time = '12:00'
        filename_interior = 'ERA5_single_step_test_interior.nc'
        filename_surface = 'ERA5_single_step_test_surface.nc'
    else:
        date = dtime[0]
        time = dtime[1]
        dtime_str = f'{date.replace('-','')}{time.replace(':','')}'
        filename_interior = f'ERA5_interior_{dtime_str}.nc'
        filename_surface = f'ERA5_surface_{dtime_str}.nc'
        
    ds_interior = xr.open_dataset(output_path+filename_interior)
    
    if 'p' in ds_interior.data_vars and 'theta' in ds_interior.data_vars:
        print('Pressure and potential temperature already computed. Skipping computation.')
        return
        
    ################# Get pressure on model levels #######################
    # Load coefficients that define model levels from surface pressure
    coeffs = np.load(output_path+'/era5_l137_coefficients.npz')
    a_half_vals = coeffs['a']  # 138 entries for half-level boundaries
    b_half_vals = coeffs['b']  # 138 entries for half-level boundaries
    
    # Load surface pressure
    ds_surface = xr.open_dataset(output_path+filename_surface)
    p_surf = np.exp(ds_surface['lnsp'])

    # Drop any dummy vertical dimensions from the 2D surface field for broadcasting
    if 'model_level' in p_surf.dims:
        p_surf = p_surf.squeeze('model_level', drop=True)
    elif 'hybrid' in p_surf.dims:
        p_surf = p_surf.squeeze('hybrid', drop=True)

    # Convert the 1D arrays into DataArrays with a explicit 'half_level' dimension
    a_half = xr.DataArray(a_half_vals, dims=['half_level'])
    b_half = xr.DataArray(b_half_vals, dims=['half_level'])

    # Compute 3-d pressure
    ## Compute the pressure at all 138 half-level interfaces
    p_half = a_half + b_half * p_surf

    # Calculate full-level pressure as mean of bounding half-levels
    p_full = np.swapaxes(
                0.5 * (p_half.isel(half_level=slice(0, 137)).values + 
                    p_half.isel(half_level=slice(1, 138)).values)
                         , 0, 1)
    
    ############ Get potential temperature on model levels ##################
    theta_full = ds_interior['t']*(p00/p_full)**kappa

    ########### Add model-level pressure and potential temperature to the dataset and save to netCDF #############
    ds_interior['p'] = xr.DataArray(
        p_full, 
        coords=ds_interior['t'].coords,
        dims=ds_interior['t'].dims
    )
    ds_interior['p'].attrs = {'units': 'Pa', 'long_name': 'Pressure'}
    
    ds_interior['theta'] = xr.DataArray(
            theta_full, 
            coords=ds_interior['t'].coords,
            dims=ds_interior['t'].dims
        )
    ds_interior['theta'].attrs = {'units': 'K', 'long_name': 'Potential temperature'}
    os.remove(output_path+filename_interior)
    ds_interior.to_netcdf(output_path+filename_interior)

    print("\n3-d Model-level pressure ('p') and potential temperature ('theta') calculated and merged")
    
    return ds_interior, ds_surface

