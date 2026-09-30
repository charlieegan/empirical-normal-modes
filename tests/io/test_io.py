import xarray as xr
import numpy as np
import enm

# Execute the functions
dtime = None #['2024-01-01','12:00']
output_path = None # '/storage/research/s2senm/cq934523/ERA5_model_level_data/'
ds_interior, ds_surface = enm.io.download_ERA5_model_level_data_and_get_pressure_and_theta(dtime=dtime,output_path=output_path,
                                                              kappa=2./7.,p00=101325.,
                                                              coeff_output_filename='era5_l137_coefficients.npz')
