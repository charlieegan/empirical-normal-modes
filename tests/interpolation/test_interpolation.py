import numpy as np
import xarray as xr
import enm
import matplotlib.pyplot as plt
import time

# Load ERA5 data
input_path = '/storage/research/s2senm/cq934523/ERA5_model_level_data/'
filename_int = 'ERA5_single_step_test_interior.nc'
filename_surf = 'ERA5_single_step_test_surface.nc'

ds_int = xr.open_dataset(input_path+filename_int)
ds_surf = xr.open_dataset(input_path+filename_surf)

vars_dict = enm.interpolation.apply_isentropic_interpolation(ds_int,ds_surf)

q = vars_dict['q']
plt.imshow(q[0,:75,:180,0])




