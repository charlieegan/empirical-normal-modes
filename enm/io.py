import xarray as xr
from ecmwf.opendata import Client

# Initialize the ECMWF Open Data Client
client = Client()

print("Downloading model level data from ECMWF...")
# 1. Download atmospheric variables on model levels (ml)
# We request 't' (temperature), 'u' (zonal wind), 'v' (meridional wind)
client.retrieve(
    type="fc",              # Forecast data
    levtype="ml",           # Model levels (Hybrid sigma-pressure coordinates)
    levelist=[130, 137],    # Example: Lower model levels (IFS has 137 levels total)
    param=["t", "u", "v"],  # Base archived parameters
    step=12,                # Forecast step hour
    target="atmosphere_ml.grib"
)

# 2. Download mandatory surface reference parameters needed for pressure/geopotential calculations
client.retrieve(
    type="fc",
    levtype="ml",
    levelist=1,             # Surface reference fields are historically assigned to level 1
    param=["lnsp", "z"],    # Logarithm of surface pressure & Surface Geopotential
    step=12,
    target="surface_ml.grib"
)
print("Downloads complete!")

# ----------------------------------------------------
# Reading the GRIB data using xarray + cfgrib engine
# ----------------------------------------------------
print("\nOpening datasets with xarray...")

# Load atmospheric data
# Using filter_by_keys protects xarray from dimension conflicts
ds_atmos = xr.open_dataset('atmosphere_ml.grib', engine='cfgrib', 
                            backend_kwargs={'filter_by_keys': {'typeOfLevel': 'hybrid'}})

# Load surface reference data
ds_surf = xr.open_dataset('surface_ml.grib', engine='cfgrib',
                           backend_kwargs={'filter_by_keys': {'typeOfLevel': 'hybrid'}})

print("\n--- Atmospheric Data Structure ---")
print(ds_atmos)

print("\n--- Surface Base Structure ---")
print(ds_surf)

# Access your downloaded U and V fields directly
u_wind = ds_atmos['u']
v_wind = ds_atmos['v']
temperature = ds_atmos['t']

print(f"\nSuccessfully read U-wind shape: {u_wind.shape}")