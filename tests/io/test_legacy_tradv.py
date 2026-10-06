import numpy as np
import os
from enm.io.legacy_tradv import legacy_to_netcdf

# run test
dtime = '2010012218'

HEAD_PATH = f"/storage/research/diamet/swrmethn/inv3/tradv/diag_2009/dihead{dtime}"
DATA_PATH = f"/storage/research/diamet/swrmethn/inv3/tradv/diag_2009/diout{dtime}"
OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), f"legacy_input.nc")

ds = legacy_to_netcdf(HEAD_PATH, DATA_PATH,OUT_PATH)
print('Legacy output saved in netCDF format at ',OUT_PATH)
