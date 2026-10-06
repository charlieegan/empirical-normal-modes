"""
Generate the Python-side output for comparison against the IDL run on
2010-02-22 18Z, for the interior fields in /storage/research/diamet/
swrmethn/inv3/tradv/diag_2009/{dihead,diout}2010022218.

Run with: python tests/interpolation/run_2010022218.py
"""
import os

import numpy as np

from enm.io.legacy_tradv import read_tradv, to_current_repo_format
from enm.interpolation import apply_isentropic_interpolation

date=2010012218

HEAD_PATH = f"/storage/research/diamet/swrmethn/inv3/tradv/diag_2009/dihead{date}"
DATA_PATH = f"/storage/research/diamet/swrmethn/inv3/tradv/diag_2009/diout{date}"
OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), f"python_output_{date}.npz")

if __name__ == "__main__":
    ds = read_tradv(HEAD_PATH, DATA_PATH, date=date)
    ds_int, ds_surf = to_current_repo_format(ds)

    vars_dict = apply_isentropic_interpolation(ds_int, ds_surf, exp_type="ERA5")

    for name, arr in vars_dict.items():
        print(f"{name}: shape={np.shape(arr)} dtype={np.asarray(arr).dtype}")

    np.savez(OUT_PATH, **vars_dict,
             latitude=ds_int.latitude.values, longitude=ds_int.longitude.values)
    print(f"Saved to {OUT_PATH}")
