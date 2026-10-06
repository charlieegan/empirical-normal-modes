"""
Compare raw/input fields between the IDL legacy code and the Python translation,
for the same TRADV timestep.

Loads:
  - idl_output_2010022218.sav : IDL `save`d variables (scipy.io.readsav; variable
    names come back lowercased).
  - python_output_2010012218.npz : arrays saved by run_python_interp_2010012218.py.

NOTE: as of writing, neither saved file actually contains u, v, t, theta, p or ps
-- the IDL `save` statement only lists the derived isentropic fields (qth, zetath,
denth, ...) and the .npz only lists apply_isentropic_interpolation's output dict.
This script will report "not found" for every variable until both sides are
updated to also save the raw fields compared here. See the bottom of this file
for exactly what each side needs to add.

Run with: python tests/interpolation/test_interpolation_against_legacy.py
"""
import os

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable

import numpy as np
from scipy.io import readsav

# Diverging pair (blue <-> red) with a neutral gray midpoint, equal weight per arm.
DIVERGING_CMAP = LinearSegmentedColormap.from_list(
    "blue_red_diverging", ["#2a78d6", "#f0efec", "#e34948"]
)


def colorbar(im, ax):
    """Colorbar matched to the height of im's axes, regardless of the image's aspect ratio
    (plt.colorbar(im, ax=ax) sizes off the full axes bounding box, not the image extent, so it
    comes out oversized whenever imshow's own aspect makes the image shorter than its axes)."""
    cax = make_axes_locatable(ax).append_axes("right", size="5%", pad=0.05)
    return plt.colorbar(im, cax=cax)


def imshow_diverging(ax, data):
    """imshow with the blue/red diverging colormap, symmetric about zero so the neutral midpoint
    always lands exactly on 0 regardless of how skewed the data's own min/max are."""
    vmax = np.nanmax(np.abs(data))
    return ax.imshow(data, cmap=DIVERGING_CMAP, vmin=-vmax, vmax=vmax)

HERE = os.path.dirname(os.path.abspath(__file__))
IDL_SAV_PATH = os.path.join(HERE, "idl_output_2010022218.sav")
PYTHON_NPZ_PATH = os.path.join(HERE, "python_output_2010012218.npz")

P00 = 1.0e5   # Pa, IDL normalises ps (and p) by p00
RADEA = 6371299.0

plot_figs = False # control whether or not to plot figures comparing outputs


def compare(name, idl_arr, py_arr, rtol=1e-3, atol=0):
    if idl_arr is None:
        print(f"[{name}] SKIPPED -- not found in IDL output")
        return
    if py_arr is None:
        print(f"[{name}] SKIPPED -- not found in Python output")
        return
    idl_arr = np.asarray(idl_arr, dtype=np.float64)
    py_arr = np.asarray(py_arr, dtype=np.float64)
    if idl_arr.shape != py_arr.shape:
        print(f"[{name}] SHAPE MISMATCH: idl={idl_arr.shape} python={py_arr.shape}")
        return
    diff = py_arr - idl_arr
    finite = np.isfinite(diff)
    if not finite.any():
        print(f"[{name}] no finite values to compare")
        return
    diff = diff[finite]
    close = np.isclose(py_arr, idl_arr, rtol=rtol, atol=atol, equal_nan=True)
    frac_close = close[finite].mean()
    print(f"[{name}] shape={idl_arr.shape}  max|diff|={np.abs(diff).max():.6g}  "
          f"mean|diff|={np.abs(diff).mean():.6g}  "
          f"fraction within rtol={rtol:g}: {frac_close:.4%}")

if not os.path.exists(IDL_SAV_PATH):
    raise FileNotFoundError(f"IDL output not found: {IDL_SAV_PATH}")
if not os.path.exists(PYTHON_NPZ_PATH):
    raise FileNotFoundError(f"Python output not found: {PYTHON_NPZ_PATH}")

idl = readsav(IDL_SAV_PATH)
py = np.load(PYTHON_NPZ_PATH)

idl['r'] = idl.pop('denth')
idl['zeta_unmod'] = idl.pop('zetath')
idl['zeta'] = idl.pop('zetamod')
idl['q'] = idl.pop('qth')

print(f"IDL variables:    {sorted(idl.keys())}")
print(f"Python variables: {sorted(py.files)}")
print()

# Compare the derived variables r, zeta and q
derived_var_names = ['r','zeta','q']

lats = py['latitude']
thlev = py['thlev']

for var_name in derived_var_names:
    
    idl_var = idl[str.lower(var_name)].copy() if str.lower(var_name) in idl else None
    py_var = py[var_name][0] if var_name in py.files else None
    compare(var_name, idl_var, py_var, atol=0)

    if plot_figs:

        idl_var[idl_var==0] = np.nan # mask out values outside of the domain

        fig, ax = plt.subplots(1,3,dpi=400,figsize=(12,4),layout='tight')
        
        vmax = np.nanmax(np.abs(idl_var[:,:,0]))
        im0 = ax[0].contourf(lats,thlev,idl_var[:,:,0],cmap=DIVERGING_CMAP, vmin=-vmax, vmax=vmax)
        ax[0].set_title('IDL '+var_name)
        colorbar(im0,ax[0])
        
        im1 = ax[1].contourf(lats,thlev,py_var[:,:,0],cmap=DIVERGING_CMAP, vmin=-vmax, vmax=vmax)
        ax[1].set_title('Python '+var_name)
        colorbar(im1,ax[1])
        
        log_rel_err = np.log(np.abs((idl_var[:,:,0]-py_var[:,:,0]))/np.abs(idl_var[:,:,0]))
        im2 = ax[2].contourf(lats,thlev,log_rel_err,cmap='coolwarm')
        ax[2].set_title(var_name+' log(relative error)')
        colorbar(im2,ax[2])
        
        if var_name == 'r':
            for axi in ax:
                axi.set_yscale('log')
    