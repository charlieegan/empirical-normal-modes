"""
Reader for legacy TRADV gridded files (``dihead<date>`` / ``diout<date>``) used by the IDL code.

Only used to build a common test input for validating the Python translation
against IDL. The production input is gridded netCDF model output.

File layout (from IDL ``read_trhead`` / ``read_trdata`` with ltracer=0)
-----------------------------------------------------------------------
header (text): date, ntrac, diagnostic type codes, nlon, nlat, latitudes (N->S),
               nlp half levels, A (aeta) and B (beta) coefficients.
data (raw float32 stream, no record markers), tracers in header order:
    type 0 (surface pressure/p00): nlat records of nlon values
    other types: for level in 0..nlev-1 (top first): for lat in 0..nlat-1: nlon values
Type codes: 0=ps/p00, 2=theta, 3=Ertel PV (PVU), 5=u, 6=v, 8=heating.

Unit conventions inferred from the IDL (verify against IDL output):
    U_idl = u_stored * cos^2(lat) / 86400   = (u/a) cos(lat)  [rad/s]
    V_idl = v_stored / 86400                = (v/a) cos(lat)  [rad/s]
    =>  u = u_stored * a cos(lat) / 86400,   v = v_stored * a / (86400 cos(lat))
A and aeta are normalised by p00 = 1e5 Pa; psurf is ps/p00.
"""
#from __future__ import annotations

import os
from dataclasses import dataclass
import re

import numpy as np
import xarray as xr

P00 = 1.0e5            # Pa, normalisation of A and ps in the legacy files
RADEA = 6371299.0      # m
KAPPA = 2.0 / 7.0


@dataclass
class TradvHeader:
    date: int
    itratyp: list
    nlon: int
    nlat: int
    nlev: int
    lat: np.ndarray    # degrees, north -> south
    aeta: np.ndarray   # half levels, A/p00, top -> surface
    beta: np.ndarray   # half levels, B


def read_tradv_header(head_path: str) -> TradvHeader:
    with open(head_path) as f:
        lines = f.read().splitlines()
    date = int(lines[0].split()[-1])
    ntrac = int(lines[2].split()[-1])
    itratyp = [int(x) for x in lines[3].split()[-ntrac:]]
    nlon = int(lines[5].split()[0])
    nlat = int(lines[6].split()[0])
    # lines[4] is the relaxation-timescale line (skipped by the IDL too).
    tok = " ".join(lines[7:]).split()
    lat = np.array(tok[:nlat], dtype=float)
    nlp = int(tok[nlat])
    start = tok.index("B:", nlat) + 1
    aeta = np.array(tok[start:start + nlp], dtype=float)
    beta = np.array(tok[start + nlp:start + 2 * nlp], dtype=float)
    if len(lat) != nlat or len(aeta) != nlp or len(beta) != nlp:
        raise ValueError("Header parsing failed: unexpected number of values")
    return TradvHeader(date, itratyp, nlon, nlat, nlp - 1, lat, aeta, beta)


def expected_nbytes(h: TradvHeader) -> int:
    n = sum(h.nlat * h.nlon if t == 0 else h.nlev * h.nlat * h.nlon for t in h.itratyp)
    return 4 * n


def _plausible_theta(x) -> bool:
    return bool(np.all(np.isfinite(x)) and 150.0 < np.median(x) < 5000.0)


def _split_fields(raw, h):
    out, pos = {}, 0
    for t in h.itratyp:
        if t == 0:
            n, shape = h.nlat * h.nlon, (h.nlat, h.nlon)
        else:
            n, shape = h.nlev * h.nlat * h.nlon, (h.nlev, h.nlat, h.nlon)
        out[t] = raw[pos:pos + n].reshape(shape)
        pos += n
    return out


def read_tradv(head_path: str, data_path: str, date: int = None, endian: str = "auto",
               p00: float = P00, a: float = RADEA, kappa: float = KAPPA) -> xr.Dataset:
    """Read one TRADV timestep into a single SI-unit xarray Dataset.

    Variables: theta, t, u, v, p (time, model_level, latitude, longitude);
    pv (PVU, as stored); ps (Pa) and lnsp (time, latitude, longitude);
    a_half, b_half (Pa and dimensionless, on half_level).
    """
    h = read_tradv_header(head_path)
    if date is None:
        m = re.search(r"(\d{10})$", os.path.basename(data_path))
        if m is None:
            raise ValueError("Cannot infer date from file name; pass date=YYYYMMDDHH. "
                             "(The date inside the header is unreliable.)")
        date = int(m.group(1))
    nbytes = os.path.getsize(data_path)
    if nbytes != expected_nbytes(h):
        raise ValueError(f"{data_path}: {nbytes} bytes, expected {expected_nbytes(h)} "
                         f"for types {h.itratyp} on {h.nlev}x{h.nlat}x{h.nlon}")

    dtypes = {"little": "<f4", "big": ">f4"}
    order = ["little", "big"] if endian == "auto" else [endian]
    for e in order:
        f = _split_fields(np.fromfile(data_path, dtype=dtypes[e]), h)
        if 2 in f and _plausible_theta(f[2]):
            break
    else:
        raise ValueError("Could not find a byte order giving plausible theta values")
    f = {k: v.astype(np.float32) for k, v in f.items()}

    lat = h.lat
    lon = np.arange(h.nlon) * 360.0 / h.nlon
    coslat = np.cos(np.deg2rad(lat))[None, :, None]

    ps = f[0].astype(np.float64) * p00                       # Pa
    a_half = h.aeta * p00                                     # Pa
    p_half = a_half[:, None, None] + h.beta[:, None, None] * ps[None]
    p_full = 0.5 * (p_half[:-1] + p_half[1:])                 # Pa, top -> surface

    theta = f[2]
    u = f[5] * a * coslat / 86400.0
    v = f[6] * a / (86400.0 * coslat)
    t = theta * (p_full / p00) ** kappa

    time = np.array([_date_to_dt64(int(date))])
    d3 = ("time", "model_level", "latitude", "longitude")
    d2 = ("time", "latitude", "longitude")
    ds = xr.Dataset(
        {
            "theta": (d3, theta[None], {"units": "K"}),
            "t": (d3, t[None].astype(np.float32), {"units": "K"}),
            "u": (d3, u[None].astype(np.float32), {"units": "m s-1"}),
            "v": (d3, v[None].astype(np.float32), {"units": "m s-1"}),
            "p": (d3, p_full[None], {"units": "Pa"}),
            "ps": (d2, ps[None], {"units": "Pa"}),
            "lnsp": (d2, np.log(ps)[None], {"units": "1"}),
            "a_half": (("half_level",), a_half, {"units": "Pa"}),
            "b_half": (("half_level",), h.beta),
        },
        coords={"time": time, "model_level": np.arange(1, h.nlev + 1),
                "latitude": lat, "longitude": lon,
                "half_level": np.arange(h.nlev + 1)},
        attrs={"source": "legacy TRADV", "p00_Pa": p00,
               "wind_conversion": "inferred from IDL; verify"},
    )
    if 3 in f:
        ds["pv"] = (d3, f[3][None], {"units": "PVU (as stored)"})
    return ds


def _date_to_dt64(date: int) -> np.datetime64:
    s = f"{date:010d}"
    return np.datetime64(f"{s[:4]}-{s[4:6]}-{s[6:8]}T{s[8:10]}:00")


def to_current_repo_format(ds: xr.Dataset):
    """Split into (ds_int, ds_surf) as used by apply_isentropic_interpolation
    (surface lnsp carries a length-1 model_level dim, as in the MARS files)."""
    ds_int = ds[["theta", "t", "u", "v", "p"]].copy()
    ds_surf = ds[["lnsp"]].expand_dims(model_level=[1], axis=1).copy()
    return ds_int, ds_surf


def legacy_to_netcdf(head_path: str, data_path: str, out_path: str) -> xr.Dataset:
    ds = read_tradv(head_path, data_path)
    ds.to_netcdf(out_path)
    return ds