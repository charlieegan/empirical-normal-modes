"""
Thin ctypes loader for the consolidated native extension
(``libenm_native.so``), built by ``setup.py`` from every ``.cpp`` file under
``_native/src/``.

Other modules should import from here rather than touching ctypes directly.
If the shared library hasn't been built (e.g. no g++/OpenMP at install
time), ``load()`` returns ``None`` and callers are expected to fall back to
their pure-NumPy implementation -- see enm.interpolation.isentropic for the
pattern.
"""
import ctypes
from pathlib import Path

_SO_PATH = Path(__file__).parent / "src" / "libenm_native.so"
_lib = None
_load_attempted = False


def load():
    """Return the loaded native library (a ctypes.CDLL), or None if it
    hasn't been built."""
    global _lib, _load_attempted
    if not _load_attempted:
        _load_attempted = True
        if _SO_PATH.exists():
            _lib = ctypes.CDLL(str(_SO_PATH))
            _configure_signatures(_lib)
    return _lib


def _configure_signatures(lib):
    """Set argtypes/restype for each native kernel we expose, so ctypes
    validates arguments instead of silently misinterpreting them."""
    f32p = ctypes.POINTER(ctypes.c_float)
    lib.enm_interp_theta_levels.argtypes = [
        f32p,          # th
        f32p,          # vars
        ctypes.c_long,  # nlev
        ctypes.c_long,  # ncols
        ctypes.c_long,  # nvars
        f32p,          # thlev
        ctypes.c_long,  # ntarget
        f32p,          # out
        ctypes.c_float,  # fill_value
    ]
    lib.enm_interp_theta_levels.restype = None
