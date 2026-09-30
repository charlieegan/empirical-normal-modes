"""
Builds src/enm/_native/src/libenm_native.so as part of `pip install`.

The native kernels are plain shared-library functions (extern "C", raw
pointers) called via ctypes -- not conventional Python C-extensions -- so
instead of a normal Extension() build we register a dummy Extension purely
to make setuptools schedule a "build_ext" step, and override that step to
compile every .cpp file under _native/src/ into ONE shared library with
g++. Add a new kernel by dropping a .cpp file in that directory; it's
picked up automatically (no changes needed here).
"""
import subprocess
import sys
from pathlib import Path

from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext as _build_ext

ROOT = Path(__file__).parent
NATIVE_DIR = ROOT / "src" / "enm" / "_native"
SRC_DIR = NATIVE_DIR / "src"
OUT = SRC_DIR / "libenm_native.so"


def compile_cmd(sources):
    return [
        "g++", "-O3", "-march=native", "-fopenmp", "-fPIC", "-shared", "-std=c++17",
        "-o", str(OUT), *sources,
    ]


class build_ext(_build_ext):
    def run(self):
        sources = sorted(str(p) for p in SRC_DIR.glob("*.cpp"))
        if not sources:
            print(
                f"[enm] No .cpp sources found under {SRC_DIR}; skipping native build.",
                file=sys.stderr,
            )
            return
        try:
            subprocess.run(compile_cmd(sources), check=True)
        except (subprocess.CalledProcessError, FileNotFoundError) as exc:
            print(
                f"\n[enm] WARNING: native build step failed ({exc}).\n"
                "[enm] Falling back to pure-NumPy implementations at runtime.\n"
                "[enm] To get the fast backend, install g++ with OpenMP support "
                "and re-run:\n"
                "[enm]   pip install -e . --no-build-isolation --force-reinstall --no-deps\n",
                file=sys.stderr,
            )

    def get_outputs(self):
        return [str(OUT)] if OUT.exists() else []


setup(
    name="enm",
    version="0.1.0",
    description="Empirical normal mode analysis of the atmosphere",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    package_data={"enm._native": ["src/*.cpp", "src/*.so"]},
    include_package_data=True,
    install_requires=["numpy", "scipy"],
    extras_require={"xarray": ["xarray"]},
    python_requires=">=3.9",
    ext_modules=[Extension("enm._native._trigger", sources=[])],
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
)
