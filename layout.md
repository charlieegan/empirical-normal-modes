empirical-normal-modes/                     <- GitHub repo root
├── pyproject.toml
├── setup.py                                <- custom build_ext, compiles the native extension
├── README.md
├── src/
│   └── enm/
│       ├── __init__.py
│       ├── io/                             <- (1) load data, save to netCDF
│       │   ├── __init__.py
│       │   └── era5_io.py
│       ├── interpolation/                  <- (2) isentropic + finite differencing
│       │   ├── __init__.py
│       │   ├── isentropic_interpolation.py
│       ├── integrals/                      <- (3) isentropic-coordinate integrals
│       │   ├── __init__.py
│       │   └── isentropic_integrals.py
│       ├── background/                     <- (4) thin adapter around the external package
│       │   ├── __init__.py
│       │   └── state.py
│       ├── wave_activity/                  <- (5)
│       │   ├── __init__.py
│       │   └── ...
│       ├── modes/                          <- (6) eigenvalue problem
│       │   ├── __init__.py
│       │   └── eigensolve.py
│       └── _native/                        <- ONE consolidated C++ extension
│           ├── __init__.py
│           └── src/
│               └── isentropic_interp.cpp
└── tests/
    ├── io/
    |   └──test_io.py
    ├── interpolation/
    |   └──test_interpolation.py
    ├── integrals/
    ├── background/
    ├── wave_activity/
    └── modes/