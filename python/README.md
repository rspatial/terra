# pyterra

Python bindings for the **terra** spatial C++ library, exposed with the same core types as the R package (`SpatRaster`, `SpatVector`, etc.).

## Build (Linux / WSL)

Use a virtual environment (recommended on Ubuntu, which uses PEP 668 for the system Python):

```bash
cd python
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
```

You need CMake, a C++17 compiler, and development packages for GDAL, GEOS, and PROJ (e.g. `libgdal-dev`, `libgeos-dev`, `libproj-dev` on Debian/Ubuntu).

## License

GPL-3.0-or-later (same as terra).
