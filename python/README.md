# terra (Python)

Python bindings for the **terra** geospatial C++ library, with the same core types and API as the R package "terra".

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

## Usage


```python
import terra as pt

r = pt.rast()                          # default global grid, like R
r = pt.rast("dem.tif")
e = pt.ext(-180, 180, -90, 90)
v = pt.vect(e, crs="EPSG:4326")
pt.crs(r)                              # WKT
pt.crs(r, "+proj=longlat")             # set CRS, returns r
pt.messages(r, "myfn")                 # raise on error / warn
```


The raw C++ interface

```python
import terra as pt

r = pt.SpatRaster(["elevation.tif"], [], [], False, [], [], [], False, False, [])
v = pt.SpatVector()
v.read("countries.gpkg")
ext = pt.SpatExtent(-10, 40, 35, 70)
opt = pt.SpatOptions()
opt.filenames = ["output.tif"]
```

## License

GPL-3.0
