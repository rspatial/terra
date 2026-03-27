# terra (Python)

Python version of the **terra** pacakage for spatial data handling.

## Build (Linux / WSL)

Use a virtual environment (recommended on Ubuntu):

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
import terra as geo

r = geo.rast()              
e = geo.ext(-100, 100, -60, 60)
geo.crop(r, e)

v = geo.vect(e, crs="EPSG:4326")
geo.crs(v, proj4=True)
```


## License

GPL-3.0
