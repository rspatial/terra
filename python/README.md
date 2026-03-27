# terra (Python)

Python bindings for the **terra** geospatial C++ library, with the same core
types and API as the R package "terra".

---

## Requirements

| Tool | Minimum version |
|------|----------------|
| Python | 3.9 |
| CMake | 3.18 |
| C++ compiler | C++17 (MSVC 2019+, GCC 9+, Clang 10+) |
| GDAL | 3.0 |
| GEOS | 3.8 |
| PROJ | 6.0 |
| pybind11 | 2.12 |
| scikit-build-core | 0.9 |

Optional: **numpy**, **pandas** (for array / data-frame interop),
**matplotlib** (for `terra.plot()`).

---

## Windows installation

### Option A — conda / mamba (recommended)

conda from [miniforge](https://github.com/conda-forge/miniforge/releases) or
[Anaconda](https://www.anaconda.com/download) provides pre-built GDAL, GEOS,
and PROJ binaries with matching CMake packages, which is by far the easiest
route on Windows.

1. **Install Visual Studio Build Tools** (free) with the
   *Desktop development with C++* workload from
   <https://visualstudio.microsoft.com/downloads/>.
   You need MSVC and the Windows SDK.

2. **Create a conda environment** and install the geospatial dependencies:

   ```powershell
   conda create -n terra python=3.12
   conda activate terra
   conda install -c conda-forge gdal geos proj pybind11 cmake ninja
   ```

3. **Install scikit-build-core** (the Python build backend):

   ```powershell
   pip install scikit-build-core
   ```

4. **Build and install** the package from the `python/` directory of the
   terra source tree:

   ```powershell
   cd C:\path\to\terra\python
   pip install -e .
   ```

   scikit-build-core runs CMake automatically and places the compiled
   `_terra.pyd` extension inside `src/terra/`.

5. **Verify**:

   ```python
   import terra as pt
   r = pt.rast()
   print(r)
   ```

> **Tip:** if CMake cannot find the conda libraries, pass the prefix
> explicitly:
> ```powershell
> pip install -e . -C cmake.define.CMAKE_PREFIX_PATH="%CONDA_PREFIX%\Library"
> ```

---

### Option B — vcpkg

[vcpkg](https://github.com/microsoft/vcpkg) is Microsoft's C++ package
manager and is pre-installed on GitHub Actions Windows runners.

1. Install Visual Studio Build Tools (see step 1 above).

2. Install vcpkg and the geospatial libraries:

   ```powershell
   git clone https://github.com/microsoft/vcpkg C:\vcpkg
   C:\vcpkg\bootstrap-vcpkg.bat
   C:\vcpkg\vcpkg install gdal geos proj --triplet x64-windows
   ```

3. Install Python build dependencies:

   ```powershell
   pip install scikit-build-core pybind11
   ```

4. Build the package, pointing CMake to the vcpkg toolchain:

   ```powershell
   cd C:\path\to\terra\python
   pip install -e . `
     -C cmake.define.CMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake `
     -C cmake.define.VCPKG_TARGET_TRIPLET=x64-windows
   ```

> **DLL search path:** when using vcpkg you must add the vcpkg DLL directory
> to `PATH` before importing terra so that `gdal.dll`, `geos_c.dll`, and
> `proj.dll` are found:
> ```powershell
> $env:PATH = "C:\vcpkg\installed\x64-windows\bin;" + $env:PATH
> ```
> This is not needed with conda, where the environment's `Library\bin` is
> already on `PATH`.

---

## Linux installation

```bash
# Debian / Ubuntu
sudo apt-get install -y libgdal-dev libgeos-dev libproj-dev cmake

# Fedora / RHEL
sudo dnf install -y gdal-devel geos-devel proj-devel cmake

cd python
python3 -m venv .venv
source .venv/bin/activate
pip install scikit-build-core pybind11
pip install -e .
```

---

## macOS installation

```bash
brew install gdal geos proj cmake

cd python
python3 -m venv .venv
source .venv/bin/activate
pip install scikit-build-core pybind11
pip install -e .
```

---

## Usage

```python
import terra as pt
import matplotlib.pyplot as plt

# Create / read a raster
r = pt.rast("elevation.tif")

# Inspect
print(r)                         # shows class, size, resolution, CRS, …

# Plot (requires matplotlib)
pt.plot(r)
plt.show()

# Work with extents
e = pt.ext(-180, 180, -90, 90)
print(e)

# Arithmetic operators (after import operators are registered automatically)
r2 = r * 2
r3 = r + r2
mask = r > 500

# Vector
v = pt.vect("countries.gpkg")
v_crop = pt.crop_vect(v, e)
```

---

## Optional dependencies

Install extras for array/data-frame interoperability and plotting:

```bash
pip install "terra[geo]"          # numpy + pandas
pip install matplotlib            # for pt.plot() and pt.plot_rgb()
```

---

## License

GPL-3.0-or-later
