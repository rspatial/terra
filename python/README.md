# terra (Python)

Python version of the **terra** package.

Note that this is an alpha (very early) version! Please report any bugs or suggestions.

---

## Requirements

| Dependency | Minimum version |
|---|---|
| Python | 3.9 |
| CMake | 3.18 |
| C++ compiler | C++17 (MSVC 2019+, GCC 9+, Clang 10+) |
| GDAL | 3.0 |
| GEOS | 3.8 |
| PROJ | 6.0 |
| pybind11 | 2.12 |
| scikit-build-core | 0.9 |
| numpy | 1.20 |
| pandas | 1.3 |
| matplotlib | 3.4 |

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

3. **Install Python build and runtime dependencies**:

   ```powershell
   pip install scikit-build-core numpy pandas matplotlib
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

3. Install Python build and runtime dependencies:

   ```powershell
   pip install scikit-build-core pybind11 numpy pandas matplotlib
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
pip install scikit-build-core pybind11 numpy pandas matplotlib
pip install -e .
```

---

## macOS installation

```bash
brew install gdal geos proj cmake

cd python
python3 -m venv .venv
source .venv/bin/activate
pip install scikit-build-core pybind11 numpy pandas matplotlib
pip install -e .
```

---

## Usage

```python
import terra
import matplotlib.pyplot as plt

# Create / read a raster
f = "https://github.com/rspatial/terra/raw/refs/heads/master/inst/ex/elev.tif"
r = terra.rast(f)
e = terra.ext(6, 6.4, 49.5, 50)
x = r.crop(e)

# Inspect
print(r)

# Plot
terra.plot(r)
plt.show()

# Work with extents
e = terra.ext(-180, 180, -90, 90)
print(e)

# Arithmetic operators
r2 = r * 2
r3 = r + r2
mask = r > 500

# Vector
ff = "https://github.com/rspatial/terra/raw/refs/heads/master/inst/ex/lux.shp"
v = terra.vect(ff)
v
```

---

## License

GPL-3.0-or-later
