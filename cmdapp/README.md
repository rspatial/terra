# cmdapp — standalone terra C++ library and demo

This directory builds the [terra](https://github.com/rspatial/terra) C++ spatial
library as a **standalone static library** (`libterra.a`) and links a small
command-line demo against it — no R or Rcpp required.

## Prerequisites

* A C++17 compiler (`g++`, `clang++`, or MSVC)
* [GDAL](https://gdal.org/) development files
* [PROJ](https://proj.org/) development files
* [GEOS](https://libgeos.org/) C API (`geos_c`)
* GNU Make (Linux / macOS / MSYS2) **or** CMake-compatible toolchain on Windows

### Linux (Debian / Ubuntu)

```bash
sudo apt install build-essential gdal-bin libgdal-dev libproj-dev libgeos-dev
```

### macOS (Homebrew)

```bash
brew install gdal geos proj pkg-config
```

### Windows

There are several options for getting the required toolchain and libraries on
Windows.

**Conda / Miniforge (recommended)**

[Miniforge](https://github.com/conda-forge/miniforge) provides a lightweight
conda distribution with the conda-forge channel pre-configured. From a
Miniforge or Anaconda prompt:

```bash
conda create -n terra-dev -c conda-forge compilers make gdal libgdal geos proj pkg-config
conda activate terra-dev
```

This gives you `gcc`/`g++`, `make`, `gdal-config`, and `pkg-config` — the
Makefile works as-is inside the activated environment. On Windows this uses the
MSVC-compatible `m2w64` or UCRT toolchain that conda-forge ships.

**MSYS2**

Install [MSYS2](https://www.msys2.org/) and open an **MSYS2 UCRT64** shell:

```bash
pacman -S mingw-w64-ucrt-x86_64-gcc \
          mingw-w64-ucrt-x86_64-gdal \
          mingw-w64-ucrt-x86_64-geos \
          mingw-w64-ucrt-x86_64-proj \
          mingw-w64-ucrt-x86_64-pkg-config \
          make
```

**Other options**

* [vcpkg](https://vcpkg.io/) or [OSGeo4W](https://trac.osgeo.org/osgeo4w/) —
  install GDAL/PROJ/GEOS, then compile with the Visual Studio developer command
  prompt (you will need to adjust include/library paths manually — see
  [Linking your own program](#linking-your-own-program)).
* WSL (Windows Subsystem for Linux) — follow the Linux instructions above
  inside the WSL distribution.

## Building

```bash
cd cmdapp
make            # build library + ./terra executable
make -j4        # parallel build
```

On Linux you can use `make -j$(nproc)`, on macOS `make -j$(sysctl -n hw.ncpu)`.

Build artefacts go into `build/`:

| Target | Output | Description |
|--------|--------|-------------|
| `make lib` | `build/libterra.a` | Static library (all `src/*.cpp` except Rcpp glue) |
| `make` or `make all` | `./terra` (or `terra.exe`) | Demo executable linked against the static lib |
| `make clean` | | Remove all build artefacts |

There is also a thin shell wrapper `compile.sh`:

```bash
./compile.sh -j4        # same as: cd cmdapp && make -j4
./compile.sh clean      # same as: cd cmdapp && make clean
```

### Notes for Windows

The same Makefile works from an **MSYS2 UCRT64** terminal or an activated
**conda** environment (see [Prerequisites](#prerequisites)). GDAL/PROJ data and
plugin paths are baked in at compile time via `gdal-config`, so the executable
is self-contained as long as the toolchain's DLLs are on `PATH` (they are by
default inside MSYS2 or an activated conda env).

To run `terra.exe` outside these shells, make sure the relevant `bin/` directory
(MSYS2 UCRT64 or conda env) is on your system `PATH`, or copy the required DLLs
alongside the executable.

## Usage

```bash
./terra show ../inst/ex/elev.tif
./terra aggregate ../inst/ex/elev.tif out.tif 10 mean 1
```

## How it works

### Library

The Makefile automatically collects every `../src/*.cpp` and excludes the three
R/Rcpp interface files (`RcppExports.cpp`, `RcppFunctions.cpp`,
`RcppModule.cpp`). Everything is compiled with `-Dstandalone` so that
R-specific code paths (guarded by `#ifdef useRcpp` in the terra sources) are
skipped. The resulting object files are packed into `build/libterra.a`.

### GDAL / PROJ initialisation

R's terra initialises GDAL in `.onLoad` via `gdal_init()`. The standalone app
mirrors this in `gdal_init.cpp` (`terra_gdal_app_init`), which:

* sets `GDAL_DATA`, `GDAL_DRIVER_PATH`, and PROJ search paths (detected at
  build time by the Makefile, overridable at runtime via environment variables)
* calls `GDALAllRegister()` / `OGRRegisterAll()`
* enables PROJ network access (PROJ >= 7.1)

The paths are baked in at compile time as `TERRA_GDAL_DATA`,
`TERRA_PROJ_DATA`, and `TERRA_GDAL_PLUGINDIR` macros. At runtime the following
environment variables take precedence when set:

| Environment variable | Purpose |
|----------------------|---------|
| `GDAL_DATA` | GDAL auxiliary data (e.g. coordinate-system definitions) |
| `GDAL_DRIVER_PATH` | Directory containing GDAL plugin drivers (`.so` / `.dll`) |
| `PROJ_LIB` or `PROJ_DATA` | PROJ data directory (datum grids, init files) |

### App source files

| File | Role |
|------|------|
| `main.cpp` | Entry point — parses arguments, initialises GDAL, dispatches commands |
| `gdal_init.cpp` / `gdal_init.h` | GDAL/PROJ driver registration and configuration |
| `show.cpp` / `show.h` | Print `SpatRaster` / `SpatVector` summaries |
| `fun.h` | Helper functions (e.g. `aggregate`) |

## Linking your own program

To use `libterra.a` in your own C++ project:

```bash
# 1. Build the library
make -C cmdapp lib

# 2. Compile and link your program
g++ -std=c++17 -Dstandalone $(gdal-config --cflags) -Isrc/ \
    my_program.cpp \
    cmdapp/build/libterra.a \
    $(gdal-config --libs) $(gdal-config --dep-libs) -lgeos_c \
    -o my_program
```

Include `gdal_init.h` from this directory and call `terra_gdal_app_init()`
before constructing any `SpatRaster` from a file. See `main.cpp` for a minimal
working example.
