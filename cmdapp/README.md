# cmdapp — standalone terra C++ example

This directory builds the [terra](https://github.com/rspatial/terra) C++ spatial
library as a **standalone static library** (`libterra_standalone.a`) and links a
small command-line demo against it — no R or Rcpp required.

## Prerequisites

* A C++17 compiler (`g++` or `clang++`)
* [GDAL](https://gdal.org/) development files and `gdal-config` on `PATH`
* [PROJ](https://proj.org/) development files and `pkg-config proj`
* [GEOS](https://libgeos.org/) (`-lgeos_c`)
* GNU Make

On Debian / Ubuntu:

```bash
sudo apt install build-essential gdal-bin libgdal-dev libproj-dev libgeos-dev
```

## Building

```bash
cd cmdapp
make            # build library + ./terra executable
make -j$(nproc) # parallel build
```

Build artefacts go into `build/`:

| Target | File | Description |
|--------|------|-------------|
| `make lib` | `build/libterra_standalone.a` | Static library (all `src/*.cpp` except Rcpp glue) |
| `make` or `make all` | `./terra` | Demo executable linked against the static lib |
| `make clean` | | Remove all build artefacts |

There is also a thin shell wrapper `compile.sh` that `cd`s here and calls `make`:

```bash
./compile.sh -j4        # same as: cd cmdapp && make -j4
./compile.sh clean      # same as: cd cmdapp && make clean
```

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
skipped. The resulting object files are packed into
`build/libterra_standalone.a`.

### GDAL / PROJ initialisation

R's terra initialises GDAL in `.onLoad` via `gdal_init()`. The standalone app
mirrors this in `gdal_init.cpp` (`terra_gdal_app_init`), which:

* sets `GDAL_DATA`, `GDAL_DRIVER_PATH`, and PROJ search paths (detected at
  build time by the Makefile, overridable at runtime via environment variables)
* calls `GDALAllRegister()` / `OGRRegisterAll()`
* enables PROJ network access

The paths are baked in at compile time as `TERRA_GDAL_DATA`,
`TERRA_PROJ_DATA`, and `TERRA_GDAL_PLUGINDIR` macros. At runtime the
environment variables `GDAL_DATA`, `PROJ_LIB` / `PROJ_DATA`, and
`GDAL_DRIVER_PATH` take precedence when set.

### App source files

| File | Role |
|------|------|
| `main.cpp` | Entry point — parses arguments, initialises GDAL, dispatches commands |
| `gdal_init.cpp` / `gdal_init.h` | GDAL/PROJ driver registration and config |
| `show.cpp` / `show.h` | Print `SpatRaster` / `SpatVector` summaries |
| `fun.h` | Helper functions (e.g. `aggregate`) |

## Linking your own program

To use `libterra_standalone.a` in your own C++ project:

```bash
# 1. Build the library
make -C cmdapp lib

# 2. Compile and link
g++ -std=c++17 -Dstandalone $(gdal-config --cflags) -Isrc/ \
    my_program.cpp \
    cmdapp/build/libterra_standalone.a \
    $(gdal-config --libs) $(gdal-config --dep-libs) -lgeos_c \
    -o my_program
```

Remember to call `terra_gdal_app_init()` (or at minimum `GDALAllRegister()`)
before constructing any `SpatRaster` from a file.
