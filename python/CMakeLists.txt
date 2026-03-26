cmake_minimum_required(VERSION 3.18)
project(pyterra CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ── pybind11 ────────────────────────────────────────────────────────────────
# Installed via:  pip install pybind11
# or via cmake:   find_package(pybind11 CONFIG REQUIRED)
find_package(pybind11 CONFIG REQUIRED)


# ── GDAL ────────────────────────────────────────────────────────────────────
find_package(GDAL REQUIRED)


# ── GEOS ────────────────────────────────────────────────────────────────────
# Try cmake config first (GEOS >= 3.9 ships one), fall back to pkg-config.
find_package(GEOS QUIET CONFIG)
if (NOT GEOS_FOUND)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(GEOS REQUIRED geos)
endif()


# ── PROJ ────────────────────────────────────────────────────────────────────
find_package(PROJ CONFIG QUIET)
if (NOT PROJ_FOUND)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(PROJ REQUIRED proj)
endif()


# ── terra C++ sources (exclude all R/Rcpp-specific files) ───────────────────
set(TERRA_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../src")

file(GLOB TERRA_SOURCES "${TERRA_SRC_DIR}/*.cpp")

# Remove the three Rcpp-dependent files; everything else is pure C++.
list(FILTER TERRA_SOURCES EXCLUDE REGEX ".*(RcppModule|RcppFunctions|RcppExports)\\.cpp$")


# ── GEOS version feature flags ───────────────────────────────────────────────
# terra uses #ifdef GEOS350 / GEOS380 / GEOS390 guards internally.
if (GEOS_FOUND)
    # Normalise version string to a comparable integer form
    if (DEFINED GEOS_VERSION)
        set(_geos_ver "${GEOS_VERSION}")
    elseif (DEFINED GEOS_VERSION_STRING)
        set(_geos_ver "${GEOS_VERSION_STRING}")
    else()
        set(_geos_ver "0.0.0")
    endif()

    if (_geos_ver VERSION_GREATER_EQUAL "3.5.0")
        set(GEOS_CXXFLAGS "-DGEOS350")
    endif()
    if (_geos_ver VERSION_GREATER_EQUAL "3.8.0")
        string(APPEND GEOS_CXXFLAGS " -DGEOS380")
    endif()
    if (_geos_ver VERSION_GREATER_EQUAL "3.9.0")
        string(APPEND GEOS_CXXFLAGS " -DGEOS390")
    endif()
    if (_geos_ver VERSION_GREATER_EQUAL "3.11.0")
        string(APPEND GEOS_CXXFLAGS " -DGEOS3110")
    endif()
    if (_geos_ver VERSION_GREATER_EQUAL "3.12.0")
        string(APPEND GEOS_CXXFLAGS " -DGEOS3120")
    endif()

    separate_arguments(GEOS_CXXFLAGS_LIST UNIX_COMMAND "${GEOS_CXXFLAGS}")
endif()


# ── PROJ version flag ────────────────────────────────────────────────────────
if (PROJ_FOUND)
    if (DEFINED PROJ_VERSION)
        set(_proj_ver "${PROJ_VERSION}")
    elseif (DEFINED PROJ_VERSION_STRING)
        set(_proj_ver "${PROJ_VERSION_STRING}")
    else()
        set(_proj_ver "0.0.0")
    endif()
    if (_proj_ver VERSION_GREATER_EQUAL "6.0")
        set(PROJ_CXXFLAGS "-DHAVE_PROJ_H")
    endif()
    separate_arguments(PROJ_CXXFLAGS_LIST UNIX_COMMAND "${PROJ_CXXFLAGS}")
endif()


# ── TBB (optional) ──────────────────────────────────────────────────────────
find_package(TBB QUIET)
if (TBB_FOUND)
    set(TBB_LINK_LIBS TBB::tbb)
    set(TBB_COMPILE_DEFS HAVE_TBB)
endif()


# ── Build the extension ──────────────────────────────────────────────────────
pybind11_add_module(_pyterra
    src/pyterra_module.cpp
    ${TERRA_SOURCES}
)

target_include_directories(_pyterra PRIVATE
    ${TERRA_SRC_DIR}
    ${GDAL_INCLUDE_DIRS}
)

# Absorb GEOS / PROJ include directories whether they came from cmake config
# packages (targets) or from pkg-config (variables).
if (TARGET GEOS::geos_c)
    target_link_libraries(_pyterra PRIVATE GEOS::geos_c)
elseif (DEFINED GEOS_INCLUDE_DIRS)
    target_include_directories(_pyterra PRIVATE ${GEOS_INCLUDE_DIRS})
    target_link_directories(_pyterra PRIVATE ${GEOS_LIBRARY_DIRS})
    target_link_libraries(_pyterra PRIVATE ${GEOS_LIBRARIES})
endif()

if (TARGET PROJ::proj)
    target_link_libraries(_pyterra PRIVATE PROJ::proj)
elseif (DEFINED PROJ_INCLUDE_DIRS)
    target_include_directories(_pyterra PRIVATE ${PROJ_INCLUDE_DIRS})
    target_link_directories(_pyterra PRIVATE ${PROJ_LIBRARY_DIRS})
    target_link_libraries(_pyterra PRIVATE ${PROJ_LIBRARIES})
endif()

target_link_libraries(_pyterra PRIVATE
    GDAL::GDAL
    ${TBB_LINK_LIBS}
)

target_compile_definitions(_pyterra PRIVATE
    standalone          # suppresses useRcpp; keeps useGDAL active
    ${GEOS_CXXFLAGS_LIST}
    ${PROJ_CXXFLAGS_LIST}
    ${TBB_COMPILE_DEFS}
)

# Place the built .so / .pyd directly inside the Python package directory so
# an editable install (pip install -e .) works without further configuration.
set_target_properties(_pyterra PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/pyterra"
)
