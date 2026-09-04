// Copyright (c) 2018-2026  Robert J. Hijmans
//
// Compatibility shims for older GDAL headers.
// Include after gdal_priv.h / gdal.h / gdal_version.h.

#ifndef TERRA_GDAL_COMPAT_H
#define TERRA_GDAL_COMPAT_H

#include "gdal_version.h"

// CSLConstList was added in GDAL 2.3 (const char * const *)
#if GDAL_VERSION_NUM < GDAL_COMPUTE_VERSION(2, 3, 0)
#ifndef CSLConstList
typedef char **CSLConstList;
#endif
#endif

#endif
