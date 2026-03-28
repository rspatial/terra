#ifndef TERRA_GDAL_INIT_H
#define TERRA_GDAL_INIT_H

/* Same role as R package .gdalinit / gdal_init (RcppFunctions.cpp): register drivers,
 * set CPL options, point PROJ at search paths. Call once before using SpatRaster.
 * gdal_plugindir: path to GDAL plugin .so files (needed on distros that ship
 * drivers like GTiff as plugins). Set BEFORE GDALAllRegister(). */
void terra_gdal_app_init(const char *gdal_data, const char *proj_search_path,
                         const char *gdal_plugindir);

#endif
