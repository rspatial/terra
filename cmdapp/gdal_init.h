#ifndef TERRA_GDAL_APP_INIT_H
#define TERRA_GDAL_APP_INIT_H

/* Same role as R package .gdalinit / gdal_init (RcppFunctions.cpp): register drivers,
 * set CPL options, point PROJ at search paths. Call once before using SpatRaster. */
void terra_gdal_app_init(const char *gdal_data, const char *proj_search_path);

#endif
