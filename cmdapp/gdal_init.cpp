#include <cstdio>
#include <vector>

#include "gdal_priv.h"
#include "ogr_api.h"
#include "ogr_spatialref.h"
#include "cpl_conv.h"
#include "cpl_error.h"
#include "proj.h"

#include "gdal_init.h"

static void terra_gdal_error_handler(CPLErr eErrClass, int err_no, const char *msg) {
	(void)err_no;
	if (eErrClass >= CE_Failure && msg != nullptr)
		std::fprintf(stderr, "%s\n", msg);
}

void terra_gdal_app_init(const char *gdal_data, const char *proj_search_path,
                         const char *gdal_plugindir) {
	CPLSetErrorHandler(terra_gdal_error_handler);

	if (gdal_data && gdal_data[0])
		CPLSetConfigOption("GDAL_DATA", gdal_data);
	if (gdal_plugindir && gdal_plugindir[0])
		CPLSetConfigOption("GDAL_DRIVER_PATH", gdal_plugindir);

	GDALAllRegister();
	OGRRegisterAll();

	CPLSetConfigOption("GDAL_MAX_BAND_COUNT", "9999999");
	CPLSetConfigOption("OGR_CT_FORCE_TRADITIONAL_GIS_ORDER", "YES");
	CPLSetConfigOption("CPL_VSIL_USE_TEMP_FILE_FOR_RANDOM_WRITE", "YES");

#if GDAL_VERSION_NUM >= 3000000
	if (proj_search_path != nullptr && proj_search_path[0] != '\0') {
		std::vector<char *> cpaths(2);
		cpaths[0] = const_cast<char *>(proj_search_path);
		cpaths[1] = nullptr;
		OSRSetPROJSearchPaths(cpaths.data());
	}
#endif

#if PROJ_VERSION_MAJOR > 7 || (PROJ_VERSION_MAJOR == 7 && PROJ_VERSION_MINOR >= 1)
#ifndef __EMSCRIPTEN__
	proj_context_set_enable_network(PJ_DEFAULT_CTX, 1);
#endif
#endif
}
