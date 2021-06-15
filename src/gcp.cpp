#include "gdalwarper.h"
#include "gdal_priv.h"
#include "cpl_string.h"
#include <gdal.h>
#include "spatRaster.h"


SpatRaster SpatRaster::applyGCP(std::vector<double> fx, std::vector<double> fy, std::vector<double> tx, std::vector<double> ty) {

	SpatRaster out;
	std::vector<double> cls = cellFromXY(fx, fy);
	std::vector<std::vector<int_64>> rc = rowColFromCell(cls);
	
    GDAL_GCP *gcps = NULL;
    gcps = (GDAL_GCP *) CPLRealloc (gcps, (fx.size()) * sizeof(GDAL_GCP));
	GDALInitGCPs(fx.size(), gcps);
    for (size_t i = 0; i < fx.size(); i++){
      gcps[i].dfGCPPixel = rc[1][i];
      gcps[i].dfGCPLine = rc[0][i];
      gcps[i].dfGCPX = tx[i];
      gcps[i].dfGCPY = ty[i];
      gcps[i].dfGCPZ = (float) 0.0;
    }

	SpatOptions opt;
    GDALDatasetH hSrcDS; //hDstDS, 

	if (!open_gdal(hSrcDS, 0, opt)) {
		out.setError("bad");
		return out;
	}
	std::string srccrs = getSRS("wkt");
    const char *projection = srccrs.c_str();
    GDALSetGCPs(hSrcDS, fx.size(), gcps, projection); 

	//if (!get_output_bounds(hSrcDS, srccrs, srccrs, out)) {
	//	GDALClose( hSrcDS );
	//	return out;
	//}
	
	return out;

}

