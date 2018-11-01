#include "spatraster.h"
#include "util.h"
using namespace std;


bool SpatRaster::compare_geom(SpatRaster x, bool lyrs, bool crs) {
	bool e1 = is_equal(x.extent.xmax, extent.xmax, 1);
	bool e2 = is_equal(x.extent.xmin, extent.xmin, 1);
	bool e3 = is_equal(x.extent.ymax, extent.ymax, 1);
	bool e4 = is_equal(x.extent.ymin, extent.ymin, 1);
	bool eOK = (e1 && e2 && e3 && e4);
	bool rcOK = (nrow == x.nrow) && (ncol == x.ncol);
	bool lyrOK = true;
	if (lyrs) {
		lyrOK = nlyr() == x.nlyr();
	} 
	bool crsOK = true;
	if (crs) {
		crsOK = getCRS() == x.getCRS();
	}
	return (rcOK && eOK && lyrOK && crsOK);
}

