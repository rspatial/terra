#include "geos_spat.h"

SpatVector SpatVector::buffer2(double dist, unsigned nQuadSegs, unsigned capstyle) {

	GEOSContextHandle_t hGEOSCtxt = CPL_geos_init();
	SpatVector out;
	
	std::vector<GeomPtr> g = geom_from_spat(this, hGEOSCtxt);

	std::vector<GeomPtr> b(size());

	std::string vt = type();
	if ((vt == "points") && (dist <= 0)) {
		dist = -1 * dist;
	}
	for (size_t i = 0; i < g.size(); i++) {
		GEOSGeometry* pt = GEOSBuffer_r(hGEOSCtxt, g[i].get(), dist, nQuadSegs);
		if (pt == NULL) {
			out.setError("GEOS exception");
			CPL_geos_finish(hGEOSCtxt);
			return(out);
		} 
		b[i] = geos_ptr(pt, hGEOSCtxt);	
	}

	out = spat_from_geom(hGEOSCtxt, b, "polygons");

//	out = spat_from_geom(hGEOSCtxt, g, "points");
	CPL_geos_finish(hGEOSCtxt);
	return out;	
}
