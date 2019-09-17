
#include "geos_spat.h"

/*
bool geos_buffer(GEOSContextHandle_t hGEOSCtxt, std::vector<GeomPtr> &g, double dist, unsigned nQuadSegs) {
	std::vector<GeomPtr> g(size());
	for (size_t i = 0; i < g.size(); i++) {
		GEOSGeometry* pt = GEOSBuffer_r(hGEOSCtxt, g[i].get(), dist, nQuadSegs);
		if (pt == NULL) {
			return false;
		} 
		g[i] = geos_ptr(pt, hGEOSCtxt);	
	}
	return true;
}
*/

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




SpatVector SpatVector::intersect(SpatVector v) {

	SpatVector out;

	GEOSContextHandle_t hGEOSCtxt = CPL_geos_init();
	std::vector<GeomPtr> x = geom_from_spat(this, hGEOSCtxt);
	std::vector<GeomPtr> y = geom_from_spat(&v, hGEOSCtxt);	
	std::vector<GeomPtr> result;
	std::vector<std::vector<unsigned>> atts(2);
	size_t nx = size();
	size_t ny = v.size();
	
		
	for (size_t i = 0; i < nx; i++) {
		for (size_t j = 0; j < ny; j++) {
			GEOSGeometry* geom = GEOSIntersection_r(hGEOSCtxt, x[i].get(), y[j].get());
			if (geom == NULL) {
				out.setError("GEOS exception");
				CPL_geos_finish(hGEOSCtxt);
				return(out);
			} 
			unsigned check = 0; //GEOSisEmpty_r(hGEOSCtxt, geom);
			if (check > 1) {
				out.setError("GEOS exception");
				CPL_geos_finish(hGEOSCtxt);
				return(out);			
			} else if (check == 1) {
				result.push_back(geos_ptr(geom, hGEOSCtxt));
				atts[0].push_back(i);
				atts[1].push_back(j);
			}
		}
	}

	// this should happen in the loop. And different types are possible!
	out = spat_from_geom(hGEOSCtxt, result, type());
	// deal with attributes
	CPL_geos_finish(hGEOSCtxt);
	return out;	
}



