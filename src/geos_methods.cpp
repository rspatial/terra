#include "geos_spat.h"


SpatVector SpatVector::allerretour() {
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	SpatVector out = vect_from_geos(g, hGEOSCtxt, type());
	geos_finish(hGEOSCtxt);
	return out;
}



SpatVectorCollection SpatVector::bienvenue() {
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	SpatVectorCollection out = coll_from_geos(g, hGEOSCtxt);
	geos_finish(hGEOSCtxt);
	return out;
}


std::vector<bool> SpatVector::geos_isvalid() {
	GEOSContextHandle_t hGEOSCtxt = geos_init2();
	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<bool> out;
	out.reserve(g.size());
	for (size_t i = 0; i < g.size(); i++) {		
		char v = GEOSisValid_r(hGEOSCtxt, g[i].get());
		out.push_back(v);
	}
	geos_finish(hGEOSCtxt);
	return {out};
}

std::vector<std::string> SpatVector::geos_isvalid_msg() {
	GEOSContextHandle_t hGEOSCtxt = geos_init2();
	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<bool> ok;
	ok.reserve(g.size());
	for (size_t i = 0; i < g.size(); i++) {		
		char v = GEOSisValid_r(hGEOSCtxt, g[i].get());
		ok.push_back(v);
	}
	geos_finish(hGEOSCtxt);
	std::vector<std::string> out(ok.size(), "");
	size_t j=0;
	for (size_t i=0; i<ok.size(); i++) {
		if (!ok[i]) {
			if (j < msgs.size()) {
				out[i] = msgs[j];
				j++;
			} else {
				out[0] = "???";	
			}
		}
	}
    msgs.resize(0); 
	return out;
}
  

SpatVector SpatVector::crop(SpatExtent e) {
	SpatVector out;
	GEOSContextHandle_t hGEOSCtxt = geos_init();

	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> p;
	p.reserve(g.size());
	std::vector<unsigned> id;
	id.reserve(g.size());
	for (size_t i = 0; i < g.size(); i++) {
		GEOSGeometry* r = GEOSClipByRect_r(hGEOSCtxt, g[i].get(), e.xmin, e.ymin, e.xmax, e.ymax);
		if (r == NULL) {
			out.setError("something bad happened");
			geos_finish(hGEOSCtxt);
			return out;
		}
		if (!GEOSisEmpty_r(hGEOSCtxt, r)) {
			p.push_back(geos_ptr(r, hGEOSCtxt));	
			id.push_back(i);
		}
	}
	if (p.size() > 0) {
		SpatVectorCollection coll = coll_from_geos(p, hGEOSCtxt);
		out = coll.get(0);
		out.msg = coll.msg;
		out.df = df.subset_rows(id);
	}
	geos_finish(hGEOSCtxt);
	return out;
}


SpatVector SpatVector::voronoi(SpatVector bnd, double tolerance, int onlyEdges) {
	SpatVector out;

#ifndef HAVE350
	out.setError("GEOS 3.5 required for voronoi")
#else 

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	SpatVector a = aggregate(false);
	std::vector<GeomPtr> g = geos_geoms(&a, hGEOSCtxt);
	std::string vt = type();
	GEOSGeometry* v;
	if (bnd.size() > 0) {
		if (bnd.type() != "polygons") {
			out.setError("boundary must be polygon");
			geos_finish(hGEOSCtxt);
			return out;
		}
		std::vector<GeomPtr> ge = geos_geoms(&bnd, hGEOSCtxt);
		v = GEOSVoronoiDiagram_r(hGEOSCtxt, g[0].get(), ge[0].get(), tolerance, onlyEdges);
	} else {
		std::vector<GeomPtr> ge = geos_geoms(&bnd, hGEOSCtxt);
		v = GEOSVoronoiDiagram_r(hGEOSCtxt, g[0].get(), NULL, tolerance, onlyEdges);
	}
	if (v == NULL) {
		out.setError("GEOS exception");
		geos_finish(hGEOSCtxt);
		return(out);
	} 
	std::vector<GeomPtr> b(1);
	b[0] = geos_ptr(v, hGEOSCtxt);	
	SpatVectorCollection coll = coll_from_geos(b, hGEOSCtxt);
	geos_finish(hGEOSCtxt);
	out = coll.get(0);
	out.msg = coll.msg;
	if (!out.hasError()) {
		out = out.disaggregate();
		// associate with attributes
	}
#endif
	return out;
}  

SpatVector SpatVector::buffer2(double dist, unsigned nQuadSegs, unsigned capstyle) {

	std::string vt = type();
	if ((vt == "points") && (dist <= 0)) {
		dist = -dist;
	}

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	SpatVector out;
//	SpatVector f = remove_holes();

	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> b(size());
	for (size_t i = 0; i < g.size(); i++) {
		GEOSGeometry* pt = GEOSBuffer_r(hGEOSCtxt, g[i].get(), dist, nQuadSegs);
		if (pt == NULL) {
			out.setError("GEOS exception");
			geos_finish(hGEOSCtxt);
			return(out);
		} 
		b[i] = geos_ptr(pt, hGEOSCtxt);	
	}
	SpatVectorCollection coll = coll_from_geos(b, hGEOSCtxt);

//	out = spat_from_geom(hGEOSCtxt, g, "points");
	geos_finish(hGEOSCtxt);
	return coll.get(0);	
}


SpatVector SpatVector::intersect(SpatVector v) {

	SpatVector out;

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> y = geos_geoms(&v, hGEOSCtxt);	
	std::vector<GeomPtr> result;
	std::vector<std::vector<unsigned>> atts(2);
	size_t nx = size();
	size_t ny = v.size();
	
		
	for (size_t i = 0; i < nx; i++) {
		for (size_t j = 0; j < ny; j++) {
			GEOSGeometry* geom = GEOSIntersection_r(hGEOSCtxt, x[i].get(), y[j].get());
			if (geom == NULL) {
				out.setError("GEOS exception");
				geos_finish(hGEOSCtxt);
				return(out);
			} 
			if (!GEOSisEmpty_r(hGEOSCtxt, geom)) {
				result.push_back(geos_ptr(geom, hGEOSCtxt));
				atts[0].push_back(i);
				atts[1].push_back(j);
			}
		}
	}

	SpatVectorCollection coll = coll_from_geos(result, hGEOSCtxt);
	// deal with attributes
	geos_finish(hGEOSCtxt);
	return coll.get(0);	
}


std::vector<bool> SpatVector::intersects(SpatVector v) {
	std::vector<bool> out;
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> y = geos_geoms(&v, hGEOSCtxt);	
	size_t nx = size();
	size_t ny = v.size();
	out.reserve(nx*ny);
	for (size_t i = 0; i < nx; i++) {
		for (size_t j = 0; j < ny; j++) {
			if (GEOSIntersects_r(hGEOSCtxt, x[i].get(), y[j].get())) {
				out.push_back(true);
			} else {
				out.push_back(false);
			}
		}
	}
	geos_finish(hGEOSCtxt);
	return out;	
}



SpatVector SpatVector::centroid() {

	SpatVector out;

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> b(size());
	for (size_t i = 0; i < g.size(); i++) {		
		GEOSGeometry* pt = GEOSGetCentroid_r(hGEOSCtxt, g[i].get());
		if (pt == NULL) {
			out.setError("NULL geom");
			geos_finish(hGEOSCtxt);
			return out;
		}
		b[i] = geos_ptr(pt, hGEOSCtxt);	
	}
	out = vect_from_geos(b, hGEOSCtxt, "points");
	geos_finish(hGEOSCtxt);

	out.df = df;
	return out;	
}

SpatVector SpatVector::unaryunion() {
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> gout(g.size());
	SpatVector out;
	for (size_t i = 0; i < g.size(); i++) {		
		GEOSGeometry* u = GEOSUnaryUnion_r(hGEOSCtxt, g[i].get());
		if (u == NULL) {
			out.setError("NULL geom");
			geos_finish(hGEOSCtxt);
			return out;
		}
		gout[i] = geos_ptr(u, hGEOSCtxt);
	}
	SpatVectorCollection coll = coll_from_geos(gout, hGEOSCtxt);
	geos_finish(hGEOSCtxt);
	out = coll.get(0);
	out.srs = srs;
	return out;
}

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

