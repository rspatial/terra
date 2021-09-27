#include <numeric>
#include "geos_spat.h"
#include "distance.h"
#include "recycle.h"
#include "string_utils.h"

	
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


std::vector<std::string> SpatVector::wkt() {
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<std::string> out;
	out.reserve(g.size());
	char * wkt;
	for (size_t i = 0; i < g.size(); i++) {	
		wkt = GEOSGeomToWKT_r(hGEOSCtxt, g[i].get());
		out.push_back(wkt);
	}
	geos_finish(hGEOSCtxt);
	return out;
}



std::vector<std::string> SpatVector::wkb() {
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<std::string> out;
	out.reserve(g.size());
	size_t len = 0;
	for (size_t i = 0; i < g.size(); i++) {	
		unsigned char *wkb = GEOSGeomToWKB_buf_r(hGEOSCtxt, g[i].get(), &len);
		std::string s( reinterpret_cast<char const*>(wkb), len) ;
		out.push_back(s);
		free(wkb);
	}
	geos_finish(hGEOSCtxt);
	return out;
}

std::vector<std::string> SpatVector::hex() {
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<std::string> out;
	out.reserve(g.size());
	size_t len = 0;
	for (size_t i = 0; i < g.size(); i++) {	
		unsigned char *hex = GEOSGeomToHEX_buf_r(hGEOSCtxt, g[i].get(), &len);
		std::string s( reinterpret_cast<char const*>(hex), len) ;		
		out.push_back(s);
		free(hex);
	}
	geos_finish(hGEOSCtxt);
	return out;
}



SpatVector SpatVector::from_hex(std::vector<std::string> x, std::string srs) {
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	size_t n = x.size();
	std::vector<GeomPtr> p;
	p.resize(n);
	for (size_t i = 0; i < n; i++) {
		const char* cstr = x[i].c_str();	
		size_t len = strlen(cstr);
		const unsigned char *hex = (const unsigned char *) cstr;
		GEOSGeometry* r = GEOSGeomFromHEX_buf_r(hGEOSCtxt, hex, len);
		p[i] = geos_ptr(r, hGEOSCtxt);
	}
	SpatVectorCollection coll = coll_from_geos(p, hGEOSCtxt);
	geos_finish(hGEOSCtxt);
	SpatVector out = coll.get(0);
	out.setSRS(srs);
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
	std::vector<std::string> out;
	out.reserve(2 * g.size());
	for (size_t i = 0; i < g.size(); i++) {	
		char v = GEOSisValid_r(hGEOSCtxt, g[i].get());
		std::string valid = {v};
		out.push_back(valid);
		if (!v) {
			char *r = GEOSisValidReason_r(hGEOSCtxt, g[i].get());
			std::string reason = r;
			free(r);
			out.push_back(reason);
		} else {
			out.push_back("");
		}
	}
	geos_finish(hGEOSCtxt);
	return {out};
}


SpatVector SpatVector::crop(SpatExtent e) {
	SpatVector out;
	out.srs = srs;
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
		} else {
			GEOSGeom_destroy_r(hGEOSCtxt, r);
		}
	}
	if (p.size() > 0) {
		SpatVectorCollection coll = coll_from_geos(p, hGEOSCtxt, id);
		out = coll.get(0);
		out.df = df.subset_rows(out.df.iv[0]);
		out.srs = srs;
	}
	geos_finish(hGEOSCtxt);
	return out;
}


SpatVector SpatVector::make_nodes() {

	SpatVector out;
	GEOSContextHandle_t hGEOSCtxt = geos_init();

	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> p;
	p.reserve(g.size());
	for (size_t i = 0; i < g.size(); i++) {
		GEOSGeometry* r = GEOSNode_r(hGEOSCtxt, g[i].get());
		if (r == NULL) {
			out.setError("something bad happened");
			geos_finish(hGEOSCtxt);
			return out;
		}
		if (!GEOSisEmpty_r(hGEOSCtxt, r)) {
			p.push_back(geos_ptr(r, hGEOSCtxt));
		} else {
			GEOSGeom_destroy_r(hGEOSCtxt, r);
		}
	}
	if (p.size() > 0) {
		SpatVectorCollection coll = coll_from_geos(p, hGEOSCtxt);
		out = coll.get(0);
		out.df = df;
	}
	geos_finish(hGEOSCtxt);
	out.srs = srs;
	return out;
}




SpatVector SpatVector::boundary() {

	SpatVector out;
	GEOSContextHandle_t hGEOSCtxt = geos_init();

	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> p;
	p.reserve(g.size());
	for (size_t i = 0; i < g.size(); i++) {
		GEOSGeometry* r = GEOSBoundary_r(hGEOSCtxt, g[i].get());
		if (r == NULL) {
			out.setError("something bad happened");
			geos_finish(hGEOSCtxt);
			return out;
		}
		if (!GEOSisEmpty_r(hGEOSCtxt, r)) {
			p.push_back(geos_ptr(r, hGEOSCtxt));
		} else {
			GEOSGeom_destroy_r(hGEOSCtxt, r);
		}
	}
	if (p.size() > 0) {
		SpatVectorCollection coll = coll_from_geos(p, hGEOSCtxt);
		out = coll.get(0);
		out.df = df;
	}
	geos_finish(hGEOSCtxt);
	out.srs = srs;
	return out;
}


SpatVector SpatVector::normalize() {

	SpatVector out;
	GEOSContextHandle_t hGEOSCtxt = geos_init();

	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> p;
	p.reserve(g.size());
	for (size_t i = 0; i < g.size(); i++) {
		GEOSGeometry* r = g[i].get();
		if (GEOSNormalize_r(hGEOSCtxt, r)) {
			g[i] = geos_ptr(r, hGEOSCtxt);
		} else {
			GEOSGeom_destroy_r(hGEOSCtxt, r);
		}
	}
	out = vect_from_geos(g, hGEOSCtxt, type());
	geos_finish(hGEOSCtxt);
	out.df = df;
	out.srs = srs;
	return out;
}



SpatVector SpatVector::line_merge() {
	SpatVector out;
	GEOSContextHandle_t hGEOSCtxt = geos_init();

	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> p;
	p.reserve(g.size());
	for (size_t i = 0; i < g.size(); i++) {
		GEOSGeometry* r = GEOSLineMerge_r(hGEOSCtxt, g[i].get());
		if (r == NULL) {
			out.setError("something bad happened");
			geos_finish(hGEOSCtxt);
			return out;
		}
		if (!GEOSisEmpty_r(hGEOSCtxt, r)) {
			p.push_back(geos_ptr(r, hGEOSCtxt));
		} else {
			GEOSGeom_destroy_r(hGEOSCtxt, r);
		}
	}
	if (p.size() > 0) {
		SpatVectorCollection coll = coll_from_geos(p, hGEOSCtxt);
		out = coll.get(0);
		out.df = df;
	}
	geos_finish(hGEOSCtxt);
	out.srs = srs;
	return out;
}



SpatVector SpatVector::simplify(double tolerance, bool preserveTopology) {
	SpatVector out;
	GEOSContextHandle_t hGEOSCtxt = geos_init();

	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> p;
	p.reserve(g.size());
	for (size_t i = 0; i < g.size(); i++) {
		GEOSGeometry* r;
		if (preserveTopology) {
			r = GEOSTopologyPreserveSimplify_r(hGEOSCtxt, g[i].get(), tolerance);			
		} else {
			r = GEOSSimplify_r(hGEOSCtxt, g[i].get(), tolerance);
		}
		if (r == NULL) {
			out.setError("something bad happened");
			geos_finish(hGEOSCtxt);
			return out;
		}
		if (!GEOSisEmpty_r(hGEOSCtxt, r)) {
			p.push_back(geos_ptr(r, hGEOSCtxt));
		} else {
			GEOSGeom_destroy_r(hGEOSCtxt, r);
		}
	}
	if (p.size() > 0) {
		SpatVectorCollection coll = coll_from_geos(p, hGEOSCtxt);
		out = coll.get(0);
		out.df = df;
	}
	geos_finish(hGEOSCtxt);
	out.srs = srs;
	return out;
}



SpatVector SpatVector::shared_paths() {
	
	if (type() == "polygons") {
		SpatVector x = as_lines();
		return x.shared_paths();
	}
	
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);

	size_t s = size();
	std::vector<long> id1, id2;
	std::vector<GeomPtr> p;
	for (size_t i=0; i<(s-1); i++) {
		for (size_t j=(i+1); j<s; j++) {
			GEOSGeometry* r = GEOSSharedPaths_r(hGEOSCtxt, x[i].get(), x[j].get());
			if (r != NULL) {
				if (!GEOSisEmpty_r(hGEOSCtxt, r)) {
					p.push_back(geos_ptr(r, hGEOSCtxt));
					id1.push_back(i+1);
					id2.push_back(j+1);				
				} else {
					GEOSGeom_destroy_r(hGEOSCtxt, r);
				}
			}
		}
	}
	
	SpatVector out;
	if (p.size() > 0) {
		SpatVectorCollection coll = coll_from_geos(p, hGEOSCtxt, std::vector<unsigned>(), false, false);
		out = coll.get(0);
		out = out.line_merge();
	}
	geos_finish(hGEOSCtxt);
	out.srs = srs;
	out.df.add_column(id1, "id1");
	out.df.add_column(id2, "id2");
	return out;
}


/*
SpatVector SpatVector::split_polygons(SpatVector lns) {
		SpatGeom glns;
		glns.gtype = lines;
		glns.setPart(SpatPart(x, y), 0);
		std::vector<double> xln = {180, 180};
		std::vector<double> yln = {-91, 91};
		glns.setPart(SpatPart(xln, yln), 1);
		SpatVector v;
		v.addGeom(glns);	
		v = v.line_merge();
		v = v.aggregate(false);
		v = v.polygonize();
		g = v.geoms[0];		
*/

SpatVector SpatVector::polygonize() {

	SpatVector out;
	out.srs = srs;
	GEOSContextHandle_t hGEOSCtxt = geos_init();

	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> p;
	p.reserve(g.size());
	size_t ngeoms = 1;
	for (size_t i = 0; i < g.size(); i++) {
		const GEOSGeometry* gi = g[i].get();		
		GEOSGeometry* r = GEOSPolygonize_r(hGEOSCtxt, &gi, ngeoms);
		if (r == NULL) {
			out.setError("something bad happened");
			geos_finish(hGEOSCtxt);
			return out;
		}
		if (!GEOSisEmpty_r(hGEOSCtxt, r)) {
			p.push_back(geos_ptr(r, hGEOSCtxt));
		} else {
			GEOSGeom_destroy_r(hGEOSCtxt, r);
		}
	}
	if (p.size() > 0) {
		SpatVectorCollection coll = coll_from_geos(p, hGEOSCtxt);
		out = coll.get(0);
		out.srs = srs;
		out.df = df;
	}
	geos_finish(hGEOSCtxt);
	return out;
}



SpatVector SpatVector::snap(double tolerance) {
	
	tolerance = std::max(0.0, tolerance);
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);

	size_t s = size();
	std::vector<long> id;
	std::vector<GeomPtr> p;
	for (size_t i=0; i<(s-1); i++) {
		GEOSGeometry* r = x[i].get();
		for (size_t j=(i+1); j<s; j++) {
			r = GEOSSnap_r(hGEOSCtxt, r, x[j].get(), tolerance);
		}
		if (r != NULL) {
			if (!GEOSisEmpty_r(hGEOSCtxt, r)) {
				x[i] = geos_ptr(r, hGEOSCtxt);
				id.push_back(i+1);
			} else {
				GEOSGeom_destroy_r(hGEOSCtxt, r);
			}
		}
	}
	SpatVector out;
	if (p.size() > 0) {
		SpatVectorCollection coll = coll_from_geos(p, hGEOSCtxt, std::vector<unsigned>(), false, false);
		out = coll.get(0);
	}
	geos_finish(hGEOSCtxt);
	out.geoms.push_back(geoms[s-1]);
	id.push_back(s);
	out.srs = srs;
	if (out.size() == size()) {
		out.df = df;
	} else {
		out.df.add_column(id, "id");
	}
	return out;
}



//GEOSPolygonizer_getCutEdges_r(GEOSContextHandle_t extHandle, const Geometry * const * g, unsigned int ngeoms)

//Geometry * GEOSPolygonize_full_r(GEOSContextHandle_t extHandle, const Geometry* g, Geometry** cuts, Geometry** dangles, Geometry** invalid)
 



SpatVector SpatVector::crop(SpatVector v) {

	SpatVector out;
	out.srs = srs;

	GEOSContextHandle_t hGEOSCtxt = geos_init();

	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
//	if ((type() != "polygons") & (type() != "mutlipolygons")) {
	if ((type() != "polygons")) {
		v = v.hull("convex");
	} else {
		v = v.aggregate(false);
	}
	std::vector<GeomPtr> y = geos_geoms(&v, hGEOSCtxt);
	std::vector<GeomPtr> result;
	std::vector<unsigned> ids;
	size_t nx = size();
	ids.reserve(nx);
	
	for (size_t i = 0; i < nx; i++) {
		GEOSGeometry* geom = GEOSIntersection_r(hGEOSCtxt, x[i].get(), y[0].get());
		if (geom == NULL) {
			out.setError("GEOS exception");
			geos_finish(hGEOSCtxt);
			return(out);
		} 
		if (!GEOSisEmpty_r(hGEOSCtxt, geom)) {
			result.push_back(geos_ptr(geom, hGEOSCtxt));
			ids.push_back(i);
		} else {
			GEOSGeom_destroy_r(hGEOSCtxt, geom);
		}
	}

//	SpatVectorCollection coll = coll_from_geos(result, hGEOSCtxt);

	if (result.size() > 0) {
//		SpatVectorCollection coll = coll_from_geos(result, hGEOSCtxt);
		SpatVectorCollection coll = coll_from_geos(result, hGEOSCtxt, ids);
		out = coll.get(0);
		out.df = df.subset_rows(out.df.iv[0]);
		out.srs = srs;
	} 
	geos_finish(hGEOSCtxt);
	return out;
}



SpatVector SpatVector::hull(std::string htype, std::string by) {

	SpatVector out;
	if (by != "") {
		SpatVector tmp = aggregate(by, false);
		if (tmp.hasError()) {
			return tmp;	
		}
		for (size_t i=0; i<tmp.size(); i++) {
			SpatVector x = tmp.subset_rows(i);
			x = x.hull(htype, "");
			if (x.hasError()) {
				return x;
			}
			if ((x.geoms.size() > 0) && (x.geoms[0].gtype == polygons)) {
				out.addGeom(x.geoms[0]);
			} else {
				SpatGeom g;
				g.gtype = polygons;
				out.addGeom(g);
			}
		}
		out.df = tmp.df;
		out.srs = out.srs;
		return out;
	}
	
	
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	SpatVector a = aggregate(false);
	std::vector<GeomPtr> g = geos_geoms(&a, hGEOSCtxt);
	//std::string vt = type();
	GEOSGeometry* h;
	if (htype == "convex") {
		h = GEOSConvexHull_r(hGEOSCtxt, g[0].get());
	} else {
		h = GEOSMinimumRotatedRectangle_r(hGEOSCtxt, g[0].get());
	}
	std::vector<GeomPtr> b(1);
	b[0] = geos_ptr(h, hGEOSCtxt);
	SpatVectorCollection coll = coll_from_geos(b, hGEOSCtxt);
	geos_finish(hGEOSCtxt);
	out = coll.get(0);
	out.srs = srs;
	return out;
}



SpatVector SpatVector::voronoi(SpatVector bnd, double tolerance, int onlyEdges) {
	SpatVector out;

#ifndef HAVE350
	out.setError("GEOS 3.5 required for voronoi")
	return out;
#else 

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	SpatVector a = aggregate(false);
	std::vector<GeomPtr> g = geos_geoms(&a, hGEOSCtxt);
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
	out.srs = srs;
	if (!out.hasError()) {
		out = out.disaggregate();
		if (bnd.size() > 0) {
			SpatDataFrame empty;
			bnd.df = empty;
			out = out.intersect(bnd);
		}
		if ((type() == "points") && (!onlyEdges)) {
			std::vector<int> atts = out.relateFirst(*this, "intersects");
			std::vector<unsigned> a;
			a.reserve(atts.size());
			for (size_t i=0; i<atts.size(); i++) {
				if (atts[i] >=0) a.push_back(atts[i]); 
			}
			if (a.size() == out.size()) {
				out.df = df.subset_rows(a);
			}
		}
	}
#endif
	return out;
}  




SpatVector SpatVector::delauny(double tolerance, int onlyEdges) {
	SpatVector out;

#ifndef HAVE350
	out.setError("GEOS 3.5 required for delauny")
	return out;
#else 

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	SpatVector a = aggregate(false);
	std::vector<GeomPtr> g = geos_geoms(&a, hGEOSCtxt);
	GEOSGeometry* v = GEOSDelaunayTriangulation_r(hGEOSCtxt, g[0].get(), tolerance, onlyEdges);
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
	out.srs = srs;
	if (!out.hasError()) {
		out = out.disaggregate();
		// associate with attributes
	}
#endif
	return out;
}


SpatGeom hullify(SpatVector b, bool ispoly) {
	if (b.nrow() == 1) return b.geoms[0];
	if (ispoly) b.addGeom(b.geoms[0]);
	SpatVector part;
	for (size_t j =0; j<(b.size()-1); j++) {
		std::vector<unsigned> range = {(unsigned)j, (unsigned)j+1};
		SpatVector g = b.subset_rows(range);
		g = g.hull("convex");
		part.addGeom(g.geoms[0]);
	}	
	part = part.aggregate(true);
	return part.geoms[0];
}


SpatVector lonlat_buf(SpatVector x, double dist, unsigned quadsegs, bool ispol, bool ishole) {


	if ((x.extent.ymin > -60) && (x.extent.ymax < 60) && ((x.extent.ymax - x.extent.ymin) < 1) && dist < 110000) {
		x.setSRS("+proj=merc");
		double f = 0.5 - (dist / 220000);
		double halfy = x.extent.ymin + f * (x.extent.ymax - x.extent.ymin);
		std::vector<double> dd = destpoint_lonlat(0, halfy, 0, dist);
		dist = dd[1] - halfy;
		return x.buffer({dist}, quadsegs);
	} 

	SpatVector tmp;
	x = x.disaggregate();
	for (size_t i =0; i<x.geoms.size(); i++) {
		SpatVector p(x.geoms[i]);
		p.srs = x.srs;
		p = p.as_points(false, true);
		std::vector<double> d(p.size(), dist);
		SpatVector b = p.point_buffer(d, quadsegs, true);
		if (b.size() <= p.size()) {	
			SpatGeom g = hullify(b, ispol);
			tmp.addGeom(g);
		} else {			
			SpatVector west, east, eastwest;
			for (size_t j =0; j<b.size(); j++) {
				if ((b.geoms[j].extent.xmin < -179.99) && (b.geoms[j].extent.xmax > 179.99)) {
					tmp.addGeom(b.geoms[j]);							
				} else if (b.geoms[j].extent.xmax < 0) {
					west.addGeom(b.geoms[j]);					
				} else {
					east.addGeom(b.geoms[j]);										
				}	
			}
			if (east.nrow() > 0) {
				SpatGeom geast = hullify(east, ispol);
				tmp.addGeom(geast);
			}
			if (west.nrow() > 0) {
				SpatGeom gwest = hullify(west, ispol);
				tmp.addGeom(gwest);
			}
		}
	}
	tmp = tmp.aggregate(true);

	if (ispol) {
		tmp = ishole ? tmp.get_holes() : tmp.remove_holes();
	}

	return tmp;
}



SpatVector SpatVector::buffer(std::vector<double> dist, unsigned quadsegs) {

	quadsegs = std::min(quadsegs, (unsigned) 180);
	SpatVector out;
	out.srs = srs;
	if (srs.is_empty()) {
		out.setError("crs not defined");
		return(out);
	}
	bool islonlat = is_lonlat();
	std::string vt = type();
	if (vt == "points" || vt == "lines" || islonlat) {
		for (size_t i=0; i<dist.size(); i++) {
			if (dist[i] <= 0) {
				dist[i] = -dist[i];
			}
		}	
	}
	recycle(dist, size());

	if (islonlat) {
		if (vt == "points") {
			return point_buffer(dist, quadsegs, false);
		} else {
			SpatVector p;
			bool ispol = vt == "polygons";
			for (size_t i =0; i<size(); i++) {
				p = subset_rows(i);
				if (ispol) {
					SpatVector h = p.get_holes();
					p = p.remove_holes();
					p = lonlat_buf(p, dist[i], quadsegs, true, false);
					if (h.size() > 0) {
						h = lonlat_buf(h, dist[i], quadsegs, true, true);
						if (h.size() > 0) {
							for (size_t j=0; j<h.geoms[0].parts.size(); j++) {
								p.geoms[0].parts[0].addHole(h.geoms[0].parts[j].x, h.geoms[0].parts[j].y);
							}
						}
					}
				} else {
					p = lonlat_buf(p, dist[i], quadsegs, false, false);
				}
				out = out.append(p, true);
			}	
			out.df = df;
			return out;	
		}
	}

	
	
	GEOSContextHandle_t hGEOSCtxt = geos_init();
//	SpatVector f = remove_holes();

	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> b(size());
	for (size_t i = 0; i < g.size(); i++) {
		GEOSGeometry* pt = GEOSBuffer_r(hGEOSCtxt, g[i].get(), dist[i], quadsegs);
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
	out = coll.get(0);
	out.srs = srs;
	out.df = df;
	
	return out;
}


SpatVector SpatVector::intersect(SpatVector v) {

	SpatVector out;
	out.srs = srs;

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	//v = v.aggregate(false);
	std::vector<GeomPtr> y = geos_geoms(&v, hGEOSCtxt);
	std::vector<GeomPtr> result;
	size_t nx = size();
	size_t ny = v.size();
	std::vector<unsigned> idx, idy;
	idx.reserve(nx);
	idy.reserve(ny);

	if (type() == "points") {
		//std::vector<bool> ixj(nx, false);
		//size_t count = 0;
		for (size_t j = 0; j < ny; j++) {
			PrepGeomPtr pr = geos_ptr(GEOSPrepare_r(hGEOSCtxt, y[j].get()), hGEOSCtxt);
			for (size_t i = 0; i < nx; i++) {
				if (GEOSPreparedIntersects_r(hGEOSCtxt, pr.get(), x[i].get())) {
					//if (!ixj[i]
					//ixj[i] = true;
					idx.push_back(i);
					idy.push_back(j);
					//count++;
				}
			}
		}
		//std::vector<int> sx;
		//sx.reserve(count);
		//for (size_t i=0; i<ixj.size(); i++) {
		//	if (ixj[i]) sx.push_back(i);	
		//}
		out = subset_rows(idx);
		
	} else {
		
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
					idx.push_back(i);
					idy.push_back(j);
				} else {
					GEOSGeom_destroy_r(hGEOSCtxt, geom);
				}
			}
		}
	//SpatVectorCollection coll = coll_from_geos(result, hGEOSCtxt);
		if (result.size() > 0) {
			SpatVectorCollection coll = coll_from_geos(result, hGEOSCtxt);
			out = coll.get(0);
			out.srs = srs;
		}
	}
	geos_finish(hGEOSCtxt);

	if (!srs.is_same(v.srs, true)) {
		out.addWarning("different crs"); 
	}

	if ((type() == "polygons") && (v.type() == "polygons") && (out.type() != "polygons")) {
		// intersection is point or line, return empty 
		out = SpatVector();
		out.addWarning("no intersection");
		out.srs = srs;
	}

	SpatDataFrame df1 = df.subset_rows(idx);
	SpatDataFrame df2 = v.df.subset_rows(idy);
	if (!df1.cbind(df2)) {
		out.addWarning("could not combine attributes");
	}
	out.df = df1;
	return out;
}


/*

std::function<char(GEOSContextHandle_t, const GEOSGeometry *, const GEOSGeometry *)> getRelateFun(const std::string rel) {
	std::function<char(GEOSContextHandle_t, const GEOSGeometry *, const GEOSGeometry *)> rfun;
	if (rel == "intersects") {
		rfun = GEOSIntersects_r;
	} else if (rel == "disjoint") {
		rfun = GEOSDisjoint_r;
	} else if (rel == "touches") {
		rfun = GEOSTouches_r;
	} else if (rel == "crosses") {
		rfun = GEOSCrosses_r;
	} else if (rel == "within") {
		rfun = GEOSWithin_r;
	} else if (rel == "contains") {
		rfun = GEOSContains_r;
	} else if (rel == "containsproperly") {
		rfun = GEOSContainsProperly_r;
	} else if (rel == "overlaps") {
		rfun = GEOSOverlaps_r;
	} else if (rel == "covers") {
		rfun = GEOSCovers_r;
	} else if (rel == "coveredby") {
		rfun = GEOSCoveredBy_r;
	}
	return rfun;
}
*/

std::function<char(GEOSContextHandle_t, const GEOSPreparedGeometry *, const GEOSGeometry *)> getPrepRelateFun(const std::string rel) {
	std::function<char(GEOSContextHandle_t, const GEOSPreparedGeometry *, const GEOSGeometry *)> rfun;
	if (rel == "intersects") {
		rfun = GEOSPreparedIntersects_r;
	} else if (rel == "disjoint") {
		rfun = GEOSPreparedDisjoint_r;
	} else if (rel == "touches") {
		rfun = GEOSPreparedTouches_r;
	} else if (rel == "crosses") {
		rfun = GEOSPreparedCrosses_r;
	} else if (rel == "within") {
		rfun = GEOSPreparedWithin_r;
	} else if (rel == "contains") {
		rfun = GEOSPreparedContains_r;
	} else if (rel == "containsproperly") {
		rfun = GEOSPreparedContainsProperly_r;
	} else if (rel == "overlaps") {
		rfun = GEOSPreparedOverlaps_r;
	} else if (rel == "covers") {
		rfun = GEOSPreparedCovers_r;
	} else if (rel == "coveredby") {
		rfun = GEOSPreparedCoveredBy_r;
	}
	return rfun;
}


int getRel(std::string &relation) {
	int pattern = 1;
	std::string rel = relation;
	std::transform(rel.begin(), rel.end(), rel.begin(), ::tolower);
	std::vector<std::string> f {"rook", "queen", "intersects", "touches", "crosses", "overlaps", "within", "contains", "covers", "coveredby", "disjoint"};
	if (std::find(f.begin(), f.end(), rel) == f.end()) {
		if (relation.size() != 9) {
			pattern = 2;
		} else {
			std::string r = relation;
			for (size_t i=0; i<9; i++) {
				if (!(r.at(i) == 'T' || r.at(i) == 'F' || r.at(i) == '0' || r.at(i) == '1' || r.at(i) == '2' || r.at(i) == '*')) {
					pattern = 2;
					break;
				}
			}
		}
	} else if (rel == "rook") {
		relation = "F***1****";
	} else if (rel == "queen") {
		relation = "F***T****";
	} else {
		pattern = 0;
		relation = rel;
	}
	return pattern;
}

std::vector<int> SpatVector::relate(SpatVector v, std::string relation) {

	std::vector<int> out;
	int pattern = getRel(relation);
	if (pattern == 2) {
		setError("'" + relation + "'" + " is not a valid relate name or pattern");
		return out;
	}

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> y = geos_geoms(&v, hGEOSCtxt);
	size_t nx = size();
	size_t ny = v.size();
	out.reserve(nx*ny);
	if (pattern == 1) {
		for (size_t i = 0; i < nx; i++) {
			for (size_t j = 0; j < ny; j++) {
				out.push_back( GEOSRelatePattern_r(hGEOSCtxt, x[i].get(), y[j].get(), relation.c_str()));
			}
		}
	} else {
		std::function<char(GEOSContextHandle_t, const GEOSPreparedGeometry *, const GEOSGeometry *)> relFun = getPrepRelateFun(relation);
		for (size_t i = 0; i < nx; i++) {
			PrepGeomPtr pr = geos_ptr(GEOSPrepare_r(hGEOSCtxt, x[i].get()), hGEOSCtxt);
			for (size_t j = 0; j < ny; j++) {
				out.push_back( relFun(hGEOSCtxt, pr.get(), y[j].get()));
			}
		} 
	}
	geos_finish(hGEOSCtxt);

	return out;
}


std::vector<int> SpatVector::relateFirst(SpatVector v, std::string relation) {

	int pattern = getRel(relation);
	if (pattern == 2) {
		setError("'" + relation + "'" + " is not a valid relate name or pattern");
		std::vector<int> out;
		return out;
	}
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> y = geos_geoms(&v, hGEOSCtxt);
	size_t nx = size();
	size_t ny = v.size();
	std::vector<int> out(nx, -1);
	if (pattern == 1) {
		for (size_t i = 0; i < nx; i++) {
			for (size_t j = 0; j < ny; j++) {
				if (GEOSRelatePattern_r(hGEOSCtxt, x[i].get(), y[j].get(), relation.c_str())) {
					out[i] = j;
					continue;
				}
			}
		}
	} else {
		//std::function<char(GEOSContextHandle_t, const GEOSGeometry *, const GEOSGeometry *)> relFun = getRelateFun(relation);
		std::function<char(GEOSContextHandle_t, const GEOSPreparedGeometry *, const GEOSGeometry *)> relFun = getPrepRelateFun(relation);

		for (size_t i = 0; i < nx; i++) {
			PrepGeomPtr pr = geos_ptr(GEOSPrepare_r(hGEOSCtxt, x[i].get()), hGEOSCtxt);
			for (size_t j = 0; j < ny; j++) {
				if (relFun(hGEOSCtxt, pr.get(), y[j].get())) {
					out[i] = j;
					continue;
				}
			}
		}
	}
	geos_finish(hGEOSCtxt);
	return out;
}



std::vector<int> SpatVector::relate(std::string relation, bool symmetrical) {

	std::vector<int> out;
	int pattern = getRel(relation);
	if (pattern == 2) {
		setError("'" + relation + "'" + " is not a valid relate name or pattern");
		return out;
	}

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);

	if (symmetrical) {
		size_t s = size();
		size_t n = ((s-1) * s)/2;
		out.reserve(n);
		if (pattern == 1) {
			for (size_t i=0; i<(s-1); i++) {
				for (size_t j=(i+1); j<s; j++) {
					out.push_back( GEOSRelatePattern_r(hGEOSCtxt, x[i].get(), x[j].get(), relation.c_str()));
				}
			}
		} else {
			std::function<char(GEOSContextHandle_t, const GEOSPreparedGeometry *, const GEOSGeometry *)> relFun = getPrepRelateFun(relation);
			for (size_t i=0; i<(s-1); i++) {
				PrepGeomPtr pr = geos_ptr(GEOSPrepare_r(hGEOSCtxt, x[i].get()), hGEOSCtxt);
				for (size_t j=(i+1); j<s; j++) {
					out.push_back( relFun(hGEOSCtxt, pr.get(), x[j].get()));
				}
			} 			
		}
	} else {
		size_t nx = size();
		out.reserve(nx*nx);
		if (pattern == 1) {
			for (size_t i = 0; i < nx; i++) {
				for (size_t j = 0; j < nx; j++) {
					out.push_back( GEOSRelatePattern_r(hGEOSCtxt, x[i].get(), x[j].get(), relation.c_str()));
				}
			}
		} else {
			std::function<char(GEOSContextHandle_t, const GEOSPreparedGeometry *, const GEOSGeometry *)> relFun = getPrepRelateFun(relation);
			for (size_t i = 0; i < nx; i++) {
				PrepGeomPtr pr = geos_ptr(GEOSPrepare_r(hGEOSCtxt, x[i].get()), hGEOSCtxt);
				for (size_t j = 0; j < nx; j++) {
					out.push_back( relFun(hGEOSCtxt, pr.get(), x[j].get()));
				}
			} 
		}
	}

	geos_finish(hGEOSCtxt);

	return out;
}



std::vector<double> SpatVector::geos_distance(SpatVector v, bool parallel) {

	std::vector<double> out;

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> y = geos_geoms(&v, hGEOSCtxt);
	size_t nx = size();
	size_t ny = v.size();
	double d;

	if (parallel) {
		bool nyone = false;
		if (nx != ny) {
			if (ny == 1) {
				nyone = true;
			} else if ((nx == 1) && (ny > 1)) {
				std::swap(x, y);
				std::swap(nx, ny);
				nyone = true;
			} else {
				// recycle?
				setError("vectors have different lengths");
				return out;
			}
		}
		if (nyone) {
			out.reserve(nx);
			for (size_t i = 0; i < nx; i++) {
				if ( GEOSDistance_r(hGEOSCtxt, x[i].get(), y[0].get(), &d)) {
					out.push_back(d);
				} else {
					out.push_back(NAN);				
				}
			}
		} else {
			out.reserve(nx);
			for (size_t i = 0; i < nx; i++) {
				if ( GEOSDistance_r(hGEOSCtxt, x[i].get(), y[i].get(), &d)) {
					out.push_back(d);
				} else {
					out.push_back(NAN);				
				}
			}
		}
	} else {
		out.reserve(nx*ny);
		for (size_t i = 0; i < nx; i++) {
			for (size_t j = 0; j < ny; j++) {
				if ( GEOSDistance_r(hGEOSCtxt, x[i].get(), y[j].get(), &d)) {
					out.push_back(d);
				} else {
					out.push_back(NAN);				
				}
			}
		}
	}
	geos_finish(hGEOSCtxt);
	return out;
 }

std::vector<double> SpatVector::geos_distance(bool sequential) {

	std::vector<double> out;

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	size_t s = size();
	double d;
	if (sequential) {
		out.reserve(s);
		out.push_back(0);
		for (size_t i=0; i<(s-1); i++) {
			if ( GEOSDistance_r(hGEOSCtxt, x[i].get(), x[i+1].get(), &d)) {
				out.push_back(d);
			} else {
				out.push_back(NAN);				
			}
		}
	} else {
		out.reserve((s-1) * s / 2);
		for (size_t i=0; i<(s-1); i++) {
			for (size_t j=(i+1); j<s; j++) {
				if ( GEOSDistance_r(hGEOSCtxt, x[i].get(), x[j].get(), &d)) {
					out.push_back(d);
				} else {
					out.push_back(NAN);				
				}
			}
		}
	}
	if (s == 1) {
		out.push_back(0);
	}	
	geos_finish(hGEOSCtxt);
	return out;
 }


SpatVector SpatVector::unite(SpatVector v) {
	SpatVector intsec = intersect(v);
	if (intsec.hasError()) {
		return intsec;
	}
	SpatVector sdif = symdif(v);
	if (sdif.hasError()) {
		return intsec;
	}
	return intsec.append(sdif, true);
}


SpatVector SpatVector::unite() {
	int n = size();

	std::vector<long> x(1, 1);
	SpatDataFrame d;
	d.add_column(x, "id_1");
	SpatVector out = subset_rows(0);
	out.df = d;
	for (int i=1; i<n; i++) {
		std::string name = "id_" + std::to_string(i+1);
		SpatDataFrame d;
		d.add_column(x, name);
		SpatVector r = subset_rows(i);
		r.df = d;
		out = out.unite(r);
	}

	for (size_t i=0; i<out.df.iv.size(); i++) {
		for (size_t j=0; j<out.df.iv[i].size(); j++) {
			if (out.df.iv[i][j] != 1) {
				out.df.iv[i][j] = 0;
			}
		}
	}

	return out;
}



SpatVector SpatVector::symdif(SpatVector v) {
	if ((type() != "polygons") || (v.type() != "polygons")) {
		SpatVector out;
		out.setError("expect two polygon geometries");
		return out;
	}
	SpatVector out = erase(v);
	out = out.append(v.erase(*this), true);
	return out;
	
/*
	SpatVector out;
	out.srs = srs;

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> y = geos_geoms(&v, hGEOSCtxt);
	std::vector<GeomPtr> result;
	std::vector<unsigned> ids;
	ids.reserve(size());
	size_t nx = size();
	size_t ny = v.size();

	
	for (size_t i = 0; i < nx; i++) {
		GEOSGeometry* geom = x[i].get();
		for (size_t j = 0; j < ny; j++) {
			geom = GEOSDifference_r(hGEOSCtxt, geom, y[j].get());
			if (geom == NULL) {
				out.setError("GEOS exception");
				geos_finish(hGEOSCtxt);
				return(out);
			} 
			if (GEOSisEmpty_r(hGEOSCtxt, geom)) {
				break;
			}
		}
		if (!GEOSisEmpty_r(hGEOSCtxt, geom)) {
			result.push_back(geos_ptr(geom, hGEOSCtxt));
			ids.push_back(i);
		} else {
			GEOSGeom_destroy_r(hGEOSCtxt, geom);
		}
	}

	if (result.size() > 0) {
		SpatVectorCollection coll = coll_from_geos(result, hGEOSCtxt);
		out = coll.get(0);
		out.srs = srs;
		out.df = df.subset_rows(ids);
	} 
	geos_finish(hGEOSCtxt);
	if (!srs.is_same(v.srs, true)) {
		out.addWarning("different crs"); 
	}

	return out.append(v, true);
	*/
}




SpatVector SpatVector::cover(SpatVector v, bool identity) {
	if (v.srs.is_empty()) {
		v.srs = srs; 
	}
	SpatVector out = erase(v);
	if (identity) {
		SpatVector insect = intersect(v);
		v = v.erase(insect);
		out = out.append(insect, true);
		out = out.append(v, true);
	} else {
		out = out.append(v, true);
	}
	return out;
}


SpatVector SpatVector::erase(SpatVector v) {

	SpatVector out;
	out.srs = srs;

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> y = geos_geoms(&v, hGEOSCtxt);
	std::vector<GeomPtr> result;
	std::vector<unsigned> ids;
	ids.reserve(size());
	size_t nx = size();
	size_t ny = v.size();

	
	for (size_t i = 0; i < nx; i++) {
		GEOSGeometry* geom = x[i].get();
		for (size_t j = 0; j < ny; j++) {
			geom = GEOSDifference_r(hGEOSCtxt, geom, y[j].get());
			if (geom == NULL) {
				out.setError("GEOS exception");
				geos_finish(hGEOSCtxt);
				return(out);
			} 
			if (GEOSisEmpty_r(hGEOSCtxt, geom)) {
				break;
			}
		}
		if (!GEOSisEmpty_r(hGEOSCtxt, geom)) {
			result.push_back(geos_ptr(geom, hGEOSCtxt));
			ids.push_back(i);
		} else {
			GEOSGeom_destroy_r(hGEOSCtxt, geom);	
		}
	}

	if (result.size() > 0) {
		SpatVectorCollection coll = coll_from_geos(result, hGEOSCtxt, ids);
		out = coll.get(0);
		out.df = df.subset_rows(out.df.iv[0]);
		out.srs = srs;
	} 
	geos_finish(hGEOSCtxt);
	if (!srs.is_same(v.srs, true)) {
		out.addWarning("different crs"); 
	}

	return out;
}



SpatVector SpatVector::nearest_point(SpatVector v, bool parallel) {
	SpatVector out;
		
	if ((size() == 0) || (v.size()==0)) {
		out.setError("empty SpatVecor(s)");
		return out;
	}
	if (!srs.is_equal(v.srs)) {
		out.setError("CRSs do not match");
		return out;
	}
	out.srs = srs;
	
	if (is_lonlat()) {
		if (type() == "points") {
			std::vector<double> nlon, nlat, dist;
			std::vector<long> id;
			std::vector<std::vector<double>> p = coordinates();
			std::vector<std::vector<double>> pv = v.coordinates();
			nearest_lonlat(id, dist, nlon, nlat, p[0], p[1], pv[0], pv[1]);
			out.setPointsGeometry(nlon, nlat);
			std::vector<long> fromid(id.size());
			std::iota(fromid.begin(), fromid.end(), 0);
			out.df.add_column(fromid, "from_id");
			out.df.add_column(p[0], "from_x");
			out.df.add_column(p[1], "from_y");	
			out.df.add_column(id, "to_id");
			out.df.add_column(nlon, "to_x");
			out.df.add_column(nlat, "to_y");
			out.df.add_column(dist, "distance");					
			return out;
		} else {
			out.setError("not yet implement for non-point lonlat vector data");
			return out;
		}		
	}

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	if (parallel) {
		if ((size() != v.size())) {
			out.setError("SpatVecors do not have the same size");
			return out;
		}		
		std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
		std::vector<GeomPtr> y = geos_geoms(&v, hGEOSCtxt);
		std::vector<GeomPtr> b(size());
		for (size_t i=0; i < x.size(); i++) {	
			GEOSCoordSequence* csq = GEOSNearestPoints_r(hGEOSCtxt, x[i].get(), y[i].get());
			GEOSGeometry* geom = GEOSGeom_createLineString_r(hGEOSCtxt, csq);
			b[i] = geos_ptr(geom, hGEOSCtxt);
		}
		out = vect_from_geos(b, hGEOSCtxt, "lines");
		
	} else {	
		SpatVector mp = v.aggregate(false);
		std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
		std::vector<GeomPtr> y = geos_geoms(&mp, hGEOSCtxt);
		std::vector<GeomPtr> b(size());
		for (size_t i = 0; i < x.size(); i++) {	
			GEOSCoordSequence* csq = GEOSNearestPoints_r(hGEOSCtxt, x[i].get(), y[0].get());
			GEOSGeometry* geom = GEOSGeom_createLineString_r(hGEOSCtxt, csq);
			b[i] = geos_ptr(geom, hGEOSCtxt);
		}
		out = vect_from_geos(b, hGEOSCtxt, "lines");
	}
	geos_finish(hGEOSCtxt);
	out.srs = srs;
	return out;
}

SpatVector SpatVector::nearest_point() {
	SpatVector out;
	if ((size() == 0)) {
		out.addWarning("empty SpatVecor");
		return out;
	}
	if ((size() == 1)) {
		out.addWarning("single geometry");
		//return *this;
	}
	size_t n = size();
	out.srs = srs;
	
	if (is_lonlat()) {
		if (type() == "points") {
			std::vector<double> nlon, nlat, dist;
			std::vector<long> id;
			std::vector<std::vector<double>> p = coordinates();
			nearest_lonlat_self(id, dist, nlon, nlat, p[0], p[1]);
			out.setPointsGeometry(nlon, nlat);
			out.df.add_column(id, "id");
			out.df.add_column(dist, "distance");		
			return out;
		} else {
			out.setError("not yet implement for non-point lonlat vector data");
			return out;
		}		
	}


	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> b(n);
	for (unsigned i = 0; i < n; i++) {	
		SpatVector xa = remove_rows({i});
		xa = xa.aggregate(false);
		std::vector<GeomPtr> y = geos_geoms(&xa, hGEOSCtxt);
		GEOSCoordSequence* csq = GEOSNearestPoints_r(hGEOSCtxt, x[i].get(), y[0].get());
		GEOSGeometry* geom = GEOSGeom_createLineString_r(hGEOSCtxt, csq);
		b[i] = geos_ptr(geom, hGEOSCtxt);
	}
	out = vect_from_geos(b, hGEOSCtxt, "lines");
	geos_finish(hGEOSCtxt);
	out.srs = srs;
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

	out.srs = srs;
	out.df = df;
	return out;
}

SpatVector SpatVector::unaryunion() {
	SpatVector out;

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> gout(g.size());
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



std::vector<double> SpatVector::width() {
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> gout(g.size());
	SpatVector tmp;
	for (size_t i = 0; i < g.size(); i++) {	
		GEOSGeometry* w = GEOSMinimumWidth_r(hGEOSCtxt, g[i].get());
		if (w == NULL) {
			setError("NULL geom");
			geos_finish(hGEOSCtxt);
			std::vector<double> out;
			return out;
		}
		gout[i] = geos_ptr(w, hGEOSCtxt);
	}
	SpatVectorCollection coll = coll_from_geos(gout, hGEOSCtxt);
	geos_finish(hGEOSCtxt);
	tmp = coll.get(0);
	tmp.srs = srs;
	return tmp.length();
}


std::vector<double> SpatVector::clearance() {
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> g = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> gout(g.size());
	SpatVector tmp;
	for (size_t i = 0; i < g.size(); i++) {	
		GEOSGeometry* w = GEOSMinimumClearanceLine_r(hGEOSCtxt, g[i].get());
		if (w == NULL) {
			setError("NULL geom");
			geos_finish(hGEOSCtxt);
			std::vector<double> out;
			return out;
		}
		gout[i] = geos_ptr(w, hGEOSCtxt);
	}
	SpatVectorCollection coll = coll_from_geos(gout, hGEOSCtxt);
	geos_finish(hGEOSCtxt);
	tmp = coll.get(0);
	tmp.srs = srs;
	return tmp.length();
}

