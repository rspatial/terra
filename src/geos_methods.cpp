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
		}
	}
	if (p.size() > 0) {
		SpatVectorCollection coll = coll_from_geos(p, hGEOSCtxt);
		out = coll.get(0);
		out.srs = srs;
		out.df = df.subset_rows(id);
	}
	geos_finish(hGEOSCtxt);
	return out;
}



SpatVector SpatVector::crop(SpatVector v) {

	SpatVector out;
	out.srs = srs;

	GEOSContextHandle_t hGEOSCtxt = geos_init();

	std::vector<GeomPtr> x;
	if ((type() != "polygons") & (type() != "mutlipolygons")) {
		out = convexhull();
		x = geos_geoms(&out, hGEOSCtxt);
	} else {
		x = geos_geoms(this, hGEOSCtxt);
	}
	std::vector<GeomPtr> y = geos_geoms(&v, hGEOSCtxt);
	std::vector<GeomPtr> result;
	std::vector<unsigned> ids;
	ids.reserve(size());
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
				ids.push_back(i);
			}
		}
	}

	SpatVectorCollection coll = coll_from_geos(result, hGEOSCtxt);
	geos_finish(hGEOSCtxt);

	if (result.size() > 0) {
		SpatVectorCollection coll = coll_from_geos(result, hGEOSCtxt);
		out = coll.get(0);
		out.srs = srs;
		out.df = df.subset_rows(ids);
	} 
	return out;
}



SpatVector SpatVector::convexhull() {
	SpatVector out;
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	SpatVector a = aggregate(false);
	std::vector<GeomPtr> g = geos_geoms(&a, hGEOSCtxt);
	std::string vt = type();
	GEOSGeometry* h = GEOSConvexHull_r(hGEOSCtxt, g[0].get());
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
	out.srs = srs;

#ifndef HAVE350
	out.setError("GEOS 3.5 required for voronoi")
	return out;
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
	out.srs = srs;
	if (!out.hasError()) {
		out = out.disaggregate();
		// associate with attributes
	}
#endif
	return out;
}  

SpatVector SpatVector::buffer2(double dist, unsigned nQuadSegs, unsigned capstyle) {

	SpatVector out;
	out.srs = srs;

	std::string vt = type();
	if ((vt == "points") && (dist <= 0)) {
		dist = -dist;
	}

	GEOSContextHandle_t hGEOSCtxt = geos_init();
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
	out = coll.get(0);
	out.srs = srs;

	return out;
}


SpatVector SpatVector::intersect(SpatVector v) {

	SpatVector out;
	out.srs = srs;

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> y = geos_geoms(&v, hGEOSCtxt);
	std::vector<GeomPtr> result;
	std::vector<std::vector<unsigned>> ids(2);
	ids[0].reserve(size());
	ids[1].reserve(size());
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
				ids[0].push_back(i);
				ids[1].push_back(j);
				//Rcpp::Rcout << "i: " << i << " j: " << j << std::endl;
			}
		}
	}

	SpatVectorCollection coll = coll_from_geos(result, hGEOSCtxt);
	geos_finish(hGEOSCtxt);

	if (result.size() > 0) {
		SpatVectorCollection coll = coll_from_geos(result, hGEOSCtxt);
		out = coll.get(0);
		out.srs = srs;
		SpatDataFrame df1 = df.subset_rows(ids[0]);
		SpatDataFrame df2 = v.df.subset_rows(ids[1]);
		if (!df1.cbind(df2)) {
			out.addWarning("could not bind attributes");
		}
		out.df = df1;
	} 
	if (!srs.is_same(v.srs, true)) {
		out.addWarning("different crs"); 
	}

	return out;
}


std::vector<int> SpatVector::relate(SpatVector v, const std::string relation) {

	std::vector<int> out;
	std::vector<std::string> f {"intersects", "touches", "crosses", "overlaps", "within", "contains", "covers", "coveredby", "disjoint"};
	if (std::find(f.begin(), f.end(), relation) == f.end()) {
		setError("invalid relation name");
		return out;
	}

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> y = geos_geoms(&v, hGEOSCtxt);
	size_t nx = size();
	size_t ny = v.size();
	out.reserve(nx*ny);
	for (size_t i = 0; i < nx; i++) {
		for (size_t j = 0; j < ny; j++) {
			if (relation == "intersects") {
				out.push_back( GEOSIntersects_r(hGEOSCtxt, x[i].get(), y[j].get()));
			} else if (relation == "disjoint") {
				out.push_back( GEOSDisjoint_r(hGEOSCtxt, x[i].get(), y[j].get()));
			} else if (relation == "touches") {
				out.push_back( GEOSTouches_r(hGEOSCtxt, x[i].get(), y[j].get()));
			} else if (relation == "crosses") {
				out.push_back( GEOSCrosses_r(hGEOSCtxt, x[i].get(), y[j].get()));
			} else if (relation == "within") {
				out.push_back( GEOSWithin_r(hGEOSCtxt, x[i].get(), y[j].get()));
			} else if (relation == "contains") {
				out.push_back( GEOSContains_r(hGEOSCtxt, x[i].get(), y[j].get()));
			} else if (relation == "overlaps") {
				out.push_back( GEOSOverlaps_r(hGEOSCtxt, x[i].get(), y[j].get()));
			} else if (relation == "covers") {
				out.push_back( GEOSCovers_r(hGEOSCtxt, x[i].get(), y[j].get()));
			} else if (relation == "coveredby") {
				out.push_back( GEOSCoveredBy_r(hGEOSCtxt, x[i].get(), y[j].get()));
			}
		}
	} 
	
	geos_finish(hGEOSCtxt);
	return out;
}



std::vector<int> SpatVector::relate(const std::string relation) {

	std::vector<int> out;
	std::vector<std::string> f {"intersects", "touches", "crosses", "overlaps", "within", "contains", "covers", "coveredby", "disjoint"};
	if (std::find(f.begin(), f.end(), relation) == f.end()) {
		setError("invalid relation name");
		return out;
	}

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	size_t nx = size();
	out.reserve(nx*nx);
	for (size_t i = 0; i < nx; i++) {
		for (size_t j = 0; j < nx; j++) {
			if (relation == "intersects") {
				out.push_back( GEOSIntersects_r(hGEOSCtxt, x[i].get(), x[j].get()));
			} else if (relation == "disjoint") {
				out.push_back( GEOSDisjoint_r(hGEOSCtxt, x[i].get(), x[j].get()));
			} else if (relation == "touches") {
				out.push_back( GEOSTouches_r(hGEOSCtxt, x[i].get(), x[j].get()));
			} else if (relation == "crosses") {
				out.push_back( GEOSCrosses_r(hGEOSCtxt, x[i].get(), x[j].get()));
			} else if (relation == "within") {
				out.push_back( GEOSWithin_r(hGEOSCtxt, x[i].get(), x[j].get()));
			} else if (relation == "contains") {
				out.push_back( GEOSContains_r(hGEOSCtxt, x[i].get(), x[j].get()));
			} else if (relation == "overlaps") {
				out.push_back( GEOSOverlaps_r(hGEOSCtxt, x[i].get(), x[j].get()));
			} else if (relation == "covers") {
				out.push_back( GEOSCovers_r(hGEOSCtxt, x[i].get(), x[j].get()));
			} else if (relation == "coveredby") {
				out.push_back( GEOSCoveredBy_r(hGEOSCtxt, x[i].get(), x[j].get()));
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
				setError("vectors have different lenghts");
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

std::vector<double> SpatVector::geos_distance() {

	std::vector<double> out;

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	size_t s = size();
	out.reserve((s-1) * s / 2);
	double d;
	for (size_t i=0; i<(s-1); i++) {
		for (size_t j=(i+1); j<s; j++) {
			if ( GEOSDistance_r(hGEOSCtxt, x[i].get(), x[j].get(), &d)) {
				out.push_back(d);
			} else {
				out.push_back(NAN);				
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
	SpatVector out = intersect(v);
	out = out.append(erase(v), true);
	out = out.append(v.erase(*this), true);
	return out;
}


SpatVector SpatVector::symdif(SpatVector v) {
	SpatVector out = erase(v);
	v = v.erase(*this);
	return out.append(v, true);
}



SpatVector SpatVector::cover(SpatVector v, bool identity) {
	if (v.srs.is_empty()) {
		v.srs = srs; 
	}
	SpatVector out = erase(v);
	out = identity ? out.append(intersect(v), true) : out.append(v, true);
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

	return out;
}



SpatVector SpatVector::nearest_point(SpatVector v, bool parallel) {
	SpatVector out;
	if ((size() == 0) || (v.size()==0)) {
		out.setError("empty SpatVecor(s)");
		return out;
	}
	if (parallel) {
		if ((size() != v.size())) {
			out.setError("SpatVecors do not have the same size");
			return out;
		}		
		GEOSContextHandle_t hGEOSCtxt = geos_init();
		std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
		std::vector<GeomPtr> y = geos_geoms(&v, hGEOSCtxt);
		std::vector<GeomPtr> b(size());
		for (size_t i=0; i < x.size(); i++) {	
			GEOSCoordSequence* csq = GEOSNearestPoints_r(hGEOSCtxt, x[i].get(), y[i].get());
			GEOSGeometry* geom = GEOSGeom_createPoint_r(hGEOSCtxt, csq);
			b[i] = geos_ptr(geom, hGEOSCtxt);
		}
		out = vect_from_geos(b, hGEOSCtxt, "points");
		geos_finish(hGEOSCtxt);
		
	} else {	
		SpatVector mp = v.aggregate(false);
		GEOSContextHandle_t hGEOSCtxt = geos_init();
		std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
		std::vector<GeomPtr> y = geos_geoms(&mp, hGEOSCtxt);
		std::vector<GeomPtr> b(size());
		for (size_t i = 0; i < x.size(); i++) {	
			GEOSCoordSequence* csq = GEOSNearestPoints_r(hGEOSCtxt, x[i].get(), y[0].get());
			GEOSGeometry* geom = GEOSGeom_createPoint_r(hGEOSCtxt, csq);
			b[i] = geos_ptr(geom, hGEOSCtxt);
		}
		out = vect_from_geos(b, hGEOSCtxt, "points");
		geos_finish(hGEOSCtxt);
	}
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
	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	SpatVector xa = aggregate(false);
	std::vector<GeomPtr> y = geos_geoms(&xa, hGEOSCtxt);
	std::vector<GeomPtr> b(size());
	for (size_t i = 0; i < x.size(); i++) {	
		GEOSCoordSequence* csq = GEOSNearestPoints_r(hGEOSCtxt, x[i].get(), y[0].get());
		GEOSGeometry* geom = GEOSGeom_createPoint_r(hGEOSCtxt, csq);
		b[i] = geos_ptr(geom, hGEOSCtxt);
	}
	out = vect_from_geos(b, hGEOSCtxt, "points");
	geos_finish(hGEOSCtxt);
	return out;
}



SpatVector SpatVector::centroid() {

	SpatVector out;
	out.srs = srs;

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
	SpatVector out;
	out.srs = srs;

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

