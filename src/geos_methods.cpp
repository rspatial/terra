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
			} else {
				GEOSGeom_destroy_r(hGEOSCtxt, geom);
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
	out.df = df;
	
	return out;
}


SpatVector SpatVector::intersect(SpatVector v) {

	SpatVector out;
	out.srs = srs;

	GEOSContextHandle_t hGEOSCtxt = geos_init();
	std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
	std::vector<GeomPtr> y = geos_geoms(&v, hGEOSCtxt);
	std::vector<GeomPtr> result;
	size_t nx = size();
	size_t ny = v.size();
	std::vector<unsigned> idx, idy;
	idx.reserve(nx);
	idy.reserve(ny);
		
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
		SpatDataFrame df1 = df.subset_rows(idx);
		SpatDataFrame df2 = v.df.subset_rows(idy);
		if (!df1.cbind(df2)) {
			out.addWarning("could not combine attributes");
		}
		out.df = df1;
	} 
	if (!srs.is_same(v.srs, true)) {
		out.addWarning("different crs"); 
	}
	geos_finish(hGEOSCtxt);

	if ((type() == "polygons") && (v.type() == "polygons") && (out.type() != "polygons")) {
		// intersection is point or line, return empty 
		out = SpatVector();
		out.addWarning("no intersection");
		out.srs = srs;
	}

	return out;
}




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
	} else if (rel == "overlaps") {
		rfun = GEOSOverlaps_r;
	} else if (rel == "covers") {
		rfun = GEOSCovers_r;
	} else if (rel == "coveredby") {
		rfun = GEOSCoveredBy_r;
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
		std::function<char(GEOSContextHandle_t, const GEOSGeometry *, const GEOSGeometry *)> relFun = getRelateFun(relation);
		for (size_t i = 0; i < nx; i++) {
			for (size_t j = 0; j < ny; j++) {
				out.push_back( relFun(hGEOSCtxt, x[i].get(), y[j].get()));
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
		std::function<char(GEOSContextHandle_t, const GEOSGeometry *, const GEOSGeometry *)> relFun = getRelateFun(relation);

		for (size_t i = 0; i < nx; i++) {
			for (size_t j = 0; j < ny; j++) {
				if (relFun(hGEOSCtxt, x[i].get(), y[j].get())) {
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
			std::function<char(GEOSContextHandle_t, const GEOSGeometry *, const GEOSGeometry *)> relFun = getRelateFun(relation);
			for (size_t i=0; i<(s-1); i++) {
				for (size_t j=(i+1); j<s; j++) {
					out.push_back( relFun(hGEOSCtxt, x[i].get(), x[j].get()));
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
			std::function<char(GEOSContextHandle_t, const GEOSGeometry *, const GEOSGeometry *)> relFun = getRelateFun(relation);
			for (size_t i = 0; i < nx; i++) {
				for (size_t j = 0; j < nx; j++) {
					out.push_back( relFun(hGEOSCtxt, x[i].get(), x[j].get()));
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
			GEOSGeometry* geom = GEOSGeom_createPoint_r(hGEOSCtxt, csq);
			b[i] = geos_ptr(geom, hGEOSCtxt);
		}
		out = vect_from_geos(b, hGEOSCtxt, "points");
		
	} else {	
		SpatVector mp = v.aggregate(false);
		std::vector<GeomPtr> x = geos_geoms(this, hGEOSCtxt);
		std::vector<GeomPtr> y = geos_geoms(&mp, hGEOSCtxt);
		std::vector<GeomPtr> b(size());
		for (size_t i = 0; i < x.size(); i++) {	
			GEOSCoordSequence* csq = GEOSNearestPoints_r(hGEOSCtxt, x[i].get(), y[0].get());
			GEOSGeometry* geom = GEOSGeom_createPoint_r(hGEOSCtxt, csq);
			b[i] = geos_ptr(geom, hGEOSCtxt);
		}
		out = vect_from_geos(b, hGEOSCtxt, "points");
	}
	geos_finish(hGEOSCtxt);
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

