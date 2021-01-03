#define GEOS_USE_ONLY_R_API
#include <geos_c.h>


#if GEOS_VERSION_MAJOR == 3
# if GEOS_VERSION_MINOR >= 5
#  define HAVE350
# endif
# if GEOS_VERSION_MINOR == 6
#  if GEOS_VERSION_PATCH >= 1
#   define HAVE361
#  endif
# endif
# if GEOS_VERSION_MINOR >= 7
#  define HAVE361
#  define HAVE370
# endif
#else
# if GEOS_VERSION_MAJOR > 3
#  define HAVE350
#  define HAVE370
#  define HAVE361
# endif
#endif

#include "spatVector.h"
#include <cstdarg> 
#include <cstring> 
#include <memory>
#include <functional>


using GeomPtr = std::unique_ptr<GEOSGeometry, std::function<void(GEOSGeometry*)> >;

static GeomPtr geos_ptr(GEOSGeometry* g, GEOSContextHandle_t hGEOSctxt) {
	auto deleter = std::bind(GEOSGeom_destroy_r, hGEOSctxt, std::placeholders::_1);
	return GeomPtr(g, deleter);
}


#include "Rcpp.h"

template <typename... Args>
inline void warnNoCall(const char* fmt, Args&&... args ) {
    Rf_warningcall(R_NilValue, tfm::format(fmt, std::forward<Args>(args)... ).c_str());
}

template <typename... Args>
inline void NORET errNoCall(const char* fmt, Args&&... args) {
    throw Rcpp::exception(tfm::format(fmt, std::forward<Args>(args)... ).c_str(), false);
}


static void __errorHandler(const char *fmt, ...) { 
	char buf[BUFSIZ], *p;
	va_list ap;
	va_start(ap, fmt);
	vsprintf(buf, fmt, ap);
	va_end(ap);
	p = buf + strlen(buf) - 1;
	if(strlen(buf) > 0 && *p == '\n') *p = '\0';
    errNoCall(buf); 
	return; 
} 

static void __warningHandler(const char *fmt, ...) {
	char buf[BUFSIZ], *p;
	va_list ap;
	va_start(ap, fmt);
	vsprintf(buf, fmt, ap);
	va_end(ap);
	p = buf + strlen(buf) - 1;
	if(strlen(buf) > 0 && *p == '\n') *p = '\0';
    warnNoCall(buf); 
	return;
}

void geos_finish(GEOSContextHandle_t ctxt) {
#ifdef HAVE350
	GEOS_finish_r(ctxt);
#else
	finishGEOS_r(ctxt);
#endif
}


GEOSContextHandle_t geos_init(void) {
#ifdef HAVE350
	GEOSContextHandle_t ctxt = GEOS_init_r();
	GEOSContext_setNoticeHandler_r(ctxt, __warningHandler);
	GEOSContext_setErrorHandler_r(ctxt, __errorHandler);
	return ctxt;
#else
	return initGEOS_r((GEOSMessageHandler) __warningHandler, (GEOSMessageHandler) __errorHandler);
#endif
}


static void __warningIgnore(const char *fmt, ...) {
	return;
}

GEOSContextHandle_t geos_init2(void) {

#ifdef HAVE350
	GEOSContextHandle_t ctxt = GEOS_init_r();
	GEOSContext_setNoticeHandler_r(ctxt, __warningIgnore);
	GEOSContext_setErrorHandler_r(ctxt, __errorHandler);
	return ctxt;
#else
	return initGEOS_r((GEOSMessageHandler) __warningIgnore, (GEOSMessageHandler) __errorHandler);
#endif
}





GEOSGeometry* geos_line(const std::vector<double> &x, const std::vector<double> &y, GEOSContextHandle_t hGEOSCtxt) {
	GEOSCoordSequence *pseq;
	size_t n = x.size();
	pseq = GEOSCoordSeq_create_r(hGEOSCtxt, n, 2);
	for (size_t i = 0; i < n; i++) {
		GEOSCoordSeq_setX_r(hGEOSCtxt, pseq, i, x[i]);
		GEOSCoordSeq_setY_r(hGEOSCtxt, pseq, i, y[i]);
	}
	GEOSGeometry* g = GEOSGeom_createLineString_r(hGEOSCtxt, pseq);
	// GEOSCoordSeq_destroy(pseq); 
	return g;
}



GEOSGeometry* geos_linearRing(const std::vector<double> &x, const std::vector<double> &y, GEOSContextHandle_t hGEOSCtxt) {
	GEOSCoordSequence *pseq;
	size_t n = x.size();
	pseq = GEOSCoordSeq_create_r(hGEOSCtxt, n, 2);
	for (size_t i = 0; i < n; i++) {
		GEOSCoordSeq_setX_r(hGEOSCtxt, pseq, i, x[i]);
		GEOSCoordSeq_setY_r(hGEOSCtxt, pseq, i, y[i]);
	}
	GEOSGeometry* g = GEOSGeom_createLinearRing_r(hGEOSCtxt, pseq);
	// GEOSCoordSeq_destroy(pseq); 
	return g;
}


GEOSGeometry* geos_polygon(const std::vector<double> &x, const std::vector<double> &y, std::vector<std::vector<double>> &hx, std::vector<std::vector<double>> &hy, GEOSContextHandle_t hGEOSCtxt) {
	GEOSGeometry* shell = geos_linearRing(x, y, hGEOSCtxt);
	size_t nh = hx.size();
	std::vector<GEOSGeometry*> holes(nh);
	for (size_t i=0; i<nh; i++) {
		holes[i] = geos_linearRing(hx[i], hy[i], hGEOSCtxt);
	}
	GEOSGeometry* g = GEOSGeom_createPolygon_r(hGEOSCtxt, shell, &holes[0], nh);
	return g;
}



void getHoles(SpatPart &p, std::vector<std::vector<double>> &hx, std::vector<std::vector<double>> &hy) {
	size_t nh = p.nHoles();
	hx.resize(0);
	hy.resize(0);
	hx.reserve(nh);
	hy.reserve(nh);
	if (nh == 0) return;	
	for (size_t i=0; i<nh; i++) {
		SpatHole h = p.getHole(i);
		hx.push_back(h.x); 
		hy.push_back(h.y); 
	}
	return;
}
		

std::vector<GeomPtr> geos_geoms(SpatVector *v, GEOSContextHandle_t hGEOSCtxt) {
	size_t n = v->size();
	std::vector<GeomPtr> g;
	g.reserve(n);
	std::string vt = v->type();
	if (vt == "points") {
		for (size_t i=0; i<n; i++) {
			SpatGeom svg = v->getGeom(i);
			size_t np = svg.size();
			GEOSCoordSequence *pseq;
			std::vector<GEOSGeometry*> geoms;
			geoms.reserve(np);
			for (size_t j = 0; j < np; j++) {
				//SpatPart svp = svg.getPart(j);
				pseq = GEOSCoordSeq_create_r(hGEOSCtxt, 1, 2);
				GEOSCoordSeq_setX_r(hGEOSCtxt, pseq, 0, svg.parts[j].x[0]);
				GEOSCoordSeq_setY_r(hGEOSCtxt, pseq, 0, svg.parts[j].y[0]);
				GEOSGeometry* pt = GEOSGeom_createPoint_r(hGEOSCtxt, pseq);
				if (pt != NULL) {
					geoms.push_back(pt);
				}
			}
			GEOSGeometry* gcol = (np == 1) ? geoms[0] :		
				GEOSGeom_createCollection_r(hGEOSCtxt, GEOS_MULTIPOINT, &geoms[0], np);			
			g.push_back( geos_ptr(gcol, hGEOSCtxt) );
		}

	} else if (vt == "lines") {
		// gp = NULL;
		for (size_t i=0; i<n; i++) {
			SpatGeom svg = v->getGeom(i);
			size_t np = svg.size();
			std::vector<GEOSGeometry*> geoms;
			geoms.reserve(np);
			for (size_t j=0; j < np; j++) {
				//SpatPart svp = svg.getPart(j);			
				GEOSGeometry* gp = geos_line(svg.parts[j].x, svg.parts[j].y, hGEOSCtxt); 
				if (gp != NULL) {
					geoms.push_back(gp);
				}
			}		
			GEOSGeometry* gcol = (np == 1) ? geoms[0] :		
				GEOSGeom_createCollection_r(hGEOSCtxt, GEOS_MULTILINESTRING, &geoms[0], np);			
			g.push_back( geos_ptr(gcol, hGEOSCtxt) );
		}

	} else { // polygons

		std::vector<std::vector<double>> hx, hy;
		for (size_t i=0; i<n; i++) {
			SpatGeom svg = v->getGeom(i);
			size_t np = svg.size();
			std::vector<GEOSGeometry*> geoms(np);
			for (size_t j=0; j < np; j++) {
				SpatPart svp = svg.getPart(j);
				getHoles(svp, hx, hy); 
				geoms[j] = geos_polygon(svp.x, svp.y, hx, hy, hGEOSCtxt); 
			}
			//Rcpp::Rcout << np << std::endl;
			GEOSGeometry* gcol = (np == 1) ? geoms[0] :		
				GEOSGeom_createCollection_r(hGEOSCtxt, GEOS_MULTIPOLYGON, &geoms[0], np);			
			g.push_back( geos_ptr(gcol, hGEOSCtxt));
		}
	}
	return g;
}



SpatVector vect_from_geos(std::vector<GeomPtr> &geoms , GEOSContextHandle_t hGEOSCtxt, std::string vt) {

	SpatVector out;
	SpatVector v;

	size_t ng = geoms.size();
	std::vector<unsigned> gid, gp, hole;
	std::vector<double> x, y;
	bool xok, yok;

	if ((vt == "points") | (vt == "lines")) {	
		for(size_t i = 0; i < ng; i++) {
			GEOSGeometry* g = geoms[i].get();
			size_t np = GEOSGetNumGeometries_r(hGEOSCtxt, g);
			for(size_t j = 0; j<np; j++) {
				const GEOSGeometry* part = GEOSGetGeometryN_r(hGEOSCtxt, g, j);
				const GEOSCoordSequence* crds = GEOSGeom_getCoordSeq_r(hGEOSCtxt, part); 		
				int npts = -1;
//				if (vt == "points") {	
				npts = GEOSGetNumCoordinates_r(hGEOSCtxt, part);
//				} else if (vt == "lines") {
//					npts = GEOSGeomGetNumPoints_r(hGEOSCtxt, part); // for lines
//				}		
				if (npts < 0) {
					out.setError("GEOS exception 9");
					return out;
				}
				double xvalue = 0;
				double yvalue = 0;
				for (int p=0; p < npts; p++) {
					xok = GEOSCoordSeq_getX_r(hGEOSCtxt, crds, p, &xvalue);				
					yok = GEOSCoordSeq_getY_r(hGEOSCtxt, crds, p, &yvalue);				
					if (xok & yok) {
						x.push_back(xvalue);
						y.push_back(yvalue);
						gid.push_back(i);			
						gp.push_back(j);			
						hole.push_back(0);
					}
				}
			}
		}
	} else { // polygons
		for(size_t i = 0; i < ng; i++) {
			GEOSGeometry* g = geoms[i].get();
			size_t np = GEOSGetNumGeometries_r(hGEOSCtxt, g);
			
			for(size_t j = 0; j<np; j++) {
				
				const GEOSGeometry* part = GEOSGetGeometryN_r(hGEOSCtxt, g, j);
				const GEOSGeometry* ring = GEOSGetExteriorRing_r(hGEOSCtxt, part);
				const GEOSCoordSequence* crds = GEOSGeom_getCoordSeq_r(hGEOSCtxt, ring); 		
				int npts = -1;
				npts = GEOSGetNumCoordinates_r(hGEOSCtxt, ring);
				if (npts < 0) {
					Rcpp::Rcout << "exception 99" << std::endl;
					continue;
				}
				double xvalue = 0;
				double yvalue = 0;
				for (int p=0; p < npts; p++) {
					xok = GEOSCoordSeq_getX_r(hGEOSCtxt, crds, p, &xvalue);				
					yok = GEOSCoordSeq_getY_r(hGEOSCtxt, crds, p, &yvalue);				
					if (xok & yok) {
						x.push_back(xvalue);
						y.push_back(yvalue);
						gid.push_back(i);			
						gp.push_back(j);			
						hole.push_back(0);
					}
				}
				int nholes = GEOSGetNumInteriorRings_r(hGEOSCtxt, part);
				for (int h=0; h < nholes; h++) {
					const GEOSGeometry* ring = GEOSGetInteriorRingN_r(hGEOSCtxt, part, h);

					const GEOSCoordSequence* crds = GEOSGeom_getCoordSeq_r(hGEOSCtxt, ring); 		
					int npts = -1;
					npts = GEOSGetNumCoordinates_r(hGEOSCtxt, ring);
					if (npts < 0) {
						Rcpp::Rcout << "exception 99" << std::endl;
						continue;
					}
					double xvalue = 0;
					double yvalue = 0;
					for (int p=0; p < npts; p++) {
						xok = GEOSCoordSeq_getX_r(hGEOSCtxt, crds, p, &xvalue);			
						yok = GEOSCoordSeq_getY_r(hGEOSCtxt, crds, p, &yvalue);				
						if (xok & yok) {
							x.push_back(xvalue);
							y.push_back(yvalue);
							gid.push_back(i);
							gp.push_back(j);
							hole.push_back(h);
						}
					}
				}
			}
		}	
	}	
	v.setGeometry(vt, gid, gp, x, y, hole);
	return v;
}



bool pointsFromGeom(GEOSContextHandle_t hGEOSCtxt, const GEOSGeometry* part, 
const unsigned i, const unsigned j, std::vector<double> &x, std::vector<double> &y, 
std::vector<unsigned> &gid, std::vector<unsigned> &gp, std::vector<unsigned> &hole, std::string &msg) {
	const GEOSCoordSequence* crds = GEOSGeom_getCoordSeq_r(hGEOSCtxt, part); 		
	int npts = -1;
	npts = GEOSGetNumCoordinates_r(hGEOSCtxt, part);
	if (npts < 0) {
		msg = "GEOS exception 9";
		return false;
	}
	double xvalue = 0;
	double yvalue = 0;
	for (int p=0; p < npts; p++) {
		bool xok = GEOSCoordSeq_getX_r(hGEOSCtxt, crds, p, &xvalue);
		bool yok = GEOSCoordSeq_getY_r(hGEOSCtxt, crds, p, &yvalue);	
		if (xok & yok) {
			x.push_back(xvalue);
			y.push_back(yvalue);
			gid.push_back(i);			
			gp.push_back(j);			
			hole.push_back(0);
		}
	}
	return true;
}


bool polysFromGeom(GEOSContextHandle_t hGEOSCtxt, const GEOSGeometry* part, 
const unsigned i, const unsigned j, std::vector<double> &x, std::vector<double> &y, 
std::vector<unsigned> &gid, std::vector<unsigned> &gp, std::vector<unsigned> &hole, std::string &msg) {
	const GEOSGeometry* ring = GEOSGetExteriorRing_r(hGEOSCtxt, part);
	const GEOSCoordSequence* crds = GEOSGeom_getCoordSeq_r(hGEOSCtxt, ring); 		
	int npts = -1;
	npts = GEOSGetNumCoordinates_r(hGEOSCtxt, ring);
	if (npts < 0) {
		msg = "exception 99";
		return false;
	}
	double xvalue = 0;
	double yvalue = 0;
	for (int p=0; p < npts; p++) {
		bool xok = GEOSCoordSeq_getX_r(hGEOSCtxt, crds, p, &xvalue);
		bool yok = GEOSCoordSeq_getY_r(hGEOSCtxt, crds, p, &yvalue);
		if (xok & yok) {
			x.push_back(xvalue);
			y.push_back(yvalue);
			gid.push_back(i);			
			gp.push_back(j);			
			hole.push_back(0);
		}
	}
	int nholes = GEOSGetNumInteriorRings_r(hGEOSCtxt, part);
	for (int h=0; h < nholes; h++) {
		const GEOSGeometry* ring = GEOSGetInteriorRingN_r(hGEOSCtxt, part, h);
		const GEOSCoordSequence* crds = GEOSGeom_getCoordSeq_r(hGEOSCtxt, ring); 		
		int npts = -1;
		npts = GEOSGetNumCoordinates_r(hGEOSCtxt, ring);
		if (npts < 0) {
			msg  = "exception 123";
			return false;
		}
		double xvalue = 0;
		double yvalue = 0;
		for (int p=0; p < npts; p++) {
			bool xok = GEOSCoordSeq_getX_r(hGEOSCtxt, crds, p, &xvalue);	
			bool yok = GEOSCoordSeq_getY_r(hGEOSCtxt, crds, p, &yvalue);	
			if (xok & yok) {
				x.push_back(xvalue);
				y.push_back(yvalue);
				gid.push_back(i);			
				gp.push_back(j);			
				hole.push_back(h);
			}
		}
	}
	return true;
}


SpatVectorCollection coll_from_geos(std::vector<GeomPtr> &geoms , GEOSContextHandle_t hGEOSCtxt) {

	SpatVectorCollection out;

	size_t ng = geoms.size();
	std::vector<unsigned> pt_gid, pt_gp, pt_hole;
	std::vector<unsigned> ln_gid, ln_gp, ln_hole;
	std::vector<unsigned> pl_gid, pl_gp, pl_hole;
	std::vector<double> pt_x, pt_y, ln_x, ln_y, pl_x, pl_y;

	std::string msg;
	//Rcpp::Rcout << ng << " geoms" << std::endl;
	for(size_t i = 0; i < ng; i++) {
		GEOSGeometry* g = geoms[i].get();
		std::string gt = GEOSGeomType_r(hGEOSCtxt, g);
		size_t np = GEOSGetNumGeometries_r(hGEOSCtxt, g);

		//Rcpp::Rcout << gt << std::endl;
		//Rcpp::Rcout << np << " parts" << std::endl;
		
		if (gt == "Point" || gt == "MultiPoint") {
			for(size_t j = 0; j<np; j++) {
				const GEOSGeometry* part = GEOSGetGeometryN_r(hGEOSCtxt, g, j);

				if (!pointsFromGeom(hGEOSCtxt, part, i, j, pt_x, pt_y, pt_gid, pt_gp, pt_hole, msg)) {
					out.setError(msg);
					return out;
				}
			}	
		} else if (gt == "LineString" || gt == "MultiLineString") {
			for(size_t j = 0; j<np; j++) {
				const GEOSGeometry* part = GEOSGetGeometryN_r(hGEOSCtxt, g, j);
				if (!pointsFromGeom(hGEOSCtxt, part, i, j, ln_x, ln_y, ln_gid, ln_gp, ln_hole, msg)) {
					out.setError(msg);
					return out;
				}
			}
		} else if (gt == "Polygon" || gt == "MultiPolygon") {
			for(size_t j = 0; j<np; j++) {
				const GEOSGeometry* part = GEOSGetGeometryN_r(hGEOSCtxt, g, j);
				if (!polysFromGeom(hGEOSCtxt, part, i, j, pl_x, pl_y, pl_gid, pl_gp, pl_hole, msg)) {
					out.setError(msg);
					return out;
				}
			}
		} else if (gt == "GeometryCollection") {
			//Rcpp::Rcout << GEOSGeom_getDimensions_r(hGEOSCtxt, g) << std::endl;
			for(size_t j = 0; j<np; j++) {
				const GEOSGeometry* gg = GEOSGetGeometryN_r(hGEOSCtxt, g, j);
				std::string ggt = GEOSGeomType_r(hGEOSCtxt, gg);
				const GEOSGeometry* part = GEOSGetGeometryN_r(hGEOSCtxt, gg, j);
				if (ggt == "Polygon" || ggt == "MultiPolygon") {
					if (!polysFromGeom(hGEOSCtxt, part, i, j, pl_x, pl_y, pl_gid, pl_gp, pl_hole, msg)) {
						out.setError(msg);
						return out;
					}
				} else if (ggt == "Point" || ggt == "MultiPoint") {
					if (!polysFromGeom(hGEOSCtxt, part, i, j, pt_x, pt_y, pt_gid, pt_gp, pt_hole, msg)) {
						out.setError(msg);
						return out;
					}
				} else if (ggt == "Line" || ggt == "MultiLine") {
					if (!polysFromGeom(hGEOSCtxt, part, i, j, pl_x, pl_y, pl_gid, pl_gp, pl_hole, msg)) {
						out.setError(msg);
						return out;
					}
				} else {
					out.addWarning("unhandeled Collection geom: " + ggt);
				}
			}
		} else {
			out.setError("what is this: " + gt + "?");
		}
	}

	if (pl_x.size() > 0) {
		SpatVector v;
		v.setGeometry("polygons", pl_gid, pl_gp, pl_x, pl_y, pl_hole);
		out.push_back(v);
		//Rcpp::Rcout << "pls" << std::endl;
	}
	if (ln_x.size() > 0) {
		SpatVector v;
		v.setGeometry("lines", ln_gid, ln_gp, ln_x, ln_y, ln_hole);
		out.push_back(v);
		//Rcpp::Rcout << "lns" << std::endl;
	}
	if (pt_x.size() > 0) {
		SpatVector v;
		v.setGeometry("points", pt_gid, pt_gp, pt_x, pt_y, pt_hole);
		out.push_back(v);
		//Rcpp::Rcout << "pts" << std::endl;
	}
	return out;
}

