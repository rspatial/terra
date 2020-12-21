//#define GEOS_USE_ONLY_R_API
//#include <geos_c.h>

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

//#include "geos_spat.h"
#include "spatVector.h"
#include <cstdarg> 
#include <cstring> 
#include <memory>

#include "Rcpp.h"


static GeomPtr geos_ptr(GEOSGeometry* g, GEOSContextHandle_t hGEOSctxt) {
	auto deleter = std::bind(GEOSGeom_destroy_r, hGEOSctxt, std::placeholders::_1);
	return GeomPtr(g, deleter);
}


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


GEOSContextHandle_t CPL_geos_init(void) {
#ifdef HAVE350
	GEOSContextHandle_t ctxt = GEOS_init_r();
	GEOSContext_setNoticeHandler_r(ctxt, __warningHandler);
	GEOSContext_setErrorHandler_r(ctxt, __errorHandler);
	return ctxt;
#else
	return initGEOS_r((GEOSMessageHandler) __warningHandler, (GEOSMessageHandler) __errorHandler);
#endif
}



GEOSGeometry* geos_line(std::vector<double> x, std::vector<double> y, GEOSContextHandle_t hGEOSCtxt) {
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



GEOSGeometry* geos_linearRing(std::vector<double> x, std::vector<double> y, GEOSContextHandle_t hGEOSCtxt) {
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


GEOSGeometry* geos_polygon(std::vector<double> &x, std::vector<double> &y, std::vector<std::vector<double>> &hx, std::vector<std::vector<double>> &hy, GEOSContextHandle_t hGEOSCtxt) {
	GEOSGeometry* shell = geos_linearRing(x, y, hGEOSCtxt);
	size_t nh = hx.size();
	Rcpp::Rcout << "holes: "<< nh << std::endl;
	std::vector<GEOSGeometry*> holes(nh);
	for (size_t i=0; i<nh; i++) {
		holes[i] = geos_linearRing(hx[i], hy[i], hGEOSCtxt);
	}
	GEOSGeometry* g = GEOSGeom_createPolygon_r(hGEOSCtxt, shell, &holes[0], nh);
	return g;
}



void getHoles(SpatPart &p, std::vector<std::vector<double>> &hx, std::vector<std::vector<double>> &hy) {
	size_t nh = p.nHoles();
	hx = std::vector<std::vector<double>>(nh,std::vector<double>(0));	
	hy = std::vector<std::vector<double>>(nh,std::vector<double>(0));	
	if (nh == 0) return;	
	for (size_t i=0; i<nh; i++) {
		SpatHole h = p.getHole(i);
		hx[i] = h.x; 
		hy[i] = h.y;
	}
	return;
}
		

std::vector<GeomPtr> SpatVector::geos_geoms(GEOSContextHandle_t hGEOSCtxt) {
	size_t n = size();
	std::vector<GeomPtr> g;
	g.reserve(n);
	std::string vt = type();
	if (vt == "points") {
		std::vector<std::vector<double>> xy = coordinates();
		std::vector<double> x = xy[0];
		std::vector<double> y = xy[1];
		GEOSCoordSequence *pseq;
		for (size_t i = 0; i < n; i++) {
			pseq = GEOSCoordSeq_create_r(hGEOSCtxt, 1, 2);
			GEOSCoordSeq_setX_r(hGEOSCtxt, pseq, 0, x[i]);
			GEOSCoordSeq_setY_r(hGEOSCtxt, pseq, 0, y[i]);
			GEOSGeometry* pt = GEOSGeom_createPoint_r(hGEOSCtxt, pseq);
			g.push_back( geos_ptr(pt, hGEOSCtxt) );
			// GEOSCoordSeq_destroy(pseq); 
		}
	} else if (vt == "lines") {
		// gp = NULL;
		for (size_t i=0; i<n; i++) {
			SpatGeom svg = getGeom(i);
			size_t np = svg.size();
			std::vector<GEOSGeometry*> geoms;
			geoms.reserve(np);
			for (size_t j=0; j < np; j++) {
				SpatPart svp = svg.getPart(j);			
				GEOSGeometry* gp = geos_line(svp.x, svp.y, hGEOSCtxt); 
				if (gp != NULL) {
					geoms.push_back(gp);
				}
			}		
			GEOSGeometry* gcol = (np == 1) ? geoms[0] :		
				GEOSGeom_createCollection_r(hGEOSCtxt, GEOS_GEOMETRYCOLLECTION, &geoms[0], np);			
			g.push_back( geos_ptr(gcol, hGEOSCtxt) );
		}

	} else { // polygons

		std::vector<std::vector<double>> hx, hy;
		for (size_t i=0; i<n; i++) {
			SpatGeom svg = getGeom(i);
			size_t np = svg.size();
			std::vector<GEOSGeometry*> geoms(np);
			for (size_t j=0; j < np; j++) {
				SpatPart svp = svg.getPart(j);
				getHoles(svp, hx, hy); 
				geoms[j] = geos_polygon(svp.x, svp.y, hx, hy, hGEOSCtxt); 
			}
			GEOSGeometry* gcol = (np == 1) ? geoms[0] :		
				GEOSGeom_createCollection_r(hGEOSCtxt, GEOS_GEOMETRYCOLLECTION, &geoms[0], np);			
			g.push_back( geos_ptr(gcol, hGEOSCtxt));
		}
	}
	return g;
}



SpatVectorCollection vect_from_geos(std::vector<GeomPtr> &geoms , GEOSContextHandle_t hGEOSCtxt, std::string vt) {

	SpatVectorCollection out;
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
				Rcpp::Rcout << "holes: " << nholes << std::endl;
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
							hole.push_back(1);
						}
					}
				}
			}
		}	
	}	
	v.setGeometry(vt, gid, gp, x, y, hole);
	out.push_back(v);
	return out;
}



SpatVector SpatVector::allerretour() {
	GEOSContextHandle_t hGEOSCtxt = CPL_geos_init();
	std::vector<GeomPtr> g = geos_geoms(hGEOSCtxt);
	SpatVectorCollection out = vect_from_geos(g, hGEOSCtxt, type());
	geos_finish(hGEOSCtxt);
	return out.get(0);
}


SpatVector SpatVector::buffer2(double dist, unsigned nQuadSegs, unsigned capstyle) {

	std::string vt = type();
	if ((vt == "points") && (dist <= 0)) {
		dist = -dist;
	}

	GEOSContextHandle_t hGEOSCtxt = CPL_geos_init();
	SpatVector out;
	SpatVector f = remove_holes();

	std::vector<GeomPtr> g = geos_geoms(hGEOSCtxt);
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
	SpatVectorCollection coll = vect_from_geos(b, hGEOSCtxt, "polygons");

//	out = spat_from_geom(hGEOSCtxt, g, "points");
	geos_finish(hGEOSCtxt);
	return coll.get(0);	
}




SpatVector SpatVector::intersect(SpatVector v) {

	SpatVector out;

	GEOSContextHandle_t hGEOSCtxt = CPL_geos_init();
	std::vector<GeomPtr> x = geos_geoms(hGEOSCtxt);
	std::vector<GeomPtr> y = geos_geoms(hGEOSCtxt);	
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
			unsigned check = 0; //GEOSisEmpty_r(hGEOSCtxt, geom);
			if (check > 1) {
				out.setError("GEOS exception");
				geos_finish(hGEOSCtxt);
				return(out);			
			} else if (check == 1) {
				result.push_back(geos_ptr(geom, hGEOSCtxt));
				atts[0].push_back(i);
				atts[1].push_back(j);
			}
		}
	}

	// this should happen in the loop. And different types are possible!
	SpatVectorCollection coll = vect_from_geos(result, hGEOSCtxt, type());
	// deal with attributes
	geos_finish(hGEOSCtxt);
	return coll.get(0);	
}




SpatVector SpatVector::centroid() {

	SpatVector out;

	GEOSContextHandle_t hGEOSCtxt = CPL_geos_init();
	std::vector<GeomPtr> g = geos_geoms(hGEOSCtxt);
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
	SpatVectorCollection coll = vect_from_geos(b, hGEOSCtxt, "points");
	geos_finish(hGEOSCtxt);

	out = coll.get(0);
	out.df = df;
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

