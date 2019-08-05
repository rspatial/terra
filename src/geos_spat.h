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


#include <cstdarg> 
#include <cstring> 
#include <memory>
#include <functional>

#include "spatVector.h"



static void __errorHandler(const char *fmt, ...) { // #nocov start
	char buf[BUFSIZ], *p;
	va_list ap;
	va_start(ap, fmt);
	vsprintf(buf, fmt, ap);
	va_end(ap);
	p = buf + strlen(buf) - 1;
	if(strlen(buf) > 0 && *p == '\n') *p = '\0';
//	error(buf);
	return; // #nocov end
} 

static void __warningHandler(const char *fmt, ...) {
	char buf[BUFSIZ], *p;
	va_list ap;
	va_start(ap, fmt);
	vsprintf(buf, fmt, ap);
	va_end(ap);
	p = buf + strlen(buf) - 1;
	if(strlen(buf) > 0 && *p == '\n') *p = '\0';
//	warning(buf);
	return;
}

void CPL_geos_finish(GEOSContextHandle_t ctxt) {
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


using GeomPtr = std::unique_ptr<GEOSGeometry, std::function<void(GEOSGeometry*)> >;

static GeomPtr geos_ptr(GEOSGeometry* g, GEOSContextHandle_t hGEOSctxt) {
	auto deleter = std::bind(GEOSGeom_destroy_r, hGEOSctxt, std::placeholders::_1);
	return GeomPtr(g, deleter);
}

std::vector<GeomPtr> make_geos_points(std::vector<double> x, std::vector<double> y, GEOSContextHandle_t hGEOSCtxt) {
	GEOSCoordSequence *pseq;
	size_t n = x.size();
	std::vector<GeomPtr> g(n);
	for (size_t i = 0; i < n; i++) {
		pseq = GEOSCoordSeq_create_r(hGEOSCtxt, 1, 2);
		GEOSCoordSeq_setX_r(hGEOSCtxt, pseq, 0, x[i]);
		GEOSCoordSeq_setY_r(hGEOSCtxt, pseq, 0, y[i]);
		GEOSGeometry* pt = GEOSGeom_createPoint_r(hGEOSCtxt, pseq);
		g[i] = geos_ptr(pt, hGEOSCtxt);
		// GEOSCoordSeq_destroy(pseq); 
	}
	return g;
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
	std::vector<GEOSGeometry*> holes(nh);
	for (size_t i=0; i<nh; i++) {
		holes[i] = geos_linearRing(hx[i], hy[i], hGEOSCtxt);
	}
	GEOSGeometry* g = GEOSGeom_createPolygon_r(hGEOSCtxt, shell, &holes[0], nh);
	return g;
}



std::vector<GeomPtr> make_geos_lines(SpatVector *v, GEOSContextHandle_t hGEOSCtxt) {
	size_t n = v->size();
	std::vector<GeomPtr> g(n);
	GEOSGeometry* gp = NULL;
	for (size_t i=0; i<n; i++) {
		SpatGeom svg = v->getGeom(i);
		size_t np = svg.size();
		std::vector<GEOSGeometry*> geoms(np);
		for (size_t j=0; j < np; j++) {
			SpatPart svp = svg.getPart(j);			
			gp = geos_line(svp.x, svp.y, hGEOSCtxt); 
			if (gp != NULL) {
				geoms[j] = gp;
			}
		}		
        GEOSGeometry* gcol = (np == 1) ? geoms[0] :		
			GEOSGeom_createCollection_r(hGEOSCtxt, GEOS_GEOMETRYCOLLECTION, &geoms[0], np);			
    	g[i] = geos_ptr(gcol, hGEOSCtxt);
	}
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
		


std::vector<GeomPtr> make_geos_polygons(SpatVector *v, GEOSContextHandle_t hGEOSCtxt) {
	size_t n = v->size();
	std::vector<GeomPtr> g(n);
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
		GEOSGeometry* gcol = (np == 1) ? geoms[0] :		
			GEOSGeom_createCollection_r(hGEOSCtxt, GEOS_GEOMETRYCOLLECTION, &geoms[0], np);			
    	g[i] = geos_ptr(gcol, hGEOSCtxt);
	}
	return g;
}



std::vector<GeomPtr> geom_from_spat(SpatVector* v, GEOSContextHandle_t hGEOSCtxt) {
	std::vector<GeomPtr> g;
	std::string vt = v->type();
	if (vt == "points") {
		std::vector<std::vector<double>> xy = v->coordinates();
		g = make_geos_points(xy[0], xy[1], hGEOSCtxt);
	} else if (vt == "lines") {
		g = make_geos_polygons(v, hGEOSCtxt);
	} else { // polygons
		g = make_geos_polygons(v, hGEOSCtxt);
	}
	return g;
}

	


SpatVector spat_from_geom(GEOSContextHandle_t hGEOSCtxt, std::vector<GeomPtr> & geoms, std::string vt) {
	SpatVector out;

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
				const GEOSGeometry* shell = GEOSGetExteriorRing_r(hGEOSCtxt, part);
				const GEOSCoordSequence* crds = GEOSGeom_getCoordSeq_r(hGEOSCtxt, shell); 		
				int npts = -1;
				npts = GEOSGetNumCoordinates_r(hGEOSCtxt, part);
				if (npts < 0) {
					out.setError("GEOS exception 99");
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
				//int nholes = GEOSGetNumInteriorRings_r(hGEOSCtxt, part);
				

			}
		}	
	}
	out.setGeometry(vt, gid, gp, x, y, hole);
	return out;
}


