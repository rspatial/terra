// Copyright (c) 2018-2023  Robert J. Hijmans
//
// This file is part of the "spat" library.
//
// spat is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// spat is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with spat. If not, see <http://www.gnu.org/licenses/>.

#define GEOS_USE_ONLY_R_API
#include <geos_c.h>

#if GEOS_VERSION_MAJOR == 3
# if GEOS_VERSION_MINOR >= 5
#  define GEOS350
# endif
# if GEOS_VERSION_MINOR == 6
#  if GEOS_VERSION_PATCH >= 1
#   define GEOS361
#  endif
# endif
# if GEOS_VERSION_MINOR >= 7
#  define GEOS361
#  define GEOS370
# endif
# if GEOS_VERSION_MINOR >= 8
#  define GEOS361
#  define GEOS370
#  define GEOS380
# endif
#else
# if GEOS_VERSION_MAJOR > 3
#  define GEOS350
#  define GEOS370
#  define GEOS361
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

using PrepGeomPtr= std::unique_ptr<const GEOSPreparedGeometry, std::function<void(const GEOSPreparedGeometry*)> >;
static PrepGeomPtr geos_ptr(const GEOSPreparedGeometry* pg, GEOSContextHandle_t hGEOSctxt) {
	auto deleter = std::bind(GEOSPreparedGeom_destroy_r, hGEOSctxt, std::placeholders::_1);
	return PrepGeomPtr(pg, deleter);
}

using TreePtr= std::unique_ptr<GEOSSTRtree, std::function<void(GEOSSTRtree*)> >;
static TreePtr geos_ptr(GEOSSTRtree* t, GEOSContextHandle_t hGEOSctxt) {
	auto deleter = std::bind(GEOSSTRtree_destroy_r, hGEOSctxt, std::placeholders::_1);
	return TreePtr(t, deleter);
}

/*
using wkbPtr = std::unique_ptr<unsigned char, std::function<void(unsigned char*)> >;
static wkbPtr geos_wkb(unsigned char wkb, GEOSContextHandle_t hGEOSctxt) {
	auto deleter = std::bind(GEOSFree_r, hGEOSctxt, std::placeholders::_1);
	return wkbPtr(wkb, deleter);
}
*/

#ifdef useRcpp
#include "Rcpp.h"

template <typename... Args>
inline void warnNoCall(const char* fmt, Args&&... args ) {
    Rf_warningcall(R_NilValue, "%s", tfm::format(fmt, std::forward<Args>(args)... ).c_str());
}

template <typename... Args>
inline void NORET errNoCall(const char* fmt, Args&&... args) {
    throw Rcpp::exception(tfm::format(fmt, std::forward<Args>(args)... ).c_str(), false);
}


static void __errorHandler(const char *fmt, ...) { 
	char buf[BUFSIZ], *p;
	va_list ap;
	va_start(ap, fmt);
	size_t n = BUFSIZ;
	vsnprintf(buf, n, fmt, ap);
//	vsprintf(buf, fmt, ap);
	va_end(ap);
	p = buf + strlen(buf) - 1;
	if (strlen(buf) > 0 && *p == '\n') *p = '\0';
    errNoCall(buf); 
	return; 
} 

static void __warningHandler(const char *fmt, ...) {
	char buf[BUFSIZ], *p;
	va_list ap;
	va_start(ap, fmt);
	size_t n = BUFSIZ;
	vsnprintf(buf, n, fmt, ap);
//	vsprintf(buf, fmt, ap);
	va_end(ap);
	p = buf + strlen(buf) - 1;
	if (strlen(buf) > 0 && *p == '\n') *p = '\0';
    warnNoCall(buf); 
	return;
}

static void __checkInterruptFn(void*) {
	R_CheckUserInterrupt();
}

static void __checkInterrupt() {
	// Adapted (in sf) from Rcpp/Interrupt.h
	if (!R_ToplevelExec(__checkInterruptFn, nullptr)) {
		GEOS_interruptRequest();
	}
}

inline GEOSContextHandle_t geos_init(void) {
#ifdef GEOS350
	GEOSContextHandle_t ctxt = GEOS_init_r();
	GEOSContext_setNoticeHandler_r(ctxt, __warningHandler);
	GEOSContext_setErrorHandler_r(ctxt, __errorHandler);
	GEOS_interruptRegisterCallback(__checkInterrupt);
	return ctxt;
#else
	return initGEOS_r((GEOSMessageHandler) __warningHandler, (GEOSMessageHandler) __errorHandler);
#endif
}


#else 

#include <iostream>
static void __errorHandler(const char *fmt, ...) { 
	char buf[BUFSIZ], *p;
	va_list ap;
	va_start(ap, fmt);
	size_t n = BUFSIZ;
	vsnprintf(buf, n, fmt, ap);
//	vsprintf(buf, fmt, ap);
	va_end(ap);
	p = buf + strlen(buf) - 1;
	if(strlen(buf) > 0 && *p == '\n') *p = '\0';
    std::cout << buf << std::endl; 
	return; 
} 

static void __warningHandler(const char *fmt, ...) {
	char buf[BUFSIZ], *p;
	va_list ap;
	va_start(ap, fmt);
	size_t n = BUFSIZ;
	vsnprintf(buf, n, fmt, ap);
//	vsprintf(buf, fmt, ap);
	va_end(ap);
	p = buf + strlen(buf) - 1;
	if(strlen(buf) > 0 && *p == '\n') *p = '\0';
    std::cout << buf << std::endl; 
	return;
}

inline GEOSContextHandle_t geos_init(void) {
#ifdef GEOS350
	GEOSContextHandle_t ctxt = GEOS_init_r();
	GEOSContext_setNoticeHandler_r(ctxt, __warningHandler);
	GEOSContext_setErrorHandler_r(ctxt, __errorHandler);
	return ctxt;
#else
	return initGEOS_r((GEOSMessageHandler) __warningHandler, (GEOSMessageHandler) __errorHandler);
#endif
}


#endif 


inline void geos_finish(GEOSContextHandle_t ctxt) {
#ifdef GEOS350
	GEOS_finish_r(ctxt);
#else
	finishGEOS_r(ctxt);
#endif
}




static void __warningIgnore(const char *fmt, ...) {
	return;
}

inline GEOSContextHandle_t geos_init2(void) {

#ifdef GEOS350
	GEOSContextHandle_t ctxt = GEOS_init_r();
	GEOSContext_setNoticeHandler_r(ctxt, __warningIgnore);
	GEOSContext_setErrorHandler_r(ctxt, __errorHandler);
	return ctxt;
#else
	return initGEOS_r((GEOSMessageHandler) __warningIgnore, (GEOSMessageHandler) __errorHandler);
#endif
}



inline GEOSGeometry* geos_line(const std::vector<double> &x, const std::vector<double> &y, GEOSContextHandle_t hGEOSCtxt) {
	GEOSCoordSequence *pseq;
	size_t n = x.size();
	if (n < 2) {
		pseq = GEOSCoordSeq_create_r(hGEOSCtxt, 0, 2);
		GEOSGeometry* g = GEOSGeom_createLineString_r(hGEOSCtxt, pseq);
		return g;	
	}
	pseq = GEOSCoordSeq_create_r(hGEOSCtxt, n, 2);
	for (size_t i = 0; i < n; i++) {
		GEOSCoordSeq_setX_r(hGEOSCtxt, pseq, i, x[i]);
		GEOSCoordSeq_setY_r(hGEOSCtxt, pseq, i, y[i]);
	}
	GEOSGeometry* g = GEOSGeom_createLineString_r(hGEOSCtxt, pseq);
	return g;
}



inline GEOSGeometry* geos_linearRing(const std::vector<double> &x, const std::vector<double> &y, GEOSContextHandle_t hGEOSCtxt) {
	GEOSCoordSequence *pseq;
	size_t n = x.size();
	if (n < 3) {
		pseq = GEOSCoordSeq_create_r(hGEOSCtxt, 0, 2);
		GEOSGeometry* g = GEOSGeom_createLinearRing_r(hGEOSCtxt, pseq);
		return g;	
	}
	pseq = GEOSCoordSeq_create_r(hGEOSCtxt, n, 2);
	for (size_t i = 0; i < n; i++) {
		GEOSCoordSeq_setX_r(hGEOSCtxt, pseq, i, x[i]);
		GEOSCoordSeq_setY_r(hGEOSCtxt, pseq, i, y[i]);
	}
	GEOSGeometry* g = GEOSGeom_createLinearRing_r(hGEOSCtxt, pseq);
	return g;
}


inline GEOSGeometry* geos_polygon(SpatPart g, GEOSContextHandle_t hGEOSCtxt) {
	GEOSGeometry* shell = geos_linearRing(g.x, g.y, hGEOSCtxt);

	if (g.hasHoles()) {
		size_t nh=0;
		std::vector<GEOSGeometry*> holes;
		holes.reserve(g.nHoles());
		for (size_t k=0; k < g.nHoles(); k++) {
			SpatHole h = g.getHole(k);
			GEOSGeometry* glr = geos_linearRing(h.x, h.y, hGEOSCtxt);
			if (glr != NULL) {
				holes.push_back(glr);
				nh++;
			}
		}
		GEOSGeometry* geom = GEOSGeom_createPolygon_r(hGEOSCtxt, shell, &holes[0], nh);
		return geom;
	} else {
		GEOSGeometry* geom = GEOSGeom_createPolygon_r(hGEOSCtxt, shell, NULL, 0);
		return geom;
	}
}


inline std::vector<GeomPtr> geos_geoms(SpatVector *v, GEOSContextHandle_t hGEOSCtxt) {
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
			GEOSGeometry* gcol = (geoms.size() == 1) ? geoms[0] :
				GEOSGeom_createCollection_r(hGEOSCtxt, GEOS_MULTILINESTRING, &geoms[0], np);
			g.push_back( geos_ptr(gcol, hGEOSCtxt) );
		}

	} else { // polygons

//		std::vector<std::vector<double>> hx, hy;
		for (size_t i=0; i<n; i++) {
			SpatGeom svg = v->getGeom(i);
			size_t np = svg.size();
			std::vector<GEOSGeometry*> geoms;
			geoms.reserve(np);
			for (size_t j=0; j < np; j++) {
				SpatPart svp = svg.getPart(j);
				GEOSGeometry* gp = geos_polygon(svp, hGEOSCtxt);

				if (gp != NULL) {
					geoms.push_back(gp);
				}
			}
			//Rcpp::Rcout << np << std::endl;
			GEOSGeometry* gcol = (geoms.size() == 1) ? geoms[0] :
				GEOSGeom_createCollection_r(hGEOSCtxt, GEOS_MULTIPOLYGON, &geoms[0], geoms.size());
			g.push_back( geos_ptr(gcol, hGEOSCtxt));
		}
	}
	return g;
}



inline SpatVector vect_from_geos(std::vector<GeomPtr> &geoms , GEOSContextHandle_t hGEOSCtxt, std::string vt) {

	SpatVector out;
	SpatVector v;

	size_t ng = geoms.size();
	std::vector<unsigned> gid, gp, hole;
	std::vector<double> x, y;
	bool xok, yok;

	if ((vt == "points") || (vt == "lines")) {	
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
#ifdef useRcpp
					Rcpp::Rcout << "exception 99" << std::endl;
#else
					std::cout <<  "exception 66" << std::endl;
#endif
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

#ifdef useRcpp
						Rcpp::Rcout << "exception 909" << std::endl;
#else 
						std::cout << "exception 606" << std::endl; 
#endif
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
							hole.push_back(h+1);
						}
					}
				}
			}
		}	
	}	
	v.setGeometry(vt, gid, gp, x, y, hole);
	return v;
}



inline bool pointsFromGeom(GEOSContextHandle_t hGEOSCtxt, const GEOSGeometry* part, 
const unsigned i, const unsigned j, std::vector<double> &x, std::vector<double> &y, 
std::vector<unsigned> &gid, std::vector<unsigned> &gp, std::vector<unsigned> &hole, std::string &msg) {

	const GEOSCoordSequence* crds = GEOSGeom_getCoordSeq_r(hGEOSCtxt, part); 		
	int npts = -1;
	npts = GEOSGetNumCoordinates_r(hGEOSCtxt, part);
	if (npts < 0) {
		msg = "GEOS exception 9";
		return false;
	} 
	if (npts == 0) { // for #813
		x.push_back(NAN);
		y.push_back(NAN);
		gid.push_back(i);			
		gp.push_back(j);			
		hole.push_back(0);
		return true;
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



inline bool polysFromGeom(GEOSContextHandle_t hGEOSCtxt, const GEOSGeometry* part, 
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

	if (npts == 0) { // for #813
		x.push_back(NAN);
		y.push_back(NAN);
		gid.push_back(i);			
		gp.push_back(j);			
		hole.push_back(0);
		return true;
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
				hole.push_back(h+1);
			}
		}
	}
	return true;
}


inline void emptyGeom(const unsigned i, std::vector<double> &x, std::vector<double> &y, 
std::vector<unsigned> &gid, std::vector<unsigned> &gp, std::vector<unsigned> &hole) {
	x.push_back(NAN);
	y.push_back(NAN);
	gid.push_back(i);			
	gp.push_back(0);
	hole.push_back(0);
}


inline SpatVectorCollection coll_from_geos(std::vector<GeomPtr> &geoms, GEOSContextHandle_t hGEOSCtxt, const std::vector<long> &ids = std::vector<long>(), bool keepnull=true, bool increment = true) {

	SpatVectorCollection out;

	std::vector<unsigned> pt_gid, pt_gp, pt_hole;
	std::vector<unsigned> ln_gid, ln_gp, ln_hole;
	std::vector<unsigned> pl_gid, pl_gp, pl_hole;
	std::vector<double> pt_x, pt_y, ln_x, ln_y, pl_x, pl_y;
	std::vector<long> pts_ids, lin_ids, pol_ids;
	
	bool track_ids = !ids.empty(); 	
	if (track_ids) {
		pol_ids.reserve(geoms.size());
	}

	std::string msg;
	size_t ng = geoms.size();
	size_t f = 0;
	for(size_t i = 0; i < ng; i++) {
		const GEOSGeometry* g = geoms[i].get();
		char* geostype = GEOSGeomType_r(hGEOSCtxt, g);
		std::string gt = geostype;
		free(geostype);
		size_t np = GEOSGetNumGeometries_r(hGEOSCtxt, g);

		if (gt == "Point" || gt == "MultiPoint") {
			if (np == 0 && keepnull) {
				emptyGeom(f, pt_x, pt_y, pt_gid, pt_gp, pt_hole);
			}
			for(size_t j = 0; j<np; j++) {
				const GEOSGeometry* part = GEOSGetGeometryN_r(hGEOSCtxt, g, j);
				if (!pointsFromGeom(hGEOSCtxt, part, f, j, pt_x, pt_y, pt_gid, pt_gp, pt_hole, msg)) {
					out.setError(msg);
					return out;
				}
			}
			if (track_ids) pts_ids.push_back(ids[i]);
			f++;
		} else if (gt == "LineString" || gt == "MultiLineString") {
			if (np == 0 && keepnull) {
				emptyGeom(f, ln_x, ln_y, ln_gid, ln_gp, ln_hole);
			}
			for(size_t j = 0; j<np; j++) {
				const GEOSGeometry* part = GEOSGetGeometryN_r(hGEOSCtxt, g, j);
				if (!pointsFromGeom(hGEOSCtxt, part, f, j, ln_x, ln_y, ln_gid, ln_gp, ln_hole, msg)) {
					out.setError(msg);
					return out;
				}
			}
			if (track_ids) lin_ids.push_back(ids[i]);
			f++;
		} else if (gt == "Polygon" || gt == "MultiPolygon") {
			if (np == 0 && keepnull) {
				emptyGeom(f, pl_x, pl_y, pl_gid, pl_gp, pl_hole);
			}
			for(size_t j = 0; j<np; j++) {
				const GEOSGeometry* part = GEOSGetGeometryN_r(hGEOSCtxt, g, j);
				if (!polysFromGeom(hGEOSCtxt, part, f, j, pl_x, pl_y, pl_gid, pl_gp, pl_hole, msg)) {
					out.setError(msg);
					return out;
				}
			}
			if (track_ids) pol_ids.push_back(ids[i]);
			f++;

		} else if (gt == "GeometryCollection") {

			bool first=true;
			size_t kk = 0; // introduced for intersect
			for(size_t j = 0; j<np; j++) {

				const GEOSGeometry* gg = GEOSGetGeometryN_r(hGEOSCtxt, g, j);

				char* geostype = GEOSGeomType_r(hGEOSCtxt, gg);
				std::string ggt = geostype;
				free(geostype);
				size_t npp = GEOSGetNumGeometries_r(hGEOSCtxt, gg);

				//Rcpp::Rcout << geostype << " " << npp << std::endl;

				if (npp == 0 && keepnull) {
					if (ggt == "Polygon" || ggt == "MultiPolygon") {
						emptyGeom(f, pl_x, pl_y, pl_gid, pl_gp, pl_hole);
						if (track_ids) pol_ids.push_back(ids[i]);
					} else if (ggt == "Point" || ggt == "MultiPoint") {
						emptyGeom(f, pt_x, pt_y, pt_gid, pt_gp, pt_hole);
						if (track_ids) pts_ids.push_back(ids[i]);
					} else if (ggt == "LineString" || ggt == "MultiLineString") {
						emptyGeom(f, ln_x, ln_y, ln_gid, ln_gp, ln_hole);
						if (track_ids) lin_ids.push_back(ids[i]);
					}
					if (increment) f++;
				}

				for(size_t k = 0; k<npp; k++) {

					const GEOSGeometry* part = GEOSGetGeometryN_r(hGEOSCtxt, gg, k);

					if (ggt == "Polygon" || ggt == "MultiPolygon") {
						if (!polysFromGeom(hGEOSCtxt, part, f, kk, pl_x, pl_y, pl_gid, pl_gp, pl_hole, msg)) {
							out.setError(msg);
							return out;
						}
						if (track_ids && first) {
							pol_ids.push_back(ids[i]);
							first=false;
						}
					} else if (ggt == "Point" || ggt == "MultiPoint") {
						if (!pointsFromGeom(hGEOSCtxt, part, f, k, pt_x, pt_y, pt_gid, pt_gp, pt_hole, msg)) {
							out.setError(msg);
							return out;
						}
						if (track_ids && first) {
							pts_ids.push_back(ids[i]);
							first=false;
						}

					} else if (ggt == "LineString" || ggt == "MultiLineString") {
						if (!pointsFromGeom(hGEOSCtxt, part, f, k, ln_x, ln_y, ln_gid, ln_gp, ln_hole, msg)) {
							out.setError(msg);
							return out;
						}
						if (track_ids && first) {
							lin_ids.push_back(ids[i]);
							first=false;
						}

					} else {
						out.addWarning("unhandeled Collection geom: " + ggt);
					}
					kk++;
				}
				if (increment) f++;
			}
			if (!increment) f++;
		} else {
			out.setError("what is this: " + gt + "?");
		}
	}

	if (!pl_x.empty()) {
		SpatVector v;
		v.setGeometry("polygons", pl_gid, pl_gp, pl_x, pl_y, pl_hole);
		if (track_ids) v.df.add_column(pol_ids, "ids");
		out.push_back(v);
		//Rcpp::Rcout << "pls" << std::endl;
	}
	if (!ln_x.empty()) {
		SpatVector v;
		v.setGeometry("lines", ln_gid, ln_gp, ln_x, ln_y, ln_hole);
		if (track_ids) v.df.add_column(lin_ids, "ids");
		out.push_back(v);
		//Rcpp::Rcout << "lns" << std::endl;
	}

	if (!pt_x.empty()) {
		SpatVector v;
		v.setGeometry("points", pt_gid, pt_gp, pt_x, pt_y, pt_hole);

		if (track_ids) v.df.add_column(pts_ids, "ids");
		out.push_back(v);
		//Rcpp::Rcout << "pts" << std::endl;
	}
	return out;
}


