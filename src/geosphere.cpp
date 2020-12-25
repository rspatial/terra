// author Robert Hijmans
// October 2009
// version 0.1
// license GPL

/*
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#ifndef PIdiv180
#define PIdiv180 (M_PI / 180)
#endif

#ifndef geod_a
#define geod_a (6378137)
#endif

#ifndef geod_f
#define geod_f (0.0033528106647474805)
#endif


#include <vector>
#include <math.h>
#include "geodesic.h"
#include "distance.h"


void DegtoRad(double &deg) {
	deg *= 0.0174532925199433;
}

void normalizeLonDeg(double &x) {
	x = fmod((x + 180), 360 - 180);
}

void normalizeLonRad(double &x) {
	x = fmod((x + M_PI), (2 * M_PI) - M_PI);
}


bool antipodal(double x1, double y1, double x2, const double &y2, const double &tol=1e-9) {
	normalizeLonDeg(x1);
	normalizeLonDeg(x2);
	double diflon = std::abs(x1 - x2);
	double diflat = std::abs(y1 - y2);
	DegtoRad(y1);
	return ((diflat < tol) && (cos(y1) * abs(fmod(diflon, 360) - 180) < tol));
}


void antipode(const double &x, const double &y, double &xa, double &ya) {
	xa = x + 180;
	normalizeLonDeg(xa);
	ya = -y;
}

double bearing(double lon1, double lat1, double lon2, double lat2 ) {
	double s12, azi1, azi2;
	struct geod_geodesic g;
	geod_init(&g, geod_a, geod_f);
	geod_inverse(&g, lat1, lon1, lat2, lon2, &s12, &azi1, &azi2);
  	return azi1;
}


void destPoint(double lon1, double lat1, double lon2, double lat2, double d, double &dlon, double &dlat) {
	double s12, azi1, azi2;
	struct geod_geodesic g;
	geod_init(&g, geod_a, geod_f);
	geod_inverse(&g, lat1, lon1, lat2, lon2, &s12, &azi1, &azi2);
	geod_direct(&g, lat1, lon1, azi1, d, &dlat, &dlon, &azi2);
}


void intermediate_points(double lon1, double lat1, double lon2, double lat2, size_t n, std::vector<double> &lon, std::vector<double> &lat) {

	lon.resize(0);
	lat.resize(0);

	if ((lon1 == lon2) && (lat1 == lat2)) {
		return;
	}

	if (antipodal(lon1, lat1, lon2, lat2)) {
		lon.resize(1, NAN);
		lat.resize(1, NAN);
		return;
	}

	lon.reserve(n);
	lat.reserve(n);
	n++;

	double d = distance_lonlat(lon1, lat1, lon2, lat2, geod_a, geod_f);
	double step = d / n;

	double dx, dy;
	for (size_t i=1; i<n; i++) {
		destPoint(lon1, lat1, lon2, lat2, step*(i+1), dx, dy);
		lon.push_back(dx);
		lat.push_back(dy);
	}

	//if (addStartEnd) {
	//	x = rbind(p[i,1:2,drop=FALSE], x, p[i,3:4,drop=FALSE])
	//}		
	//if (breakAtDateLine) {
	//	res[[i]] <- .breakAtDateLine(x)
	//}

}

*/
