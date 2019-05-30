// Copyright (c) 2018-2019  Robert J. Hijmans
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



#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#ifndef M_2PI
#define M_2PI (6.28318530717958647692)
#endif

#include <vector>
#include <math.h>
#include "GeographicLib_geodesic.h"
#include "recycle.h"


double distance_lonlat(double lon1, double lat1, double lon2, double lat2, double a, double f) {
	double s12, azi1, azi2;
	struct geod_geodesic g;
	geod_init(&g, a, f);
	geod_inverse(&g, lat1, lon1, lat2, lon2, &s12, &azi1, &azi2);
  	return s12;
}


std::vector<double> distance_lonlat(std::vector<double> &lon1, std::vector<double> &lat1, std::vector<double> &lon2, std::vector<double> &lat2, double a, double f) {
    recycle(lon1, lon2);
    recycle(lat1, lat2);
	std::vector<double> r (lon1.size());
	double azi1, azi2;
	struct geod_geodesic g;
	geod_init(&g, a, f);
	size_t n = lat1.size();
  	for (size_t i=0; i < n; i++) {
		geod_inverse(&g, lat1[i], lon1[i], lat2[i], lon2[i], &r[i], &azi1, &azi2);
	}
  	return r;
}


std::vector<double> distance_lonlat_vd(std::vector<double> &lon1, std::vector<double> &lat1, double lon2, double lat2) {
	double a = 6378137;
	double f = 1 / 298.257223563;
	std::vector<double> vlon2(lon1.size(), lon2);
	std::vector<double> vlat2(lat1.size(), lat2);
    return distance_lonlat(lon1, lat1, vlon2, vlat2, a, f);
}

std::vector<double> distanceToNearest_lonlat(std::vector<double> lon1, std::vector<double> lat1, std::vector<double> lon2, std::vector<double> lat2, double a, double f) {
	double azi1, azi2, s12;
	size_t n = lon1.size();
	size_t m = lon2.size();
	std::vector<double> r(n);

	struct geod_geodesic g;
	geod_init(&g, a, f);
  	for (size_t i=0; i < n; i++) {
		geod_inverse(&g, lat1[i], lon1[i], lat2[0], lon2[0], &r[i], &azi1, &azi2);
		for (size_t j=1; j<m; j++) {
			geod_inverse(&g, lat1[i], lon1[i], lat2[j], lon2[j], &s12, &azi1, &azi2);
			if (s12 < r[i]) {
				r[i] = s12;
			}
		}
	}
  	return r;
}


double distance_plane(double x1, double y1, double x2, double y2) {
	return( sqrt(pow((x2-x1),2) + pow((y2-y1), 2)) );
}


std::vector<double> distance_plane(std::vector<double> &x1, std::vector<double> &y1, std::vector<double> &x2, std::vector<double> &y2) {
    recycle(x1, x2);
    recycle(y1, y2);
	std::vector<double> r (x1.size());
	size_t n = x1.size();
  	for (size_t i=0; i < n; i++) {
		r[i] = sqrt(pow((x2[i]-x1[i]),2) + pow((y2[i]-y1[i]), 2));
	}
  	return r;
}

std::vector<double> distance_plane_vd(std::vector<double> &x1, std::vector<double> &y1, double x2, double y2) {
	std::vector<double> vx2(x1.size(), x2);
	std::vector<double> vy2(y1.size(), y2);
  	return distance_plane(x1, y1, vx2, vy2);
}

std::vector<double> distanceToNearest_plane(std::vector<double> x1, std::vector<double> y1, std::vector<double> x2, std::vector<double> y2) {
	size_t n = x1.size();
	size_t m = x2.size();
	std::vector<double> r(n);
	double d;
  	for (size_t i=0; i < n; i++) {
		r[i] = sqrt(pow((x2[0]-x1[i]),2) + pow((y2[0]-y1[i]), 2));
		for (size_t j=1; j < m; j++) {
			d = sqrt(pow((x2[j]-x1[i]),2) + pow((y2[j]-y1[i]), 2));
			if (d < r[i]) {
				r[i] = d;
			}
		}
	}
  	return r;
}



// Convert degrees to radians
double toRad(double deg) {
	return( deg * 0.0174532925199433 );
}

double toDeg(double rad) {
	return( rad * 57.2957795130823 );
}



double direction_lonlat(double lon1, double lat1, double lon2, double lat2, bool degrees, double a, double f) {
	double s12, azi1, azi2;
	struct geod_geodesic g;
	geod_init(&g, a, f);
	geod_inverse(&g, lat1, lon1, lat2, lon2, &s12, &azi1, &azi2);
	if (!degrees) {
		return(toRad(azi1));
	}
	return( azi1) ;
}

std::vector<double> direction_lonlat(std::vector<double> lon1, std::vector<double> lat1, std::vector<double> lon2, std::vector<double> lat2, bool degrees, double a, double f) {
// lonlat1 and lonlat2 should have the same length

	std::vector<double> azi1(lon1.size());
	double s12, azi2;
	struct geod_geodesic g;
	geod_init(&g, a, f);
	size_t n = lat1.size();
	if (degrees) {
		for (size_t i=0; i < n; i++) {
			geod_inverse(&g, lat1[i], lon1[i], lat2[i], lon2[i], &s12, &azi1[i], &azi2);
		}
	} else {
		for (size_t i=0; i < n; i++) {
			geod_inverse(&g, lat1[i], lon1[i], lat2[i], lon2[i], &s12, &azi1[i], &azi2);
			azi1[i] = toRad(azi1[i]);
		}
	}
  	return azi1;
}


std::vector<double> directionToNearest_lonlat(std::vector<double> lon1, std::vector<double> lat1, std::vector<double> lon2, std::vector<double> lat2, bool degrees, bool from, double a, double f) {
	double azi1, azi2, s12, dist;
	size_t n = lon1.size();
	size_t m = lon2.size();
	std::vector<double> azi(n);

	struct geod_geodesic g;
	geod_init(&g, a, f);

	if (from) {
		for (size_t i=0; i < n; i++) {
			geod_inverse(&g, lat2[0], lon2[0], lat1[i], lon1[i], &dist, &azi1, &azi2);
			azi[i] = azi1;
			for (size_t j=1; j<m; j++) {
				geod_inverse(&g, lat2[j], lon2[j], lat1[i], lon1[i], &s12, &azi1, &azi2);
				if (s12 < dist) {
					azi[i] = azi1;
				}
			}
			if (!degrees) {
				azi[i] = toRad(azi[i]);
			}
		}

	} else {
		for (size_t i=0; i < n; i++) {
			geod_inverse(&g, lat1[i], lon1[i], lat2[0], lon2[0], &dist, &azi1, &azi2);
			azi[i] = azi1;
			for (size_t j=1; j<m; j++) {
				geod_inverse(&g, lat1[i], lon1[i], lat2[j], lon2[j], &s12, &azi1, &azi2);
				if (s12 < dist) {
					azi[i] = azi1;
				}
			}
			if (!degrees) {
				azi[i] = toRad(azi[i]);
			}
		}
	}
  	return azi;
}



double direction_plane(double x1, double y1, double x2, double y2, bool degrees) {
	double a;
	a = fmod(atan2( x2 - x1, y2 - y1), M_2PI);
	a = (a < 0 ? a + M_2PI : a );
	return (degrees ? toDeg(a) : a);
}



std::vector<double> direction_plane(std::vector<double> x1, std::vector<double> y1, std::vector<double> x2, std::vector<double> y2, bool degrees) {
// xy1 and xy2 should have the same length
	std::vector<double> r (x1.size());
	//double a;
	size_t n = x1.size();
  	for (size_t i=0; i < n; i++) {
		r[i] = direction_plane(x1[i], y1[i], x2[i], y2[i], degrees);
	}
  	return r;
}



std::vector<double> directionToNearest_plane(std::vector<double> x1, std::vector<double> y1, std::vector<double> x2, std::vector<double> y2, bool degrees, bool from) {
	size_t n = x1.size();
	size_t m = x2.size();
	std::vector<double> r(n);
	double d, mind;
	size_t minj;
	if (from) {
		for (size_t i = 0; i < n; i++) {
			mind = distance_plane(x1[i], y1[i], x2[0], y2[0]);
			minj = 0;
			for (size_t j = 1; j < m; j++) {
				d = distance_plane(x1[i], y1[i], x2[j], y2[j]);
				if (d < mind) {
					mind = d;
					minj = j;
				}
			}
			r[i] = direction_plane(x2[minj], y2[minj], x1[i], y1[i], degrees);
		}
	} else {
		for (size_t i = 0; i < n; i++) {
			mind = distance_plane(x1[i], y1[i], x2[0], y2[0]);
			minj = 0;
			for (size_t j = 1; j < m; j++) {
				d = distance_plane(x1[i], y1[i], x2[j], y2[j]);
				if (d < mind) {
					mind = d;
					minj = j;
				}
			}
			r[i] = direction_plane(x1[i], y1[i], x2[minj], y2[minj], degrees);
		}
	}
  	return r;
}





std::vector<double> destpoint_lonlat(double longitude, double latitude, double  bearing, double distance, double a, double f) {
	struct geod_geodesic g;
	geod_init(&g, a, f);
	double lat2, lon2, azi2;
	geod_direct(&g, latitude, longitude, bearing, distance, &lat2, &lon2, &azi2);
	std::vector<double> out = {lon2, lat2, azi2 };
	return out;
}


std::vector<std::vector<double> > destpoint_lonlat(std::vector<double> longitude, std::vector<double> latitude, std::vector<double> bearing, std::vector<double> distance, double a, double f) {
	struct geod_geodesic g;
	geod_init(&g, a, f);
	size_t n = longitude.size();
	std::vector<std::vector<double> > out;
	out.reserve(n);
	double lat2, lon2, azi2;
	for (size_t i=0; i < n; i++) {
		geod_direct(&g, latitude[i], longitude[i], bearing[i], distance[i], &lat2, &lon2, &azi2);
		out.push_back( {lon2, lat2, azi2 });
	}
	return out;
}


std::vector<double> destpoint_plane(double x, double y, double bearing, double distance) {
	bearing = bearing * M_PI / 180;
	x += distance * cos(bearing);
	y += distance * sin(bearing);
	std::vector<double> out = {x, y};
	return(out);
}


std::vector<std::vector<double> > destpoint_plane(std::vector<double>  x, std::vector<double>  y, std::vector<double>  bearing, std::vector<double>  distance) {
	size_t n = x.size();
	std::vector<std::vector<double> > out;
	out.reserve(n);
	double xd, yd, b;
	for (size_t i=0; i < n; i++) {
		b = bearing[i] * M_PI / 180;
		xd = x[i] + distance[i] * cos(b);
		yd = y[i] + distance[i] * sin(b);
		out.push_back( {xd, yd });
	}
	return(out);
}

