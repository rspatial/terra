// Copyright (c) 2018-2021  Robert J. Hijmans
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

#include "spatVector.h"
#include "geodesic.h"
#include "recycle.h"

#ifndef M_PI
#define M_PI  (3.1415926535897932384626433)
#endif

#ifndef M_2PI
#define M_2PI (3.1415926535897932384626433 * 2.0)
#endif

#ifndef M_PI_2 
#define M_PI_2 (3.1415926535897932384626433 / 2)
#endif
 
#ifndef WGS84_a 
#define WGS84_a 6378137.0;
#endif

#ifndef WGS84_f 
#define WGS84_f 1/298.257223563;;
#endif

 
inline void normLon(double &lon) {
	lon = fmod(lon + 180, 360.) - 180;
}
 
inline void normLonRad(double &lon) {
	lon = fmod(lon + M_PI, M_2PI) - M_PI;
}

inline double get_sign(double x) {
	if (x > 0.0) return 1.0;
	if (x < 0.0) return -1.0;
	return x;
}


// [[Rcpp::export]]
double dist_lonlat(const double &lon1, const double &lat1, const double &lon2, const double &lat2) {
	double a = 6378137.0;
	double f = 1/298.257223563;
	double s12, azi1, azi2;
	struct geod_geodesic g;
	geod_init(&g, a, f);
	geod_inverse(&g, lat1, lon1, lat2, lon2, &s12, &azi1, &azi2);
  	return s12;
}


void dest_lonlat(double slon, double slat, double sazi, double dist, double &dlon, double &dlat, double &dazi) {
	double a = 6378137.0;
	double f = 1/298.257223563;
	struct geod_geodesic g;
	geod_init(&g, a, f);
	geod_direct(&g, slat, slon, sazi, dist, &dlat, &dlon, &dazi);
}


// [[Rcpp::export]]
double dir_lonlat(double lon1, double lat1, double lon2, double lat2) {
	double a = 6378137.0;
	double f = 1/298.257223563;

	double s12, azi1, azi2;
	struct geod_geodesic g;
	geod_init(&g, a, f);
	geod_inverse(&g, lat1, lon1, lat2, lon2, &s12, &azi1, &azi2);
	return( azi1) ;
}


// [[Rcpp::export]]
double dist2track(double lon1, double lat1, double lon2, double lat2, double plon, double plat, bool sign) {
	double a = 1;
	double f = 0;
	struct geod_geodesic geod;
	geod_init(&geod, a, f);
	double r = 6378137.0;
	
	double d, b2, b3, azi;
	geod_inverse(&geod, lat1, lon1, lat2, lon2, &d, &b2, &azi);
	geod_inverse(&geod, lat1, lon1, plat, plon, &d, &b3, &azi);
	double toRad = M_PI / 180.;
	b2 *= toRad;
	b3 *= toRad;
	double xtr = asin(sin(b3-b2) * sin(d)) * r;
	if (!sign) xtr = fabs(xtr);
	return xtr;
}


// [[Rcpp::export]]
double alongTrackDistance(double lon1, double lat1, double lon2, double lat2, double plon, double plat) {
	double a = 1;
	double f = 0;
	struct geod_geodesic geod;
	geod_init(&geod, a, f);
	double r = 6378137.0;
	
	double d, b2, b3, azi;
	geod_inverse(&geod, lat1, lon1, lat2, lon2, &d, &b2, &azi);
	geod_inverse(&geod, lat1, lon1, plat, plon, &d, &b3, &azi);
	double toRad = M_PI / 180.;
	b2 *= toRad;
	b3 *= toRad;
	double xtr = asin(sin(b3-b2) * sin(d));

	double bsign = get_sign(cos(b2-b3));  
	return fabs(bsign * acos(cos(d) / cos(xtr)) * r);
}



#include "Rcpp.h"

// [[Rcpp::export]]
double dist2segment(double plon, double plat, double lon1, double lat1, double lon2, double lat2) {
			
// the alongTrackDistance is the length of the path along the great circle to the point of intersection
// there are two, depending on which node you start
// we want to use the min, but the max needs to be < segment length
	double seglength = dist_lonlat(lon1, lat1, lon2, lat2);
	double trackdist1 = alongTrackDistance(lon1, lat1, lon2, lat2, plon, plat);
	double trackdist2 = alongTrackDistance(lon2, lat2, lon1, lat1, plon, plat);
	if ((trackdist1 >= seglength) || (trackdist2 >= seglength)) {
		double d1 = dist_lonlat(lon1, lat1, plon, plat);
		double d2 = dist_lonlat(lon2, lat2, plon, plat);
		return d1 < d2 ? d1 : d2;
	}
	return dist2track(lon1, lat1, lon2, lat2, plon, plat, false);
}


// [[Rcpp::export]]
double dist2segmentPoint(double plon, double plat, double lon1, double lat1, double lon2, double lat2, double &ilon, double &ilat) {
			
	double seglength = dist_lonlat(lon1, lat1, lon2, lat2);
	double trackdist1 = alongTrackDistance(lon1, lat1, lon2, lat2, plon, plat);
	double trackdist2 = alongTrackDistance(lon2, lat2, lon1, lat1, plon, plat);
	if ((trackdist1 >= seglength) || (trackdist2 >= seglength)) {
		double d1 = dist_lonlat(lon1, lat1, plon, plat);
		double d2 = dist_lonlat(lat2, lat2, plon, plat);
		if (d1 < d2) {
			ilon = lon1;
			ilat = lat1;			
			return d1;	
		} else {
			ilon = lon2;
			ilat = lat2;			
			return d2;	
		}
	}
	double azi;
	double crossd = dist2track(lon1, lat1, lon2, lat2, plon, plat, false);
	if (trackdist1 < trackdist2) {
		double bear = dir_lonlat(lon1, lat1, lon2, lat2);
		dest_lonlat(lon1, lat1, bear, trackdist1, ilon, ilat, azi);
	} else {
		double bear = dir_lonlat(lon2, lat2, lon1, lat1);
		dest_lonlat(lon2, lat2, bear, trackdist2, ilon, ilat, azi);
	}
	return(crossd);
}


std::vector<double> SpatVector::linedistLonLat(SpatVector x) {
	
	std::vector<std::vector<double>> pxy = x.coordinates();
	size_t np = pxy[0].size();
	size_t ng = size();

	std::vector<double> d, dd;
	dd.reserve(np*ng);
	d.resize(np);
	
	bool poly = type() == "polygons";
	if (poly) {
		SpatVector pg;
		pg.srs = srs;
		std::vector<int> insect;
		for (size_t g=0; g<ng; g++) {
			pg.geoms = { geoms[g] };
			insect = pg.relate(x, "intersects");
			std::vector<std::vector<double>> xy = geoms[g].coordinates();
			size_t nseg = xy[0].size() - 1;
			for (size_t i=0; i<np; i++) {
				if (insect[i]) {
					d[i] = 0;
				} else {
					d[i] = dist2segment(pxy[0][i], pxy[1][i], xy[0][0], xy[1][0], xy[0][1], xy[1][1]);
					for (size_t j=1; j<nseg; j++) {
						d[i] = std::min(d[i], 
						dist2segment(pxy[0][i], pxy[1][i], xy[0][j], xy[1][j], xy[0][j+1], xy[1][j+1]));
					}
				}
			}
		}
		dd.insert(dd.end(), d.begin(), d.end());
	} else {
		for (size_t g=0; g<ng; g++) {
			std::vector<std::vector<double>> xy = geoms[g].coordinates();
			size_t nseg = xy[0].size() - 1;
			for (size_t i=0; i<np; i++) {
				d[i] = dist2segment(pxy[0][i], pxy[1][i], xy[0][0], xy[1][0], xy[0][1], xy[1][1]);
				for (size_t j=1; j<nseg; j++) {
					d[i] = std::min(d[i], 
					dist2segment(pxy[0][i], pxy[1][i], xy[0][j], xy[1][j], xy[0][j+1], xy[1][j+1]));
				}
			}
		}
		dd.insert(dd.end(), d.begin(), d.end());
	}
	return dd;
}
 


// [[Rcpp::export(name = "intermediate")]]
std::vector<std::vector<double>> intermediate(double lon1, double lat1, double lon2, double lat2, int n, double distance) {
	double a = 6378137.0;
	double f = 1/298.257223563;
	struct geod_geodesic geod;
	geod_init(&geod, a, f);
	double d, azi1, azi2;

	std::vector<std::vector<double>> out(2);	
	if (n <= 0) {
		if (distance <= 0) {
			out[0] = {lon1, lon2};
			out[1] = {lon1, lon2};
			return out;
		} else {
			geod_inverse(&geod, lat1, lon1, lat2, lon2, &d, &azi1, &azi2);
			n = std::round(d / distance);
			if (n < 2) {
				out[0] = {lon1, lon2};
				out[1] = {lon1, lon2};
				return out;				
			}
			distance = d / n;
		}
	} else if (n == 1) {
		out[0] = {lon1, lon2};
		out[1] = {lon1, lon2};					
		return out;
	} else {	
		geod_inverse(&geod, lat1, lon1, lat2, lon2, &d, &azi1, &azi2);
		//distance = d / n;
	}
	out[0].resize(n+1);
	out[1].resize(n+1);
	out[0][0] = lon1;
	out[1][0] = lat1;
	for (int i=1; i<n; i++) {
		geod_direct(&geod, lat1, lon1, azi1, distance*i, &out[1][i], &out[0][i], &azi2);
	}
	out[0][n] = lon2;
	out[1][n] = lat2;
	return out;
}


void make_dense_lonlat(std::vector<double> &lon, std::vector<double> &lat, double &interval, bool &adjust, geod_geodesic &g) {
	size_t np = lon.size();
	if (np < 2) {
		return;
	}
	size_t sz = lon.size() * 5;
	std::vector<double> xout, yout;
	xout.reserve(sz);
	yout.reserve(sz);
	for (size_t i=0; i<(np-1); i++) {
		if (xout.size() > sz) {
			sz += (np-i) * 10;
			xout.reserve(sz);
			yout.reserve(sz);			
		}
		double d, azi1, azi2;
		geod_inverse(&g, lat[i], lon[i], lat[i+1], lon[i+1], &d, &azi1, &azi2);
		size_t n = floor(d / interval);
		xout.push_back(lon[i]);
		yout.push_back(lat[i]);
		if (n < 2) {
			continue;
		}
		double step = adjust ? d / n : d;
		double newlat, newlon;
		for (size_t j=1; j<n; j++) {
			geod_direct(&g, lat[i], lon[i], azi1, step*j, &newlat, &newlon, &azi2);
			xout.push_back(newlon);
			yout.push_back(newlat);			
		}
	}
	xout.push_back(lon[np-1]);
	yout.push_back(lat[np-1]);
	lon = std::move(xout);
	lat = std::move(yout);
}

void make_dense_planar(std::vector<double> &x, std::vector<double> &y, double &interval, bool &adjust) {
	size_t np = x.size();
	if (np < 2) {
		return;
	}
	size_t sz = x.size() * 5;
	std::vector<double> xout, yout;
	xout.reserve(sz);
	yout.reserve(sz);

	double pi2 = M_PI * 2;

	for (size_t i=0; i<(np-1); i++) {
		if (xout.size() > sz) {
			sz += (np-i) * 10;
			xout.reserve(sz);
			yout.reserve(sz);			
		}
		double d = sqrt(pow((x[i+1] - x[i]),2) + pow((y[i+1] - y[i]), 2));
		size_t n = floor(d / interval);
		xout.push_back(x[i]);
		yout.push_back(y[i]);
		if (n < 2) {
			continue;
		}
			
		double a = fmod(atan2(x[i+1]-x[i], y[i+1]-y[i]), pi2);
		double step = adjust ? d / n : interval;
		double distx = step * sin(a);
		double disty = step * cos(a);
		for (size_t j=1; j<n; j++) {
			xout.push_back(x[i] + distx * j);
			yout.push_back(y[i] + disty * j);
		}
	}
	xout.push_back(x[np-1]);
	yout.push_back(y[np-1]);
	x = std::move(xout);
	y = std::move(yout);
}


SpatVector SpatVector::densify(double interval, bool adjust) {

	SpatVector out;
	if (type() == "points") {
		out.setError("cannot densify points");
		return out;
	}
	if (interval <= 0) {
		out.setError("the interval must be > 0");
		return out;
	}

	out.srs = srs;
	if (srs.is_empty()) {
		out.setError("crs not defined");
		return(out);
	}
	if (is_lonlat()) {
		double a = 6378137.0;
		double f = 1/298.257223563;
		struct geod_geodesic geod;
		geod_init(&geod, a, f);

		for (size_t i=0; i<size(); i++) {
			SpatGeom g = geoms[i];
			for (size_t j=0; j < geoms[i].size(); j++) {
				make_dense_lonlat(g.parts[j].x, g.parts[j].y, interval, adjust, geod);	
				if (g.parts[j].hasHoles()) {
					for (size_t k=0; k < g.parts[j].nHoles(); k++) {
						make_dense_lonlat(g.parts[j].holes[k].x, g.parts[j].holes[k].y, interval, adjust, geod);
					}
				}
			}
			out.addGeom(g);
		}
	} else {
		
		for (size_t i=0; i<size(); i++) {
			SpatGeom g = geoms[i];
			for (size_t j=0; j < geoms[i].size(); j++) {
				make_dense_planar(g.parts[j].x, g.parts[j].y, interval, adjust);	
				if (g.parts[j].hasHoles()) {
					for (size_t k=0; k < g.parts[j].nHoles(); k++) {
						make_dense_planar(g.parts[j].holes[k].x, g.parts[j].holes[k].y, interval, adjust);
					}
				}
			}
			out.addGeom(g);
		}
	}
	return out;
}
 
 

 
std::vector<bool> antipodal(std::vector<double> lon1, std::vector<double> lat1, std::vector<double> lon2, std::vector<double> lat2, double tol=1e-9) {
	recycle(lon1, lon2);
	recycle(lat1, lat2);
	std::vector<bool> out;
	out.reserve(lon1.size());
	double Pi180 = M_PI / 180.;
	for (size_t i=0; i<lon1.size(); i++){ 
		normLon(lon1[i]);
		normLon(lon2[i]);		
		double diflon = fabs(lon1[i] - lon2[i]);
		double diflat = fabs(lat1[i] + lat2[i]);
		out.push_back(
			(diflat < tol) && ((cos(lat2[i] * Pi180) * fabs(fmod(diflon, 360.) - 180)) < tol)
		);
	}
	return out;
}


void antipodes(std::vector<double> &lon, std::vector<double> &lat) {
	size_t n=lon.size();
	for (size_t i=0; i<n; i++) { 
		lon[i] = lon[i] + 180;
		normLon(lon[i]);
		lat[i] = -lat[i];
	}
}

