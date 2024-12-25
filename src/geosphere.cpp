// Copyright (c) 2018-2025  Robert J. Hijmans
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
#define WGS84_a 6378137.0
#endif

#ifndef WGS84_f
#define WGS84_f 1/298.257223563
#endif


#include "Rcpp.h"

inline void normLon(double &lon) {
	lon = fmod(lon + 180, 360.) - 180;
}

inline void normLonRad(double &lon) {
	lon = fmod(lon + M_PI, M_2PI) - M_PI;
}


inline void deg2rad(std::vector<double> &x) {
	const double f = 0.0174532925199433;
	for (double& d : x) d *= f;
}

inline void deg2rad(double &x) {
	const double f = 0.0174532925199433;
	x *= f;
}


inline double get_sign(const double &x) {
	return (x > 0.0) ? 1.0 : (x < 0.0) ? -1.0 : 0;
}


// [[Rcpp::export]]
double dist_lonlat(const double &lon1, const double &lat1, const double &lon2, const double &lat2) {
	//double a = 6378137.0;
	//double f = 1/298.257223563;
	double s12, azi1, azi2;
	struct geod_geodesic g;
	geod_init(&g, WGS84_a, WGS84_f);
	geod_inverse(&g, lat1, lon1, lat2, lon2, &s12, &azi1, &azi2);
  	return s12;
}


// [[Rcpp::export]]
double dist_cosine(double lon1, double lat1, double lon2, double lat2, const double &r = 6378137) {
	deg2rad(lon1);
	deg2rad(lon2);
	deg2rad(lat1);
	deg2rad(lat2);
	return r * acos((sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1-lon2)));
}

// [[Rcpp::export]]
double dist_cosine_rad(const double &lon1, const double &lat1, const double &lon2, const double &lat2, const double &r = 6378137) {
	return r * acos((sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1-lon2)));
}


// [[Rcpp::export]]
void dest_lonlat(double slon, double slat, double sazi, double dist, double &dlon, double &dlat, double &dazi) {
	//double a = 6378137.0;
	//double f = 1/298.257223563;
	struct geod_geodesic g;
	geod_init(&g, WGS84_a, WGS84_f);
	geod_direct(&g, slat, slon, sazi, dist, &dlat, &dlon, &dazi);
}


// [[Rcpp::export]]
double dir_lonlat(double lon1, double lat1, double lon2, double lat2) {
	//double a = 6378137.0;
	//double f = 1/298.257223563;

	double s12, azi1, azi2;
	struct geod_geodesic g;
	geod_init(&g, WGS84_a, WGS84_f);
	geod_inverse(&g, lat1, lon1, lat2, lon2, &s12, &azi1, &azi2);
	return( azi1) ;
}


// [[Rcpp::export]]
double dir_rad(double lon1, double lat1, double lon2, double lat2) {

	if ((lon1 == lon2) && (lat1 == lat2)) return 0; // NAN?
	
	double dLon = lon2 - lon1;
    double y = sin(dLon)  * cos(lat2); 
    double x = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dLon); 

    double azm = atan2(y, x);
	azm = fmod(azm+M_PI, M_PI);
	return azm > M_PI ? -(M_PI - azm) : azm;
}



// [[Rcpp::export]]
double dist2track(double lon1, double lat1, double lon2, double lat2, double plon, double plat, bool sign, double r=6378137) {
	double a = 1;
	double f = 0;
	struct geod_geodesic geod;
	geod_init(&geod, a, f);
	

	double d, b2, b3, azi;
	geod_inverse(&geod, lat1, lon1, lat2, lon2, &d, &b2, &azi);
	geod_inverse(&geod, lat1, lon1, plat, plon, &d, &b3, &azi);
	double toRad = M_PI / 180.;
	b2 *= toRad;
	b3 *= toRad;
	double xtr = asin(sin(b3-b2) * sin(d)) * r;
	return sign ? xtr : fabs(xtr);
}


// [[Rcpp::export]]
double dist2track_cosine_rad(double lon1, double lat1, double lon2, double lat2, double plon, double plat, bool sign, double r=6378137) {
	// same results as above but probably much faster
	double b2 = dir_rad(lon1, lat1, lon2, lat2);
	double b3 = dir_rad(lon1, lat1, plon, plat);
	double d = dist_cosine_rad(lon1, lat1, plon, plat, 1);
	double xtr = asin(sin(b3-b2) * sin(d)) * r;
	return sign ? xtr : fabs(xtr);
}




// [[Rcpp::export]]
double alongTrackDistance(double lon1, double lat1, double lon2, double lat2, double plon, double plat, double r=6378137) {
	double a = 1;
	double f = 0;
	struct geod_geodesic geod;
	geod_init(&geod, a, f);
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



// [[Rcpp::export]]
double alongTrackDistance_rad(double lon1, double lat1, double lon2, double lat2, double plon, double plat, double r=6378137) {

	double tc = dir_rad(lon1, lat1, lon2, lat2); // * toRad
	double tcp = dir_rad(lon1, lat1, plon, plat); // * toRad
    double dp = dist_cosine_rad(lon1, lat1, plon, plat, 1);
	double xtr = asin(sin(tcp-tc) * sin(dp));

// +1/-1 for ahead/behind [lat1,lon1]
	double bearing = get_sign(cos(tc - tcp));
	double angle = cos(dp) / cos(xtr);

// Fixing limits for the angle between [-1, 1] to avoid NaNs from acos
	angle = angle > 1 ? 1 : angle < -1 ? -1 : angle;
	double dist = bearing * acos(angle) * r;
	return fabs(dist);
}



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
double dist2segment_cosine_rad(double plon, double plat, double lon1, double lat1, double lon2, double lat2, double r=6378137) {
// the alongTrackDistance is the length of the path along the great circle to the point of intersection
// there are two, depending on which node you start
// we want to use the min, but the max needs to be < segment length
	double seglength = dist_cosine_rad(lon1, lat1, lon2, lat2, r);
	double trackdist1 = alongTrackDistance_rad(lon1, lat1, lon2, lat2, plon, plat, r);
	double trackdist2 = alongTrackDistance_rad(lon2, lat2, lon1, lat1, plon, plat, r);
	if ((trackdist1 >= seglength) || (trackdist2 >= seglength)) {
		double d1 = dist_cosine_rad(lon1, lat1, plon, plat, r);
		double d2 = dist_cosine_rad(lon2, lat2, plon, plat, r);
		return d1 < d2 ? d1 : d2;
	}
	return dist2track_cosine_rad(lon1, lat1, lon2, lat2, plon, plat, false, r);
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


std::vector<double> SpatVector::linedistLonLat(SpatVector x, std::string unit) {

	std::vector<double> d;

// x is points
	if ((x.type() != "points") || (type() == "points")) {
		setError("wrong place!");
		return d;
	}

	double r = 6378137;
	if (unit == "km") {
		r = 6378.137;
	} else if (unit == "ha") {
		r = 63781.37;	
	}

	std::vector<std::vector<double>> pxy = x.coordinates();
	deg2rad(pxy[0]);
	deg2rad(pxy[1]);
	size_t np = pxy[0].size();
	size_t ng = size();

//	dd.reserve(np*ng);
	double inf = std::numeric_limits<double>::infinity();
	d.resize(np, inf);

	std::vector<double> px, py;

	if (type() == "polygons") {
		std::vector<int> insect = relate(x, "intersects", true, true);
		for (size_t g=0; g<ng; g++) {
			size_t nparts = geoms[g].size();
			for (size_t h=0; h<nparts; h++) {
				px = geoms[g].parts[h].x;
				py = geoms[g].parts[h].y;
				deg2rad(px);
				deg2rad(py);
				size_t nseg = px.size() - 1;
				for (size_t i=0; i<np; i++) {
					if (d[i] != 0) {
						if (insect[i] == 1) {
							d[i] = 0;
						}
						for (size_t j=0; j<nseg; j++) {
							d[i] = std::min(d[i],
								dist2segment_cosine_rad(pxy[0][i], pxy[1][i], px[j], py[j], px[j+1], py[j+1], r));
						}
					}
				}
				size_t nh = geoms[g].parts[h].nHoles();
				for (size_t k=0; k < nh; k++) {
					px = geoms[g].parts[h].holes[k].x;
					py = geoms[g].parts[h].holes[k].y;
					deg2rad(px);
					deg2rad(py);
					size_t nseg = px.size() - 1;
					for (size_t i=0; i<np; i++) {
						if ((d[i] != 0) && (insect[i] == 0)) {
							for (size_t j=0; j<nseg; j++) {
								d[i] = std::min(d[i],
									dist2segment_cosine_rad(pxy[0][i], pxy[1][i], px[j], py[j], px[j+1], py[j+1], r));
							}
						}
					}
				}
			}
		}
	} else {
		for (size_t g=0; g<ng; g++) {
			size_t nparts = geoms[g].size();
			for (size_t h=0; h<nparts; h++) {
				px = geoms[g].parts[h].x;
				py = geoms[g].parts[h].y;
				deg2rad(px);
				deg2rad(py);
				size_t nseg = px.size() - 1;
				for (size_t i=0; i<np; i++) {
					for (size_t j=0; j<nseg; j++) {
						d[i] = std::min(d[i],
							dist2segment_cosine_rad(pxy[0][i], pxy[1][i], px[j], py[j], px[j+1], py[j+1], r));
					}
				}
			}
		}
	}
	return d;
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


void make_dense_lonlat(std::vector<double> &lon, std::vector<double> &lat, const double &interval, const bool &adjust, geod_geodesic &g) {
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
		//double hlat = lat[i] + (lat[i+1] - lat[i])/2;
		//double hlon = lon[i] + (lon[i+1] - lon[i])/2;
		//geod_inverse(&g, lat[i], lon[i], hlat, hlon, &d1, &azi1, &azi2);
		//geod_inverse(&g, hlat, hlon, lat[i+1], lon[i+1], &d2, &azi1, &azi2);
		//double d = d1 + d2;
		geod_inverse(&g, lat[i], lon[i], lat[i+1], lon[i+1], &d, &azi1, &azi2);
		size_t n = floor(d / interval);
		xout.push_back(lon[i]);
		yout.push_back(lat[i]);
		if (n < 2) {
			continue;
		}
		double step = adjust ? d / n : interval;
		double newlat, newlon;
		for (size_t j=1; j<n; j++) {
			geod_direct(&g, lat[i], lon[i], azi1, step*j, &newlat, &newlon, &azi2);
			// avoid -180 to 180 jumps
			if ((lon[i] == -180) && (newlon == 180)) newlon = -180;
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


SpatVector SpatVector::densify(double interval, bool adjust, bool ignorelonlat) {

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
	size_t n = size();
	out.reserve(n);
	if (is_lonlat() && (!ignorelonlat)) {
		double a = 6378137.0;
		double f = 1/298.257223563;
		struct geod_geodesic geod;
		geod_init(&geod, a, f);

		for (size_t i=0; i<n; i++) {
			SpatGeom g = geoms[i];
			for (size_t j=0; j < geoms[i].size(); j++) {
				make_dense_lonlat(g.parts[j].x, g.parts[j].y, interval, adjust, geod);
				if (g.parts[j].hasHoles()) {
					for (size_t k=0; k < g.parts[j].nHoles(); k++) {
						make_dense_lonlat(g.parts[j].holes[k].x, g.parts[j].holes[k].y, interval, adjust, geod);
					}
				}
			}
			g.computeExtent();
			out.addGeom(g);
		}
	} else {

		for (size_t i=0; i<n; i++) {
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
	out.df = df;
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

