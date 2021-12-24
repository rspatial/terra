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
 
 