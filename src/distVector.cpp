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

#include <functional>
#include "spatVector.h"
#include "distance.h"
#include "geosphere.h"
#include "geodesic.h"
#include "vecmath.h"


std::vector<double> SpatVector::nearestDistLonLat(std::vector<double> x, std::vector<double> y, std::string unit, std::string method) {

// for use with rasterize 
	std::vector<double> d;
	double r = 6378137;
	double m = 1;
	if (unit == "km") {
		r = 6378.137;
		m = 0.001;
	}

	std::vector<int> inside;
	if (type() == "polygons") {
//		std::vector<int> insect = relate(x, "intersects", true, true);
		inside = pointInPolygon(x, y);
	}

	std::function<double(double,double,double,double,double,double,double)> d2seg;
	if (method != "geo") {
		deg2rad(x);
		deg2rad(y);
		d2seg = dist2segment_cos;
	} else {
		d2seg = dist2segment_geo;		
	}
	size_t np = x.size();
	size_t ng = size();

	double inf = std::numeric_limits<double>::infinity();
	d.resize(np, inf);

	std::vector<double> vx, vy;

	if (type() == "polygons") {
		for (size_t g=0; g<ng; g++) {
			size_t nparts = geoms[g].size();
			for (size_t h=0; h<nparts; h++) {
				vx = geoms[g].parts[h].x;
				vy = geoms[g].parts[h].y;
				if (method != "geo") {
					deg2rad(vx);
					deg2rad(vy);
				}
				size_t nseg = vx.size() - 1;
				for (size_t i=0; i<np; i++) {
					if (d[i] != 0) {
						if (inside[i] == 1) {
							d[i] = 0;
						} else {
							for (size_t j=0; j<nseg; j++) {
								d[i] = std::min(d[i],
									d2seg(x[i], y[i], vx[j], vy[j], vx[j+1], vy[j+1], r));
							}
						}
					}
				}
				size_t nh = geoms[g].parts[h].nHoles();
				for (size_t k=0; k < nh; k++) {
					vx = geoms[g].parts[h].holes[k].x;
					vy = geoms[g].parts[h].holes[k].y;
					if (method != "geo") {
						deg2rad(vx);
						deg2rad(vy);
					}
					size_t nseg = vx.size() - 1;
					for (size_t i=0; i<np; i++) {
						if (d[i] != 0) { // && (inside[i] == 0)) {
							for (size_t j=0; j<nseg; j++) {
								d[i] = std::min(d[i],
									d2seg(x[i], y[i], vx[j], vy[j], vx[j+1], vy[j+1], r));
							}
						}
					}
				}
			}
			if ((method == "geo") && (m != 1)) {
				for (double& v : d) v *= m;
			}
		}
	} else if (type() == "lines") {
		for (size_t g=0; g<ng; g++) {
			size_t nparts = geoms[g].size();
			for (size_t h=0; h<nparts; h++) {
				vx = geoms[g].parts[h].x;
				vy = geoms[g].parts[h].y;
				if (method != "geo") {
					deg2rad(vx);
					deg2rad(vy);
				}
				size_t nseg = vx.size() - 1;
				for (size_t i=0; i<np; i++) {
					for (size_t j=0; j<nseg; j++) {
						d[i] = std::min(d[i],
							d2seg(x[i], y[i], vx[j], vy[j], vx[j+1], vy[j+1], r));
					}
				}
			}
		}
		if ((method == "geo") && (m != 1)) {
			for (double& v : d) v *= m;
		}
		
	} else { // if (type() == "points") {
		std::vector<std::vector<double>> pts = coordinates();
		if (method != "geo") {
			deg2rad(pts[0]);
			deg2rad(pts[1]);
		}
		d = pointdistance(x, y, pts[0], pts[1], false, m, true, method);
	}

	return d;
}



std::vector<double> SpatVector::distance(SpatVector x, bool pairwise, std::string unit, const std::string method) {

	std::vector<double> d;

	if (srs.is_empty() || x.srs.is_empty()) {
		setError("crs not defined");
		return(d);
	}
	if (! srs.is_same(x.srs, false) ) {
		setError("crs do not match");
		return(d);
	}

	size_t s = size();
	size_t sx = x.size();
	if ((s == 0) || (sx == 0)) {
		setError("empty SpatVector");
		return(d);
	}

	if (pairwise && (s != sx ) && (s > 1) && (sx > 1))  {
		setError("For pairwise distance, the number of geometries must match, or one should have a single geometry");
		return(d);
	}

	bool lonlat = is_lonlat();
	double m=1;
	if (!srs.m_dist(m, lonlat, unit)) {
		setError("invalid unit");
		return(d);
	}

	if ((method != "geo") && (method != "cosine")) {
		setError("invalid method. Must be 'geo' or 'cosine'");
		return(d);
	}

	std::string gtype = type();
	std::string xtype = x.type();

	if ((gtype == "points") && (xtype == "points")) {
		std::vector<std::vector<double>> p = coordinates();
		std::vector<std::vector<double>> px = x.coordinates();
		return pointdistance(p[0], p[1], px[0], px[1], pairwise, m, lonlat, method);
	} else if ((gtype == "points") || (xtype == "points")) {
		if (lonlat) {
			// not ok for multi-points
			if (gtype == "points") {
				std::vector<std::vector<double>> xy = coordinates();
				return x.nearestDistLonLat(xy[0], xy[1], unit, method);					
			} else {
				std::vector<std::vector<double>> xy = x.coordinates();
				return nearestDistLonLat(xy[0], xy[1], unit, method);					
			}
		} else {
			return geos_distance(x, pairwise, "", m);
		}
	} else {
		if (lonlat) {
			size_t n = size() * x.size();
			d.reserve(n);

			for (size_t i=0; i<size(); i++) {
				SpatVector tmp1 = subset_rows({(int)i});
				std::vector<std::vector<double>> xy1 = tmp1.coordinates();
				for (size_t j=0; j<x.size(); j++) {
					SpatVector tmp2 = x.subset_rows( {(int)j} );
					std::vector<double> d1 = tmp2.nearestDistLonLat(xy1[0], xy1[1], unit, method);

					std::vector<std::vector<double>> xy2 = tmp2.coordinates();
					std::vector<double> d2 = tmp1.nearestDistLonLat(xy2[0], xy2[1], unit, method);
					
					d.push_back(std::min(vmin(d1, false), vmin(d2, false)));
				}
			}
		} else {
			d = geos_distance(x, pairwise, "", m);
		}
	}
	return d;
}



std::vector<double> SpatVector::distance(bool sequential, std::string unit, const std::string method) {

	std::vector<double> d;
	if (srs.is_empty()) {
		setError("crs not defined");
		return(d);
	}

	bool lonlat = is_lonlat(); 
	double m=1;
	if (!srs.m_dist(m, lonlat, unit)) {
		setError("invalid unit");
		return(d);
	}
	std::string gtype = type();
	if (gtype == "points") {
		if (lonlat) {
			std::function<double(double, double, double, double)> dfun;
			if (method == "haversine") {
				dfun = distHaversine;
			} else if (method == "cosine") {
				dfun = distCosine;			
			} else if (method == "geo") {
				dfun = distLonlat;
			} else {
				setError("invalid lonlat distance method. Should be 'geo', 'cosine', or 'haversine'");
				return(d);	
			}

			if (sequential) {
				std::vector<std::vector<double>> p = coordinates();
				size_t n = p[0].size();
				d.reserve(n);
				d.push_back(0);
				n -= 1;
				if (lonlat) {
					for (size_t i=0; i<n; i++) {
						d.push_back(
							dfun(p[0][i], p[1][i], p[0][i+1], p[1][i+1]) *  m
						);
					}
				} else {
					for (size_t i=0; i<n; i++) {
						d.push_back(
							distance_plane(p[0][i], p[1][i], p[0][i+1], p[1][i+1]) * m
						);
					}
				}

			} else {
				size_t s = size();
				size_t n = ((s-1) * s)/2;
				d.reserve(n);
				std::vector<std::vector<double>> p = coordinates();
				if (lonlat) {			
					for (size_t i=0; i<(s-1); i++) {
						for (size_t j=(i+1); j<s; j++) {
							d.push_back(
								dfun(p[0][i], p[1][i], p[0][j], p[1][j]) * m
							);
						}
					}
				} else {
					for (size_t i=0; i<(s-1); i++) {
						for (size_t j=(i+1); j<s; j++) {
							d.push_back(
								distance_plane(p[0][i], p[1][i], p[0][j], p[1][j]) * m
							);
						}
					}
				}
			}
		} 
	} else {
		if (lonlat) {
			std::function<double(double, double, double, double)> dfun;
			if (method == "haversine") {
				dfun = distHaversine;
			} else if (method == "cosine") {
				dfun = distCosine;			
			} else if (method == "geo") {
				dfun = distLonlat;
			} else {
				setError("invalid lonlat distance method. Should be 'geo', 'cosine', or 'haversine'");
				return(d);	
			}

			size_t n = size();
			d.reserve(n);
			if (sequential) {
				n -= 1;
				SpatVector tmp1 = subset_rows({0});
				std::vector<std::vector<double>> xy1 = tmp1.coordinates();
				for (size_t i=0; i<n; i++) {
					SpatVector tmp2 = subset_rows( {(int)i+1} );
					std::vector<double> d1 = tmp2.nearestDistLonLat(xy1[0], xy1[1], unit, method);
					std::vector<std::vector<double>> xy2 = tmp2.coordinates();
					std::vector<double> d2 = tmp1.nearestDistLonLat(xy2[0], xy2[1], unit, method);
					d.push_back(std::min(vmin(d1, false), vmin(d2, false)));
					tmp1 = tmp2;
					xy1 = xy2;
				}
			} else {
				size_t s = size();
				size_t n = ((s-1) * s)/2;
				d.reserve(n);
				for (size_t i=0; i<(s-1); i++) {
					SpatVector tmp1 = subset_rows({(int)i});
					std::vector<std::vector<double>> xy1 = tmp1.coordinates();
					for (size_t j=(i+1); j<s; j++) {
						SpatVector tmp2 = subset_rows( {(int)j} );
						std::vector<double> d1 = tmp2.nearestDistLonLat(xy1[0], xy1[1], unit, method);
						std::vector<std::vector<double>> xy2 = tmp2.coordinates();
						std::vector<double> d2 = tmp1.nearestDistLonLat(xy2[0], xy2[1], unit, method);
						d.push_back(std::min(vmin(d1, false), vmin(d2, false)));
					}
				}
			}
		} else {
			return geos_distance(sequential, "", m);
		}
	}

	return d;
}


std::vector<double> SpatVector::pointdistance(const std::vector<double>& px, const std::vector<double>& py, const std::vector<double>& sx, const std::vector<double>& sy, bool pairwise, double m, bool lonlat, const std::string method) {

	std::vector<double> d;

	size_t szp = px.size();
	size_t szs = sx.size();
	if ((szp == 0) || (szs == 0)) {
		setError("empty SpatVector");
		return(d);
	}

	if (pairwise && (szp != szs ) && (szs > 1) && (szp > 1))  {
		setError("Can only do pairwise distance if geometries match, or if one is a single geometry");
		return(d);
	}

//	std::vector<std::vector<double>> p = coordinates();
//	std::vector<std::vector<double>> px = x.coordinates();


	size_t n = pairwise ? std::max(szs,szp) : szp*szs;
	d.reserve(n);

	std::function<double(double, double, double, double)> dfun;
	if (lonlat) {
		if (method == "haversine") {
			dfun = distHaversine;
		} else if (method == "cosine") {
			dfun = distCosine;			
		} else if (method == "geo") {
			dfun = distLonlat;
		} else {
			setError("invalid lonlat distance method. Should be 'geo', 'cosine', or 'haversine'");
			return(d);	
		}
	}

	if (pairwise) {
		if (szp == szs) {
			if (lonlat) {
				for (size_t i = 0; i < szs; i++) {
					d.push_back( dfun(px[i], py[i], sx[i], sy[i]) * m);
				}
			} else { // not reached
				for (size_t i = 0; i < szs; i++) {
					d.push_back( distance_plane(px[i], py[i], sx[i], sy[i]) * m);
				}
			}
		} else if (szp == 1) {  // to avoid recycling.
			if (lonlat) {
				for (size_t i = 0; i < szs; i++) {
					d.push_back( dfun(px[0], py[0], sx[i], sy[i]) * m);
				}
			} else { // not reached
				for (size_t i = 0; i < szs; i++) {
					d.push_back( distance_plane(px[0], py[0], sx[i], sy[i]) * m);
				}
			}
		} else { // if (szs == 1) {
			if (lonlat) {
				for (size_t i = 0; i < szp; i++) {
					d.push_back(dfun(px[i], py[i], sx[0], sy[0]) * m);
				}
			} else { // not reached
				for (size_t i = 0; i < szp; i++) {
					d.push_back(  distance_plane(px[i], py[i], sx[0], sy[0]) * m);
				}
			}
		}
	} else {
		if (lonlat) {
			for (size_t i=0; i<szp; i++) {
				for (size_t j=0; j<szs; j++) {
					d.push_back(dfun(px[i], py[i], sx[j], sy[j]) * m);
				}
			}
		} else { // not reached
			for (size_t i=0; i<szp; i++) {
				for (size_t j=0; j<szs; j++) {
					d.push_back(distance_plane(px[i], py[i], sx[j], sy[j]) * m);
				}
			}
		}
	}

	return d;
}

/*
std::vector<double> SpatVector::pointdistance_seq(const std::vector<double>& px, const std::vector<double>& py, double m, bool lonlat) {

	std::vector<double> d;
	size_t szp = px.size();
	d.reserve(szp);
	d.push_back(0);
	szp -= 1;

	if (lonlat) {
		for (size_t i = 0; i < szp; i++) {
			d.push_back( distance_lonlat(px[i], py[i], px[i+1], py[i+1]) );
		}
	} else { // not reached
		for (size_t i = 0; i < szs; i++) {
			d.push_back( distance_plane(px[i], py[i], px[i+1], py[i+1]) * m);
		}
	}
	return d;
}
*/




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

