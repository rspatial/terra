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
#include "crs.h"
//#include "sort.h"
#include "math_utils.h"

#if defined(USE_TBB)
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#endif 


double polDistLonLat(SpatVector &p1, SpatVector &p2, std::string unit, std::string method) {

	std::vector<int> inside = p1.relate(p2, "intersects", true, true);
	if (inside[0]) return 0;

	std::vector<std::vector<double>> xy = p2.coordinates();
	std::vector<double> x = xy[0];
	std::vector<double> y = xy[1];

	size_t ng = p1.size();
	size_t np = x.size();
	double d = std::numeric_limits<double>::infinity();

	double r = 6378137;
	double m = 1;
	if (unit == "km") {
		r = 6378.137;
		m = 0.001;
	}
	std::function<double(double,double,double,double,double,double,double)> d2seg;
	if (method != "geo") {
		deg2rad(x);
		deg2rad(y);
		d2seg = dist2segment_cos;
	} else {
		d2seg = dist2segment_geo;		
	}
	
	std::vector<double> dout;
	std::vector<double> vx, vy;


	if (p1.type() == "polygons") {
		for (size_t g=0; g<ng; g++) {
			size_t nparts = p1.geoms[g].size();
			for (size_t h=0; h<nparts; h++) {
				vx = p1.geoms[g].parts[h].x;
				vy = p1.geoms[g].parts[h].y;
				if (method != "geo") {
					deg2rad(vx);
					deg2rad(vy);
				}
				size_t nseg = vx.size() - 1;
				for (size_t i=0; i<np; i++) {
					if (d != 0) {
						for (size_t j=0; j<nseg; j++) {
							d = std::min(d, d2seg(x[i], y[i], vx[j], vy[j], vx[j+1], vy[j+1], r));
						}
					}
				}
				size_t nh = p1.geoms[g].parts[h].nHoles();
				for (size_t k=0; k < nh; k++) {
					vx = p1.geoms[g].parts[h].holes[k].x;
					vy = p1.geoms[g].parts[h].holes[k].y;
					if (method != "geo") {
						deg2rad(vx);
						deg2rad(vy);
					}
					size_t nseg = vx.size() - 1;
					for (size_t i=0; i<np; i++) {
						if (d != 0) {
							for (size_t j=0; j<nseg; j++) {
								d = std::min(d, d2seg(x[i], y[i], vx[j], vy[j], vx[j+1], vy[j+1], r));
							}
						}
					}
				}
			}
		}
	} else if (p1.type() == "lines") {
		for (size_t g=0; g<ng; g++) {
			size_t nparts = p1.geoms[g].size();
			for (size_t h=0; h<nparts; h++) {
				vx = p1.geoms[g].parts[h].x;
				vy = p1.geoms[g].parts[h].y;
				if (method != "geo") {
					deg2rad(vx);
					deg2rad(vy);
				}
				size_t nseg = vx.size() - 1;
				for (size_t i=0; i<np; i++) {
					for (size_t j=0; j<nseg; j++) {
						d = std::min(d, d2seg(x[i], y[i], vx[j], vy[j], vx[j+1], vy[j+1], r));
					}
				}
			}
		}


	} else { // if (type() == "points") {
		return -1;
	}



	if ((method == "geo") && (m != 1)) {
		d *= m;
	}

	return d;
}





std::vector<double> SpatVector::distLonLat(SpatVector p, std::string unit, std::string method, bool transp) {

	std::vector<std::vector<double>> xy = p.coordinates();
	std::vector<double> x = xy[0];
	std::vector<double> y = xy[1];

	size_t np = x.size();
	size_t ng = size();
	double inf = std::numeric_limits<double>::infinity();
	std::vector<std::vector<double>> d(np, std::vector<double>(ng, inf));

/*
	std::vector<int> inside = relate(p, "intersects", true, true);
	Rcpp::Rcout << inside.size() << " " << ng << " " << np << std::endl;
	
	for (size_t i=0; i<ng; i++) {
		for (size_t j=0; j<np; j++) {
			if (inside[i*np+j]) {
				d[j][i] = 0;
			}
		}
	}
*/
	if (type() == "polygons") {
		std::vector<int> inside = pointInPolygon(x, y);
		for (size_t i=0; i<ng; i++) {
			for (size_t j=0; j<np; j++) {
				if (inside[i*np+j]) {
					d[j][i] = 0;
				}
			}
		}
	}


	double r = 6378137;
	double m = 1;
	if (unit == "km") {
		r = 6378.137;
		m = 0.001;
	}
	std::function<double(double,double,double,double,double,double,double)> d2seg;
	if (method != "geo") {
		deg2rad(x);
		deg2rad(y);
		d2seg = dist2segment_cos;
	} else {
		d2seg = dist2segment_geo;		
	}
	
	std::vector<double> dout;
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
					if (d[i][g] != 0) {
						for (size_t j=0; j<nseg; j++) {
							d[i][g] = std::min(d[i][g],
								d2seg(x[i], y[i], vx[j], vy[j], vx[j+1], vy[j+1], r));
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
						if (d[i][g] != 0) {
							for (size_t j=0; j<nseg; j++) {
								d[i][g] = std::min(d[i][g],
									d2seg(x[i], y[i], vx[j], vy[j], vx[j+1], vy[j+1], r));
							}
						}
					}
				}
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
						d[i][g] = std::min(d[i][g],
							d2seg(x[i], y[i], vx[j], vy[j], vx[j+1], vy[j+1], r));
					}
				}
			}
		}


	} else { // if (type() == "points") {
		std::vector<std::vector<double>> pts = coordinates();
		if (method != "geo") {
			deg2rad(pts[0]);
			deg2rad(pts[1]);
		}
		return pointdistance(x, y, pts[0], pts[1], false, m, true, method);
	}


		dout.reserve(np*ng);
		if (transp) {
			for (size_t i=0; i<d[0].size(); i++) {
				for (size_t j=0; j<d.size(); j++) {
					dout.push_back(d[j][i]);				
				}
			}				
		} else {
			size_t j = 0;
			for (size_t i=0; i<d.size(); i++) {
				dout.insert(dout.begin()+j, d[i].begin(), d[i].end());				
				j += d[i].size();
			}	
		}

	if ((method == "geo") && (m != 1)) {
		for (double& v : dout) v *= m;
	}

	return dout;
}


/*
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
}
*/

std::vector<double> SpatVector::distance(SpatVector x, bool pairwise, std::string unit, const std::string method, bool use_nodes, SpatOptions &opt) {

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
	bool lonlat = is_lonlat();
	double m=1;
	if (!srs.m_dist(m, lonlat, unit)) {
		setError("invalid unit");
		return(d);
	}
	if (pairwise && (s != sx ) && (s > 1) && (sx > 1))  {
		setError("For pairwise distance, the number of geometries must match, or one should have a single geometry");
		return(d);
	}


	if ((method != "geo") && (method != "cosine") && (method != "haversine")) {
		setError("invalid method. Must be 'geo', 'haversine' or 'cosine'");
		return(d);
	}

	if (lonlat && use_nodes) {
		SpatVector p = as_points(true, false);
		SpatVector xp = x.as_points(true, false);
		std::vector<double> out;
		if (pairwise) {
			out.reserve(p.size());
			for (size_t i=0; i<p.size(); i++) {
				SpatVector pi = p.subset_rows(i);	
				SpatVector xpi = xp.subset_rows(i);	
				std::vector<std::vector<double>> pc = pi.coordinates();
				std::vector<std::vector<double>> xpc = xpi.coordinates();
				std::vector<double> d = pointdistance(pc[0], pc[1], xpc[0], xpc[1], false, m, lonlat, method);
				auto min_it = std::min_element(d.begin(), d.end());
				out.push_back(*min_it);
			}			
		} else {
			out.reserve(p.size() * xp.size());
			for (size_t i=0; i<p.size(); i++) {
				SpatVector pi = p.subset_rows(i);	
				std::vector<std::vector<double>> pc = pi.coordinates();
				for (size_t j=0; j<xp.size(); j++) {
					SpatVector xpj = xp.subset_rows(j);	
					std::vector<std::vector<double>> xpc = xpj.coordinates();
					std::vector<double> d = pointdistance(pc[0], pc[1], xpc[0], xpc[1], false, m, lonlat, method);
					auto min_it = std::min_element(d.begin(), d.end());
					out.push_back(*min_it);
				}
			}
		}
		return out;
	}



	std::string gtype = type();
	std::string xtype = x.type();

	if ((gtype == "points") && (xtype == "points")) {
		std::vector<std::vector<double>> p = coordinates();
		std::vector<std::vector<double>> px = x.coordinates();
		return pointdistance(p[0], p[1], px[0], px[1], pairwise, m, lonlat, method);
	} else if ((gtype == "points") || (xtype == "points")) {
		if (lonlat) {
			// not yet ok for multi-points
			if (gtype == "points") {
//				std::vector<std::vector<double>> xy = coordinates();
				if (pairwise) {
					std::vector<double> out;
					out.reserve(size());
					for (size_t i=0; i<size(); i++) {
						SpatVector p = subset_rows(i);	
						SpatVector xp = x.subset_rows(i);	
						std::vector<double> d = xp.distLonLat(p, unit, method, false);	
						auto min_it = std::min_element(d.begin(), d.end());
						out.push_back(*min_it);
					} 
					return out;					
				} else {
					return x.distLonLat(*this, unit, method, false);	
				}
			} else {
//				std::vector<std::vector<double>> xy = x.coordinates();
				if (pairwise) {
					std::vector<double> out;
					out.reserve(size());
					for (size_t i=0; i<size(); i++) {
						SpatVector p = subset_rows(i);	
						SpatVector xp = x.subset_rows(i);	
						std::vector<double> d = p.distLonLat(xp, unit, method, false);	
						auto min_it = std::min_element(d.begin(), d.end());
						out.push_back(*min_it);
					} 
					return out;					
				} else {
					return distLonLat(x, unit, method, true);					
				}
			}
		} else {
			return geos_distance(x, pairwise, "", m);
		}
	} else {
		if (lonlat) {

			if (pairwise) {
				d.reserve(size());
				for (size_t i=0; i<size(); i++) {
					SpatVector p = subset_rows(i);	
					SpatVector xp = x.subset_rows(i);	
					double d1 = polDistLonLat(p, xp, unit, method);	
					double d2 = polDistLonLat(xp, p, unit, method);
					d.push_back(std::min(d1, d2));
				} 
				return d;
			}

			size_t n = size() * x.size();
			d.reserve(n);

/*
			std::vector<std::vector<double>> e1, e2;
			e1.reserve(n);
			for (size_t g=0; g<n; g++) {
				e1.push_back(geoms[g].extent.asVector());
			}
			e2.reserve(x.size());
			for (size_t g=0; g<x.size(); g++) {
				e2.push_back(x.geoms[g].extent.asVector());
			}
			std::vector<std::vector<size_t>> idx = get_index(e1, e2);
*/

			for (size_t i=0; i<size(); i++) {
				SpatVector tmp1 = subset_rows( (long)i);
//				std::vector<std::vector<double>> xy1 = tmp1.coordinates();
				for (size_t j=0; j<x.size(); j++) {
					SpatVector tmp2 = x.subset_rows(long(j));
//					std::vector<double> d1 = tmp2.distLonLat(xy1[0], xy1[1], unit, method, false);
//					std::vector<std::vector<double>> xy2 = tmp2.coordinates();
//					std::vector<double> d2 = tmp1.distLonLat(xy2[0], xy2[1], unit, method, false);
//					std::vector<double> d1 = tmp2.distLonLat(tmp1, unit, method, false);
//					std::vector<double> d2 = tmp1.distLonLat(tmp2, unit, method, false);
//					d.push_back(std::min(vmin(d1, false), vmin(d2, false)));

					double d1 = polDistLonLat(tmp2, tmp1, unit, method);	
					double d2 = polDistLonLat(tmp1, tmp2, unit, method);
					
					d.push_back(std::min(d1, d2));
				}
			}
		} else {
			d = geos_distance(x, pairwise, "", m);
		}
	}
	return d;
}



// distance to self
std::vector<double> SpatVector::distance(bool sequential, std::string unit, const std::string method, bool use_nodes, SpatOptions &opt) {

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

	if (lonlat && use_nodes) {
		SpatVector p = as_points(true, false);
		std::vector<double> out;
		if (sequential) {
			out.reserve(p.size()-1);
			SpatVector pi = p.subset_rows(0);	
			std::vector<std::vector<double>> pic = pi.coordinates();
			for (size_t i=0; i<(p.size()-1); i++) {
				SpatVector pj = p.subset_rows(i+1);	
				std::vector<std::vector<double>> pjc = pj.coordinates();
				std::vector<double> d = pointdistance(pic[0], pic[1], pjc[0], pjc[1], false, m, lonlat, method);
				auto min_it = std::min_element(d.begin(), d.end());
				out.push_back(*min_it);
				pic = pjc;
			}
		} else {
			out.reserve((p.size() * (p.size()-1)) / 2);
			for (size_t i=0; i<(p.size()-1); i++) {
				SpatVector pi = p.subset_rows(i);	
				std::vector<std::vector<double>> pic = pi.coordinates();
				for (size_t j=(i+1); j<p.size(); j++) {
					SpatVector pj = p.subset_rows(j);	
					std::vector<std::vector<double>> pjc = pj.coordinates();
					std::vector<double> d = pointdistance(pic[0], pic[1], pjc[0], pjc[1], false, m, lonlat, method);
					auto min_it = std::min_element(d.begin(), d.end());
					out.push_back(*min_it);
				}
			}
		}
		return out;
	}
	
	std::string gtype = type();
	std::function<double(double, double, double, double)> dfun;
	if (gtype == "points") {
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
		if (sequential) {
			std::vector<std::vector<double>> p = coordinates();
			size_t n = p[0].size();
			if (lonlat) {
#if defined(USE_TBB)
				if (opt.parallel) {
					d.resize(n);
					tbb::parallel_for(tbb::blocked_range<size_t>(0, n-1),
					[&](const tbb::blocked_range<size_t>& range) {
						for (size_t i = range.begin(); i != range.end(); i++) {
							d[i+1] = dfun(p[0][i], p[1][i], p[0][i+1], p[1][i+1]) * m;
						}
					});
				} else {
					d.reserve(n);
					d.push_back(0);
					n -= 1;
					for (size_t i=0; i<n; i++) {
						d.push_back( dfun(p[0][i], p[1][i], p[0][i+1], p[1][i+1]) *  m );
					}
				}
#else
				d.reserve(n);
				d.push_back(0);
				n -= 1;
				for (size_t i=0; i<n; i++) {
					d.push_back( dfun(p[0][i], p[1][i], p[0][i+1], p[1][i+1]) *  m );
				}
#endif
			} else {
				d.reserve(n);
				d.push_back(0);
				n -= 1;
				for (size_t i=0; i<n; i++) {
					d.push_back( distance_plane(p[0][i], p[1][i], p[0][i+1], p[1][i+1]) * m );
				}
			}
		} else {
			size_t s = size();
			size_t n = ((s-1) * s)/2;
			d.reserve(n);
			std::vector<std::vector<double>> p = coordinates();
			if (lonlat) {			
#if defined(USE_TBB)
				if (opt.parallel) {
					d.resize(n);
					tbb::parallel_for(tbb::blocked_range<size_t>(0, s-2),
					[&](const tbb::blocked_range<size_t>& range) {
						for (size_t i = range.begin(); i != range.end(); i++) {
							size_t k = 0;
							for (size_t j=0; j<i; j++) {
								k += s-1-j;
							}
							for (size_t j=(i+1); j<s; j++) {
								d[k+j-i-1] = dfun(p[0][i], p[1][i], p[0][j], p[1][j]) * m;
							}
						}
					});
				} else {
					for (size_t i=0; i<(s-1); i++) {
						for (size_t j=(i+1); j<s; j++) {
							d.push_back(dfun(p[0][i], p[1][i], p[0][j], p[1][j]) * m);
						}
					}
				}
#else
				for (size_t i=0; i<(s-1); i++) {
					for (size_t j=(i+1); j<s; j++) {
						d.push_back( dfun(p[0][i], p[1][i], p[0][j], p[1][j]) * m );
					}
				}
#endif
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
	} else {
		if (lonlat) {
			if (method == "cosine") {
				dfun = distCosine;			
			} else if (method == "haversine") {
				dfun = distHaversine;
			} else if (method == "geo") {
				dfun = distLonlat;
			} else {
				setError("invalid lonlat distance method. Should be 'geo' or 'cosine' for lines and polygons");
				return(d);	
			}
			size_t n = size();
			d.reserve(n);


			if (sequential) {

				n -= 1;
//				std::vector<std::vector<size_t>> idx;
#if defined(USE_TBB)
				if (opt.parallel) {
					d.resize(n);			
					tbb::parallel_for(tbb::blocked_range<size_t>(0, n),
					[&](const tbb::blocked_range<size_t>& range) {
						for (size_t i = range.begin(); i != range.end(); i++) {
							SpatVector tmp1 = subset_rows((long)i);
							SpatVector tmp2 = subset_rows((long)i+1);
							double d1 = polDistLonLat(tmp2, tmp1, unit, method);
							double d2 = polDistLonLat(tmp1, tmp2, unit, method);
							d[i] = std::min(d1, d2);
						}
					});
				} else {
					SpatVector tmp1 = subset_rows(0);
					for (size_t i=0; i<n; i++) {
						SpatVector tmp2 = subset_rows( (long)i+1 );
						double d1 = polDistLonLat(tmp2, tmp1, unit, method);
						double d2 = polDistLonLat(tmp1, tmp2, unit, method);
						d.push_back(std::min(d1, d2));
						tmp1 = tmp2;
					}
				}
#else
				SpatVector tmp1 = subset_rows(0);
				for (size_t i=0; i<n; i++) {
					SpatVector tmp2 = subset_rows( (long)i+1 );
					double d1 = polDistLonLat(tmp2, tmp1, unit, method);
					double d2 = polDistLonLat(tmp1, tmp2, unit, method);
					d.push_back(std::min(d1, d2));
					tmp1 = tmp2;
				}
#endif				
			} else {  // not sequential
				size_t s = size();
				size_t n = ((s-1) * s)/2;
				d.reserve(n);
								
				std::vector<double> dst;
				for (size_t i=0; i<(s-1); i++) {
					SpatVector tmp1 = subset_rows(long(i));
					dst.resize(s-i-1);

#if defined(USE_TBB)
					if (opt.parallel) {
						tbb::parallel_for(tbb::blocked_range<size_t>((i+1), s),
						[&](const tbb::blocked_range<size_t>& range) {
							for (size_t j = range.begin(); j != range.end(); j++) {
								SpatVector tmp2 = subset_rows( long(j) );
								double d1 = polDistLonLat(tmp2, tmp1, unit, method);
								double d2 = polDistLonLat(tmp1, tmp2, unit, method);
								dst[j-i-1] = std::min(d1, d2);
							}
						});
					} else {
						for (size_t j=(i+1); j<s; j++) {
							SpatVector tmp2 = subset_rows( long(j) );
							double d1 = polDistLonLat(tmp2, tmp1, unit, method);
							double d2 = polDistLonLat(tmp1, tmp2, unit, method);
							dst[j-i-1] = std::min(d1, d2);
						}						
					}
#else
					for (size_t j=(i+1); j<s; j++) {
						SpatVector tmp2 = subset_rows( long(j) );
						double d1 = polDistLonLat(tmp2, tmp1, unit, method);
						double d2 = polDistLonLat(tmp1, tmp2, unit, method);
						dst[j-i-1] = std::min(d1, d2);
					}
#endif
					d.insert(d.end(), dst.begin(), dst.end());

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




void split_dateline(SpatVector &v) {
	SpatExtent e1 = {-1,  180, -91, 91};
	SpatExtent e2 = {180, 361, -91, 91};
	SpatVector ve(e1, "");
	SpatVector ve2(e2, "");
	ve = ve.append(ve2, true);
	v = v.intersect(ve, true);
	ve = v.subset_rows(1);
	ve = ve.shift(-360, 0);
	v.geoms[1] = ve.geoms[0];
	v = v.aggregate(false);
}




bool fix_date_line(SpatGeom &g, std::vector<double> &x, const std::vector<double> &y) {

	SpatPart p(x, y);
	double minx = vmin(x, false);
	double maxx = vmax(x, false);
	// need a better check but this should work for all normal cases
	if ((maxx - minx) > 180) {
		for (size_t i=0; i<x.size(); i++) {
			if (x[i] < 0) {
				x[i] += 360;
			}
		}
		double minx2 = vmin(x, false);
		double maxx2 = vmax(x, false);
		if ((maxx - minx) < (maxx2 - minx2)) {
			g.reSetPart(p);
			return false;
		}
		p.x = x;
		g.reSetPart(p);
		SpatVector v(g);
		split_dateline(v);
		g = v.geoms[0];
		return true;
	}
	g.reSetPart(p);
	return false;
}


SpatVector SpatVector::point_buffer(std::vector<double> d, unsigned quadsegs, bool no_multipolygons, bool wrap) {

	SpatVector out;
	out.reserve(size());
	std::string vt = type();
	if (vt != "points") {
		out.setError("geometry must be points");
		return out;
	}

	size_t npts = size();

	size_t n = quadsegs * 4;
	double step = 360.0 / n;
	SpatGeom g(polygons);
	g.addPart(SpatPart(0, 0));

//  not good for multipoints
//	std::vector<std::vector<double>> xy = coordinates();
	
	if (is_lonlat()) {
		
		std::vector<double> gptx = std::vector<double> {-180,  0, 180, 180, 180,   0, -180, -180, -180};
		std::vector<double> gpty = std::vector<double> {  90, 90,  90,   0, -90, -90,  -90,    0,   90};
		SpatGeom ggeom(polygons);
		ggeom.addPart(SpatPart(gptx, gpty));
		SpatVector glob;
		glob.addGeom(ggeom);	

		
		std::vector<double> brng(n);
		for (size_t i=0; i<n; i++) {
			brng[i] = i * step;
		}
		double a = 6378137.0;
		double f = 1/298.257223563;
		struct geod_geodesic gd;
		geod_init(&gd, a, f);
		double lat, lon, azi, s12, azi2;

		for (size_t p=0; p<npts; p++) { 
			std::vector<std::vector<double>> xy = geoms[p].coordinates();
			SpatVector tmp;
			for (size_t i=0; i<xy[0].size(); i++) {
				if (std::isnan(xy[0][i]) || std::isnan(xy[1][i]) || (xy[1][i] > 90) || (xy[1][i] < -90)) {
					tmp.addGeom(SpatGeom(polygons));
				} else if (d[p] > 20003931) {
					tmp = glob;
					break;
				} else {
					std::vector<double> ptx;
					std::vector<double> pty;
					ptx.reserve(n+1);
					pty.reserve(n+1);
					if (wrap) {
						for (size_t j=0; j < n; j++) {
							geod_direct(&gd, xy[1][i], xy[0][i], brng[j], d[p], &lat, &lon, &azi);
							ptx.push_back(lon);
							pty.push_back(lat);
						}
					} else {
						for (size_t j=0; j < n; j++) {
							geod_direct(&gd, xy[1][i], 0, brng[j], d[p], &lat, &lon, &azi);
							ptx.push_back(lon+xy[0][i]);
							pty.push_back(lat);
						}
					}

					geod_inverse(&gd, xy[1][i], xy[0][i],  90, xy[0][i], &s12, &azi, &azi2);
					bool npole = s12 < d[p];
					geod_inverse(&gd, xy[1][i], xy[0][i], -90, xy[0][i], &s12, &azi, &azi2);
					bool spole = s12 < d[p];

					if (npole && spole) {
						ptx.push_back(ptx[0]);
						pty.push_back(pty[0]);
						bool split = false;
						if (vmax(ptx, true) >= 0) {
							for (size_t i=0; i<ptx.size(); i++) {
								if (ptx[i] < 0) {
									ptx[i] += 360;
									split = true;
								}
							}
						}
						g.reSetPart(SpatPart(ptx, pty));
						tmp.addGeom(g);
						if (split) {
							split_dateline(tmp);
						} 	
						tmp = glob.erase(tmp);
					} else if (npole) {
						sort_unique_2d(ptx, pty);
						if (ptx[ptx.size()-1] < 180) {
								ptx.push_back(180);
								pty.push_back(pty[pty.size()-1]);
						}
						ptx.push_back(180);
						pty.push_back(90);
						ptx.push_back(-180);
						pty.push_back(90);
						if (ptx[0] > -180) {
							ptx.push_back(-180);
							pty.push_back(pty[0]);
						}
						ptx.push_back(ptx[0]);
						pty.push_back(pty[0]);
						g.reSetPart(SpatPart(ptx, pty));
						tmp.addGeom(g);
					} else if (spole) {
						sort_unique_2d(ptx, pty);
						if (ptx[ptx.size()-1] < 180) {
							ptx.push_back(180);
							pty.push_back(pty[pty.size()-1]);
						}
						ptx.push_back(180);
						pty.push_back(-90);
						ptx.push_back(-180);
						pty.push_back(-90);
						if (ptx[0] > -180) {
							ptx.push_back(-180);
							pty.push_back(pty[0]);
						}
						ptx.push_back(ptx[0]);
						pty.push_back(pty[0]);
						g.reSetPart(SpatPart(ptx, pty));
						tmp.addGeom(g);
					} else {
						ptx.push_back(ptx[0]);
						pty.push_back(pty[0]);
						if (wrap) {
							bool split = false;
							try {
								split = fix_date_line(g, ptx, pty);
							} catch(...) {}
							
							if (split & no_multipolygons) {
								for (size_t j=0; j<g.parts.size(); j++) {
									SpatGeom gg(g.parts[j], polygons);
									tmp.addGeom(gg);
								}
							} else {
								tmp.addGeom(g);
							}
						} else {
							g.reSetPart(SpatPart(ptx, pty));
							tmp.addGeom(g);		
						}
					}	
				}
			}
			if (tmp.size() > 1) {
				tmp = tmp.aggregate(true);
			}
			out.addGeom(tmp.geoms[0]);
		}
		
	} else { // not used (GEOS used for planar). Would need to be fixed for multipoints
		std::vector<std::vector<double>> xy = coordinates();

		std::vector<double> cosb(n);
		std::vector<double> sinb(n);
		std::vector<double> px(n+1);
		std::vector<double> py(n+1);
		for (size_t i=0; i<n; i++) {
			double brng = i * step;
			brng = toRad(brng);
			cosb[i] = d[i] * cos(brng);
			sinb[i] = d[i] * sin(brng);
		}
		for (size_t i=0; i<npts; i++) {
			if (std::isnan(xy[0][i]) || std::isnan(xy[1][i])) {
				out.addGeom(SpatGeom(polygons));
			} else {
				for (size_t j=0; j<n; j++) {
					px[j] = xy[0][i] + cosb[j];
					py[j] = xy[1][i] + sinb[j];
				}
				px[n] = px[0];
				py[n] = py[0];
				g.setPart(SpatPart(px, py), 0);
				out.addGeom(g);
			}
		}
	}
	out.srs = srs;
	out.df = df;
	return(out);
}



SpatGeom hullify(SpatVector b, bool ispoly) {
	if (b.nrow() == 1) return b.geoms[0];
	if (ispoly) b.addGeom(b.geoms[0]);
	SpatVector part;
	part.reserve(b.size());
	for (size_t j =0; j<(b.size()-1); j++) {
		std::vector<size_t> range = {j, j+1};
		SpatVector g = b.subset_rows(range);
		g = g.hull("convex");
		part.addGeom(g.geoms[0]);
	}
	part = part.aggregate(true);
	return part.geoms[0];
}


SpatVector lonlat_buf(SpatVector x, double dist, unsigned quadsegs, bool ispol, bool ishole) {

/*
	if ((x.extent.ymin > -60) && (x.extent.ymax < 60) && 
			((x.extent.ymax - x.extent.ymin) < 1) && dist < 110000) {
				
		SpatSRS insrs = x.srs;
		x.setSRS("+proj=merc");
		double f = 0.5 - (dist / 220000);
		double halfy = x.extent.ymin + f * (x.extent.ymax - x.extent.ymin);
		std::vector<double> dd = destpoint_lonlat(0, halfy, 0, dist);
		dist = dd[1] - halfy;
		if (ishole) dist = -dist;
		x = x.buffer({dist}, quadsegs, "", "", NAN, false);	
		x.srs = insrs;
		return x;
	}
*/
	x = x.disaggregate(false);
	SpatVector tmp;
	tmp.reserve(x.size());
	//Rcpp::Rcout << x.geoms.size() << std::endl;

//	double interval = std::max(1000.0, std::max(x.extent.ymax - x.extent.ymin, x.extent.xmax - x.extent.xmin) * 100);  // m
	for (size_t i=0; i<x.geoms.size(); i++) {
		SpatVector p(x.geoms[i]);
//		p = p.densify(interval, true, false);

		p.srs = x.srs;
		p = p.as_points(false, true);
		std::vector<double> d(p.size(), dist);
		SpatVector b = p.point_buffer(d, quadsegs, true, false);
		if (b.size() <= p.size()) {
			SpatExtent e = b.extent;
			if ((e.xmin > -180) || (e.xmax < 180)) {
				SpatGeom g = hullify(b, ispol);
				tmp.addGeom(g);
			} else {
				b = b.aggregate(true);
				b = b.remove_holes();
				tmp.addGeom(b.geoms[0]);				
			}
		} else {
			SpatVector west, east, eastwest;
			for (size_t j =0; j<b.size(); j++) {
				if ((b.geoms[j].extent.xmin < -179.99) && (b.geoms[j].extent.xmax > 179.99)) {
					tmp.addGeom(b.geoms[j]);
				} else if (b.geoms[j].extent.xmax < 0) {
					west.addGeom(b.geoms[j]);
				} else {
					east.addGeom(b.geoms[j]);
				}
			}
			if (east.nrow() > 0) {
				SpatGeom geast = hullify(east, ispol);
				tmp.addGeom(geast);
			}
			if (west.nrow() > 0) {
				SpatGeom gwest = hullify(west, ispol);
				tmp.addGeom(gwest);
			}
		}
	}
	tmp = tmp.aggregate(true);

	tmp.fix_lonlat_overflow();
	
	if (ispol) {
		if (dist < 0) {
			tmp = !ishole ? tmp.get_holes() : tmp.remove_holes();
		} else {
			tmp = ishole ? tmp.get_holes() : tmp.remove_holes();
		}
	}
	return tmp;
}


SpatVector SpatVector::buffer_lonlat(std::string vt, std::vector<double> d, unsigned quadsegs) {

	SpatVector out;
	std::vector<size_t> keep;
	keep.reserve(size());
	if (vt == "points") {
		return point_buffer(d, quadsegs, false, true);
	} else if (vt == "polygons") {
		for (size_t i =0; i<size(); i++) {
			SpatVector p(geoms[i]);
			p = p.disaggregate(false);
			SpatVector tmp;
			for (size_t j =0; j<p.size(); j++) {
				SpatVector pp(p.geoms[j]);			
				pp.srs = srs;
				SpatVector h = pp.get_holes();
				pp = pp.remove_holes();
				pp = lonlat_buf(pp, d[i], quadsegs, true, false);
				if (!(pp.empty() || h.empty())) {
					h = lonlat_buf(h, d[i], quadsegs, true, true);
					if (!h.empty()) {
						if (d[i] < 0) {
							pp = pp.erase(h);
							if (pp.empty()) continue;
							h = h.crop(pp);
							if (h.empty()) continue;
						}
						for (size_t k=0; k<h.geoms[0].parts.size(); k++) {
							pp.geoms[0].parts[0].addHole(h.geoms[0].parts[k].x, 
														 h.geoms[0].parts[k].y);
						}
					}
				}
				tmp = tmp.append(pp, true);
			}
			if (!tmp.empty()) {
				tmp = tmp.aggregate(true);
				keep.push_back(i);
				out = out.append(tmp, true);
			}
		}
		if (keep.size() < size()) {
			out.df = df.subset_rows(keep);
		} else {
			out.df = df;
		}
	} else {
		for (size_t i =0; i<size(); i++) {
			SpatVector p(geoms[i]);
			p.srs = srs;
			p = lonlat_buf(p, d[i], quadsegs, false, false);
			out = out.append(p, true);
		}
		out.df = df;
	}
	out.srs = srs;
	return out;
}



double area_polygon_lonlat(geod_geodesic &g, const std::vector<double> &lon, const std::vector<double> &lat) {
	struct geod_polygon p;
	geod_polygon_init(&p, 0);
	size_t n = lat.size();
	for (size_t i=0; i < n; i++) {
		//double lat = lat[i] > 90 ? 90 : lat[i] < -90 ? -90 : lat[i];
		// for #397
		double flat = lat[i] < -90 ? -90 : lat[i];
		geod_polygon_addpoint(&g, &p, flat, lon[i]);
	}
	double area, P;
	geod_polygon_compute(&g, &p, 0, 1, &area, &P);
	return(area < 0 ? -area : area);
}



double area_polygon_plane(std::vector<double> x, std::vector<double> y) {
// based on http://paulbourke.net/geometry/polygonmesh/source1.c
	size_t n = x.size();
	double area = x[n-1] * y[0];
	area -= y[n-1] * x[0];
	for (size_t i=0; i < (n-1); i++) {
		area += x[i] * y[i+1];
		area -= x[i+1] * y[i];
	}
	area /= 2;
	return(area < 0 ? -area : area);
}


double area_lonlat(geod_geodesic &g, const SpatGeom &geom) {
	double area = 0;
	if (geom.gtype != polygons) return area;
	for (size_t i=0; i<geom.parts.size(); i++) {
		area += area_polygon_lonlat(g, geom.parts[i].x, geom.parts[i].y);
		for (size_t j=0; j < geom.parts[i].holes.size(); j++) {
			area -= area_polygon_lonlat(g, geom.parts[i].holes[j].x, geom.parts[i].holes[j].y);
		}
	}
	return area;
}


double area_plane(const SpatGeom &geom) {
	double area = 0;
	if (geom.gtype != polygons) return area;
	for (size_t i=0; i < geom.parts.size(); i++) {
		area += area_polygon_plane(geom.parts[i].x, geom.parts[i].y);
		for (size_t j=0; j < geom.parts[i].holes.size(); j++) {
			area -= area_polygon_plane(geom.parts[i].holes[j].x, geom.parts[i].holes[j].y);
		}
	}
	return area;
}


std::vector<double> SpatVector::area(std::string unit, bool transform, std::vector<double> mask) {

	if (type() != "polygons") { // area is zero
		std::vector<double>	out(nrow(), 0);
		return out;
	}
	if (nrow() == 0) {
		std::vector<double>	out(1, 0);	
		return out;
	}

	size_t s = size();
	size_t m = mask.size();
	bool domask = false;
	if (m > 0) {
		if (s != mask.size()) {
			addWarning("mask size is not correct");
		} else {
			domask = true;
		}
	}

	std::vector<double> ar;
	ar.reserve(s);

	std::vector<std::string> ss {"m", "km", "ha"};
	if (std::find(ss.begin(), ss.end(), unit) == ss.end()) {
		setError("invalid unit");
		return {NAN};
	}
	double adj = unit == "m" ? 1 : unit == "km" ? 1000000 : 10000;

	if (srs.wkt.empty()) {
		addWarning("unknown CRS. Results can be wrong");
		if (domask) {
			for (size_t i=0; i<s; i++) {
				if (std::isnan(mask[i])) {
					ar.push_back(NAN);
				} else {
					ar.push_back(area_plane(geoms[i]));
				}
			}
		} else {
			for (size_t i=0; i<s; i++) {
				ar.push_back(area_plane(geoms[i]));
			}
		}
	} else {
		if (!srs.is_lonlat()) {
			if (transform) {
				transform = can_transform(srs.wkt, "EPSG:4326");
			}	
			if (transform) {
				SpatVector v = project("EPSG:4326", false);
				return v.area(unit, false, mask);
			} else {
				transform = false;
				double m = srs.to_meter();
				adj *= std::isnan(m) ? 1 : m * m;
				if (domask) {
					for (size_t i=0; i<s; i++) {
						if (std::isnan(mask[i])) {
							ar.push_back(NAN);
						} else {
							ar.push_back(area_plane(geoms[i]));
						}
					}
				} else {
					for (size_t i=0; i<s; i++) {
						ar.push_back(area_plane(geoms[i]));
					}
				}
			}

		} else {

			struct geod_geodesic g;
			double a = 6378137;
			double f = 1 / 298.257223563;
			geod_init(&g, a, f);
			if (domask) {
				for (size_t i=0; i<s; i++) {
					if (std::isnan(mask[i])) {
						ar.push_back(NAN);
					} else {
						ar.push_back(area_lonlat(g, geoms[i]));
					}
				}
			} else {
				for (size_t i=0; i<s; i++) {
					ar.push_back(area_lonlat(g, geoms[i]));
				}
			}
		}
	}

	if (adj != 1) {
		for (double& i : ar) i /= adj;
	}
	return ar;
}



double length_line_lonlat(geod_geodesic &g, const std::vector<double> &lon, const std::vector<double> &lat) {
	size_t n = lat.size();
	double length = 0;
	for (size_t i=1; i < n; i++) {
		length += distance_lonlat(lon[i-1], lat[i-1], lon[i], lat[i]);
	}
	return (length);
}


double length_line_plane(std::vector<double> x, std::vector<double> y) {
	size_t n = x.size();
	double length = 0;
	for (size_t i=1; i<n; i++) {
		length += sqrt(pow(x[i-1] - x[i], 2) + pow(y[i-1] - y[i], 2));
	}
	return (length);
}


double length_lonlat(geod_geodesic &g, const SpatGeom &geom) {
	double length = 0;
	if (geom.gtype == points) return length;
	for (size_t i=0; i<geom.parts.size(); i++) {
		length += length_line_lonlat(g, geom.parts[i].x, geom.parts[i].y);
		for (size_t j=0; j<geom.parts[i].holes.size(); j++) {
			length += length_line_lonlat(g, geom.parts[i].holes[j].x, geom.parts[i].holes[j].y);
		}
	}
	return length;
}


double length_plane(const SpatGeom &geom) {
	double length = 0;
	if (geom.gtype == points) return length;
	for (size_t i=0; i < geom.parts.size(); i++) {
		length += length_line_plane(geom.parts[i].x, geom.parts[i].y);
		for (size_t j=0; j < geom.parts[i].holes.size(); j++) {
			length += length_line_plane(geom.parts[i].holes[j].x, geom.parts[i].holes[j].y);
		}
	}
	return length;
}


std::vector<double> SpatVector::length() {

	size_t s = size();
	std::vector<double> r;
	r.reserve(s);

	double m = srs.to_meter();
	m = std::isnan(m) ? 1 : m;

	if (srs.wkt.empty()) {
		addWarning("unknown CRS. Results can be wrong");
	}

	if (m == 0) {
		struct geod_geodesic g;
		double a = 6378137;
		double f = 1 / 298.257223563;
		geod_init(&g, a, f);
		for (size_t i=0; i<s; i++) {
			r.push_back(length_lonlat(g, geoms[i]));
		}
	} else {
		for (size_t i=0; i<s; i++) {
			r.push_back(length_plane(geoms[i]) * m);
		}
	}
	return r;
}


double edges_geom(const SpatGeom &geom) {
	size_t edges = 0;
	if (geom.gtype == points) return edges;
	for (size_t i=0; i<geom.parts.size(); i++) {
		edges += geom.parts[i].y.size();
		for (size_t j=0; j<geom.parts[i].holes.size(); j++) {
			edges += geom.parts[i].holes[j].y.size()-1;
		}
	}
	return edges-1;
}

std::vector<size_t> SpatVector::nseg() {

	size_t s = size();
	std::vector<size_t> r;
	r.reserve(s);
	for (size_t i=0; i<s; i++) {
		r.push_back(edges_geom(geoms[i]));
	}
	return r;
}


void SpatVector::fix_lonlat_overflow() {

	if (! ((extent.xmin < -180) || (extent.xmax > 180))) { return; }
	SpatExtent world(-180, 180, -90, 90);

	std::string vt = type();
	if (vt == "points") {
		for (size_t i=0; i<geoms.size(); i++) {
			SpatGeom g = geoms[i];
			for (size_t j=0; j<g.parts.size(); j++) {
				for (size_t k=0; k<g.parts[j].x.size(); k++) {
					if (geoms[i].parts[j].x[k] < -180) { geoms[i].parts[j].x[k] += 360; }
					if (geoms[i].parts[j].x[k] > 180) { geoms[i].parts[j].x[k] -= 360; }
				}
			}
		}
	} else {
		SpatExtent east(-360, -180, -180, 180);
		SpatExtent west(180, 360, -180, 180);

		for (size_t i=0; i<geoms.size(); i++) {
			if (geoms[i].extent.xmin < -180) {
				SpatVector v(geoms[i]);
				if (geoms[i].extent.xmax <= -180) {
					v = v.shift(360, 0);
				} else {
					SpatVector add = v.crop(east, false);
					add = add.shift(360, 0);
					v = v.crop(world, false);
					v.geoms[i].addPart(add.geoms[0].parts[0]);
				}
				replaceGeom(v.geoms[0], i);
			}
			if (geoms[i].extent.xmax > 180) {
				SpatVector v(geoms[i]);
				if (geoms[i].extent.xmin >= 180) {
					v = v.shift(-360, 0);
				} else {
					SpatVector add = v.crop(west, false);
					add = add.shift(-360, 0);
					v = v.crop(world, false);
					v.geoms[i].addPart(add.geoms[0].parts[0]);
				}
				replaceGeom(v.geoms[0], i);
			}
		}
	}

	if ((extent.ymax > 90) || (extent.ymin < -90)) {
		SpatVector out = crop(world, false);
		geoms = out.geoms;
		extent = out.extent;
		df = out.df;
		srs = out.srs;
	}
	return;
}
