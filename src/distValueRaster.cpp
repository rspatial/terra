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

#include "spatRaster.h"
#include "distance.h"
#include "geodesic.h"
#include "sort.h"
#include "geosphere.h"


void dist_bounds_values(const std::vector<double>& vx, const std::vector<double>& vy, const std::vector<double>& v, const std::vector<double>& rx, const double& ry, size_t& first, size_t& last, const bool& lonlat, const std::string &method, std::vector<double> &d, std::vector<double> &dv) {

	d = std::vector<double>(rx.size(), std::numeric_limits<double>::max());
	dv = std::vector<double>(rx.size(), NAN);
	size_t oldfirst = first;
	first = vx.size();
	last = 0;


	if (lonlat) {
		std::function<double(double, double, double, double)> dfun;
		if (method == "haversine") {
			dfun = distance_hav;
		} else if (method == "cosine") {
			dfun = distance_cos;			
		} else {
			dfun = distance_geo;
		}

		
		for (size_t i=0; i<rx.size(); i++) {
			size_t thisone = 0;
			for (size_t j=oldfirst; j<vx.size(); j++) {
				double dd = dfun(rx[i], ry, vx[j], vy[j]);
				if (dd < d[i]) {
					d[i] = dd;
					dv[i] = v[j];
					thisone = j;
				}
			}
			first = std::min(thisone, first);
			last  = std::max(thisone, last);
		}
	} else {
		for (size_t i=0; i<rx.size(); i++) {
			size_t thisone = 0;
			for (size_t j=oldfirst; j<vx.size(); j++) {
				double dd = distance_plane(rx[i], ry, vx[j], vy[j]);
				if (dd < d[i]) {
					d[i] = dd;
					dv[i] = v[j];
					thisone = j;
				}
			}
			first = std::min(thisone, first);
			last  = std::max(thisone, last);
		}
	}
	last += 1;
}


void dist_only_vals(std::vector<double> &d, std::vector<double> &dv, const std::vector<double>& vx, const std::vector<double>& vy, const std::vector<double>& vv, const std::vector<double>& rx, const std::vector<double>& ry, const size_t& first, const size_t& last, const bool& lonlat, const std::vector<double>& dlast, const std::vector<double>& dvlast, bool skip, const std::vector<double>& v, const std::string& method, bool setNA) {

	size_t rxs = rx.size();
	d.reserve(rxs + dlast.size());
	dv.reserve(rxs + dlast.size());

	double inf = std::numeric_limits<double>::infinity();

	if (lonlat) {

		if (method == "geo") {

			double dd, azi1, azi2;
			struct geod_geodesic g;
			// get a and f from crs?
			double a = 6378137.0;
			double f = 1/298.257223563;
			geod_init(&g, a, f);

			if (skip) {
				for (size_t i=0; i<rxs; i++) {
					if (std::isnan(v[i])) {
						d.push_back(inf);
						dv.push_back(NAN);
						for (size_t j=first; j<last; j++) {
							geod_inverse(&g, ry[i], rx[i], vy[j], vx[j], &dd, &azi1, &azi2);
							if (dd < d[i]) {
								d[i] = dd;
								dv[i] = vv[j];
							}
						}
					} else {
						d.push_back(0);
						dv.push_back(NAN);
					}
				}
			} else { // no skip

				for (size_t i=0; i<rxs; i++) {
					d.push_back(inf);
					dv.push_back(NAN);
					for (size_t j=first; j<last; j++) {
						geod_inverse(&g, ry[i], rx[i], vy[j], vx[j], &dd, &azi1, &azi2);
						if (dd < d[i]) {
							d[i] = dd;
							dv[i] = vv[j];
						}
					}
				} 
			}

		} else { // not geo
	
			std::function<double(double, double, double, double)> dfun;
			if (method == "haversine") {
				dfun = distance_hav;
			} else if (method == "cosine") {
				dfun = distance_cos;
			} 
			if (skip) {
				
				for (size_t i=0; i<rxs; i++) {
					if (std::isnan(v[i])) {
						d.push_back(inf);
						dv.push_back(NAN);
						for (size_t j=first; j<last; j++) {
							double dd = dfun(rx[i], ry[i], vx[j], vy[j]);
							if (dd < d[i]) {
								d[i] = dd;
								dv[i] = vv[j];
							}
						}
					} else {
						d.push_back(0);
						dv.push_back(NAN);
					}
				}		
			} else {	
				for (size_t i=0; i<rxs; i++) {
					d.push_back(inf);
					dv.push_back(NAN);
					for (size_t j=first; j<last; j++) {
						double dd = dfun(rx[i], ry[i], vx[j], vy[j]);
						if (dd < d[i]) {
							d[i] = dd;
							dv[i] = vv[j];
						}
					}
				}
			} 
		}
		
	} else { // not lonlat
		if (skip) {
			for (size_t i=0; i<rxs; i++) {
				if (std::isnan(v[i])) {
					d.push_back(inf);
					dv.push_back(NAN);
					for (size_t j=first; j<last; j++) {
						double dd = distance_plane(rx[i], ry[i], vx[j], vy[j]);
						if (dd < d[i]) {
							d[i] = dd;
							dv[i] = vv[j];
						}
					}
				} else {
					d.push_back(0);
					dv.push_back(NAN);
				}
			}
		} else {
			for (size_t i=0; i<rxs; i++) {
				d.push_back(inf);
				dv.push_back(NAN);
				for (size_t j=first; j<last; j++) {
					double dd = distance_plane(rx[i], ry[i], vx[j], vy[j]);
					if (dd < d[i]) {
						d[i] = dd;
						dv[i] = vv[j];
					}
				}
			}
		}
	}

	d.insert(d.end(), dlast.begin(), dlast.end());
	dv.insert(dv.end(), dvlast.begin(), dvlast.end());


	if (skip) {
		for (size_t i=rxs; i< v.size(); i++) {
			if (!std::isnan(v[i])) {
				d[i] = 0;
				dv[i] = v[i];
			}
		}
		if (setNA) {
			double mxval = std::numeric_limits<double>::max();
			for (size_t i=0; i< v.size(); i++) {
				if (v[i] == mxval) {
					d[i] = NAN;
					dv[i] = NAN;
				}
			}
		}
	}
}


SpatRaster SpatRaster::distance_crds_vals(std::vector<double>& x, std::vector<double>& y, const std::vector<double>& v, const std::string& method, bool skip, bool setNA, std::string unit, double maxdist, SpatOptions &opt) {

	SpatRaster out = geometry();
	if (x.empty()) {
		out.setError("no locations to compute distance from");
		return(out);
	}
	const double toRad = 0.0174532925199433;
	std::vector<std::size_t> pm = sort_order_d(y);
	permute(x, pm);
	permute(y, pm);

	bool lonlat = is_lonlat(); 
	double m=1;
	if (!source[0].srs.m_dist(m, lonlat, unit)) {
		out.setError("invalid unit");
		return(out);
	}


	unsigned nc = ncol();
	if (nrow() > 1000) {
		opt.steps = std::max(opt.steps, (size_t) 4);
		opt.progress = opt.progress * 1.5;
	}
	
 	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}
	std::vector<double> cells;
	std::vector<double> dlast, dvlast;

	std::vector<int_64> cols;
	cols.resize(ncol());
	std::iota(cols.begin(), cols.end(), 0);
	std::vector<double> tox = xFromCol(cols);

	if (lonlat && (method != "geo")) {
		for (double &d : x) d *= toRad;
		for (double &d : y) d *= toRad;
		for (double &d : tox) d *= toRad;
	}

	double oldfirst = 0;
	size_t first = 0;
	size_t last  = x.size();

	std::vector<double> rv;
	if (skip) {
		if (!readStart()) {
			out.setError(getError());
			return(out);
		}
		for (size_t i = 0; i < out.bs.n; i++) {
			cells.resize((out.bs.nrows[i] -1) * nc) ;
			std::iota(cells.begin(), cells.end(), out.bs.row[i] * nc);
			std::vector<std::vector<double>> rxy = xyFromCell(cells);
			double toy = yFromRow(out.bs.row[i] + out.bs.nrows[i] - 1);
			if (lonlat && (method != "geo")) {
				toy *= toRad;
				for (double &d : rxy[0]) d *= toRad;
				for (double &d : rxy[1]) d *= toRad;
			}
			readBlock(rv, out.bs, i);
			dist_bounds_values(x, y, v, tox, toy, first, last, lonlat, method, dlast, dvlast);
			std::vector<double> d, dv;
			dist_only_vals(d, dv, x, y, v, rxy[0], rxy[1], oldfirst, last, lonlat, dlast, dvlast, true, rv, method, setNA);
			oldfirst = first;
			if (maxdist > 0) {
				if (m != 1) {
					for (size_t j=0; j<dv.size(); j++) {
						if ((d[j] / m) > maxdist) dv[j] = NAN;
					}
				} else {
					for (size_t j=0; j<dv.size(); j++) {
						if (d[j] > maxdist) dv[j] = NAN;
					}
				}
			}	
			if (!out.writeBlock(dv, i)) return out;
		}
		readStop();
	} else {
		for (size_t i = 0; i < out.bs.n; i++) {
			double toy = yFromRow(out.bs.row[i] + out.bs.nrows[i] - 1);
			cells.resize((out.bs.nrows[i] -1) * nc) ;
			std::iota(cells.begin(), cells.end(), out.bs.row[i] * nc);
			std::vector<std::vector<double>> rxy = xyFromCell(cells);
			if (lonlat && (method != "geo")) {
				toy *= toRad;
				for (double &d : rxy[0]) d *= toRad;
				for (double &d : rxy[1]) d *= toRad;
			}
			dist_bounds_values(x, y, v, tox, toy, first, last, lonlat, method, dlast, dvlast);
			std::vector<double> d, dv;
			dist_only_vals(d, dv, x, y, v, rxy[0], rxy[1], oldfirst, last, lonlat, dlast, dvlast, false, rv, method, setNA);
			oldfirst = first;
			if (maxdist > 0) {
				if (m != 1) {
					for (size_t j=0; j<dv.size(); j++) {
						if ((d[j] / m) > maxdist) dv[j] = NAN;
					}
				} else {
					for (size_t j=0; j<dv.size(); j++) {
						if (d[j] > maxdist) dv[j] = NAN;
					}
				}
			}	
			if (!out.writeBlock(dv, i)) return out;
		}
	}
	out.writeStop();
	return(out);
}


/*

SpatRaster SpatRaster::distanceValues(double target, double exclude, bool keepNA, bool remove_zero, const std::string method, SpatOptions &opt) {

	SpatRaster out = geometry(1);
	if (!hasValues()) {
		out.setError("SpatRaster has no values");
		return out;
	}

	std::vector<double> v;
	SpatOptions ops(opt);
	size_t nl = nlyr();
	if (nl > 1) {
		std::vector<std::string> nms = getNames();
		if (ops.names.size() == nms.size()) {
			nms = opt.names;
		}
		out.source.resize(nl);
		for (size_t i=0; i<nl; i++) {
			std::vector<size_t> lyr = {i};
			SpatRaster r = subset(lyr, ops);
			ops.names = {nms[i]};
			r = r.distanceValues(target, exclude, keepNA, remove_zero, method, ops);
			out.source[i] = r.source[0];
		}
		if (!opt.get_filename().empty()) {
			out = out.writeRaster(opt);
		}
		return out;
	}

	bool setNA = false;
	std::vector<std::vector<double>> p;
	if (!std::isnan(exclude)) {
		SpatRaster x;
		if (std::isnan(target)) {
			x = replaceValues({exclude}, {target}, 1, false, NAN, false, ops);
			x = x.edges(false, "inner", 8, 1, ops);
			p = x.as_points_value(1, ops);
			if (p.empty()) {
				return out.init({0}, opt);
			}
			std::vector<std::vector<double>> vv = extractXY(p[0], p[1], "", false);
			return distance_crds_vals(p[0], p[1], vv[0], method, true, setNA, unit, opt);

		} else {
			x = replaceValues({exclude, target}, {NAN, NAN}, 1, false, NAN, false, ops);
			x = x.edges(false, "inner", 8, 1, ops);
			p = x.as_points_value(1, ops);
			out = replaceValues({NAN, exclude, target}, {target, NAN, NAN}, 1, false, NAN, false, ops);
		}
	} else if (!std::isnan(target)) {
		SpatRaster x = replaceValues({target}, {NAN}, 1, false, NAN, false, ops);
		x = x.edges(false, "inner", 8, 0, ops);
		p = x.as_points_value(1, ops);
		out = replaceValues({NAN, target}, {std::numeric_limits<double>::max(), NAN}, 1, false, NAN, false, ops);
		setNA = true;
	} else {
		out = edges(false, "inner", 8, 0, ops);
		p = out.as_points_value(1, ops);
		std::vector<std::vector<double>> vv = extractXY(p[0], p[1], "", false);
		return distance_crds_vals(p[0], p[1], vv[0], method, true, setNA, unit, opt);
		
	}
	if (p.empty()) {
		return out.init({0}, opt);
	}
	std::vector<std::vector<double>> vv = extractXY(p[0], p[1], "", false);
	return out.distance_crds_vals(p[0], p[1], vv[0], method, true, setNA, unit, opt);

}

*/
