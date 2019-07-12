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

#include "spatRaster.h"
#include "distance.h"


std::vector<double> shortDistPoints(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &px, const std::vector<double> &py, const bool& lonlat) {
	std::vector<double> out;
	if (lonlat) {
		double a = 6378137.0;
		double f = 1/298.257223563;
		out = distanceToNearest_lonlat(x, y, px, py, a, f);
	} else {
		out = distanceToNearest_plane(x, y, px, py);
	}
	return out;
}



SpatRaster SpatRaster::distance(SpatVector p, SpatOptions &opt) {

	SpatRaster out = geometry();
	if (crs == "") {
		out.setError("CRS not defined");
		return(out);
	}
	
	std::string gtype = p.type();
	if (gtype != "points") {
		SpatOptions ops;
		SpatRaster x = rasterize(p, NAN, ops);
		if (gtype == "polygons") {
			std::string etype = "inner";
			x = x.edges(false, etype, 8, ops);
		}
		p = x.as_points(false, true);
	}

	bool lonlat = is_lonlat();
	unsigned nc = ncol();
 	if (!out.writeStart(opt)) { return out; }
	readStart();
	for (size_t i = 0; i < out.bs.n; i++) {
		double s = out.bs.row[i] * nc;
		std::vector<double> cells(out.bs.nrows[i] * nc) ;
		std::iota (cells.begin(), cells.end(), s);
		std::vector<std::vector<double>> xy = xyFromCell(cells);
		std::vector<std::vector<double>> pxy = p.coordinates();
		std::vector<double> d = shortDistPoints(xy[0], xy[1], pxy[0], pxy[1], lonlat);
		if (!out.writeValues(d, out.bs.row[i], out.bs.nrows[i], 0, nc)) return out;
	}
	out.writeStop();
	readStop();
	return(out);
}


SpatRaster SpatRaster::distance(SpatOptions &opt) {
	SpatRaster out = geometry();
	SpatOptions ops;
	if (nlyr() > 1) {
		out.addWarning("distance computations can only be done for one layer at a time --- to be improved");
		std::vector<unsigned> lyr = {0};
		subset(lyr, ops);
	}
	std::string etype = "inner";
	SpatRaster e = edges(false, etype, 8, ops);
	SpatVector p = e.as_points(false, true);
	out = out.distance(p, opt);
	return out;
}


SpatRaster SpatRaster::buffer(double d, SpatOptions &opt) {
	SpatRaster out = geometry();
	if (d <= 0) {
		out.setError("buffer size <= 0; nothing to compute");
		return out;
	}
	SpatOptions ops;
	if (nlyr() > 1) {
		out.addWarning("buffer computations can only be done for one layer at a time --- to be improved");
		std::vector<unsigned> lyr = {0};
		subset(lyr, ops);
	}
	std::string etype = "inner";
	SpatRaster e = edges(false, etype, 8, ops);
	SpatVector p = e.as_points(false, true);
	out = out.distance(p, ops);
	out = out.arith(d, "<=", false, opt);
	return out;
}




SpatDataFrame SpatVector::distance() {
	SpatDataFrame out;
	std::string gtype = type();
	if (gtype != "points") {
		out.setError("only inmplemented for points --- to be improved");
		return(out);
	}
	std::string crs = getCRS();
	if (crs == "") {
		out.setError("CRS not defined");
		return(out);
	}
	bool lonlat = is_lonlat();

	size_t s = size();
	size_t n = ((s-1) * s)/2;
	std::vector<double> d(n);
	size_t k = 0;
	std::vector<std::vector<double>> p = coordinates();
	if (lonlat) {
		double a = 6378137.0;
		double f = 1/298.257223563;		
		for (size_t i=0; i<(s-1); i++) {
			for (size_t j=(i+1); j<s; j++) {
				d[k] = distance_lonlat(p[0][i], p[1][i], p[0][j], p[1][j], a, f);
				k++;
			}
		}
	} else {
		for (size_t i=0; i<(s-1); i++) {
			for (size_t j=(i+1); j<s; j++) {
				d[k] = distance_plane(p[0][i], p[1][i], p[0][j], p[1][j]);
				k++;
			}
		}
	}
	out.add_column(d, "distance");
	return out;
}


SpatDataFrame SpatVector::distance(SpatVector x, bool pairwise) {

	SpatDataFrame out;
	std::string gtype = type();
	std::string xtype = x.type();
	if ((gtype != "points") || (xtype != "points")) {
		out.setError("only inmplemented for points --- to be improved");
		return(out);
	}
	std::string crs = getCRS();
	std::string xcrs = x.getCRS();
	if (crs == "") {
		out.setError("CRS not defined");
		return(out);
	}
	if (crs != xcrs) {
		out.setError("CRS do not match");
		return(out);
	}
	bool lonlat = is_lonlat();

	size_t s = size();
	size_t sx = x.size();
	if (s != sx) {
		pairwise = false;
	}
	size_t n = pairwise ? s : s*sx;
	std::vector<double> d(n);		
	std::vector<std::vector<double>> p = coordinates();
	std::vector<std::vector<double>> px = x.coordinates();
	
	if (pairwise) {
		if (lonlat) {
			double a = 6378137.0;
			double f = 1/298.257223563;		
			for (size_t i = 0; i < s; i++) {
				d[i] = distance_lonlat(p[0][i], p[1][i], px[0][i], px[1][i], a, f);
			}
		} else {
			for (size_t i = 0; i < s; i++) {
				d[i] = distance_plane(p[0][i], p[1][i], px[0][i], px[1][i]);
			}			
		} 
	} else {		
		if (lonlat) {
			double a = 6378137.0;
			double f = 1/298.257223563;		
			for (size_t i=0; i<s; i++) {
				size_t k = i * sx;
				for (size_t j=0; j<sx; j++) {
					d[k+j] = distance_lonlat(p[0][i], p[1][i], px[0][j], px[1][j], a, f);
				}
			}
		} else {
			for (size_t i=0; i<s; i++) {
				size_t k = i * sx;
				for (size_t j=0; j<sx; j++) {
					d[k+j] = distance_plane(p[0][i], p[1][i], px[0][j], px[1][j]);
				}
			}
		} 
	}

	out.add_column(d, "distance");
	return out;
}

