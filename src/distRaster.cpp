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

#include "spatRaster.h"
#include "distance.h"
#include <limits>
#include <cmath>
#include "ggeodesic.h"
#include "recycle.h"
#include "math_utils.h"


void shortDistPoints(std::vector<double> &d, const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &px, const std::vector<double> &py, const bool& lonlat, const double &lindist) {
	if (lonlat) {
		distanceToNearest_lonlat(d, x, y, px, py);
	} else {
		distanceToNearest_plane(d, x, y, px, py, lindist);
	}
}



SpatRaster SpatRaster::distance(SpatVector p, SpatOptions &opt) {

	SpatRaster out = geometry();
	if (source[0].srs.wkt == "") {
		out.setError("CRS not defined");
		return(out);
	}

	double m = source[0].srs.to_meter();
	m = std::isnan(m) ? 1 : m;

	SpatRaster x;
	std::string gtype = p.type();
	if (gtype != "points") {
		SpatOptions ops;
		std::vector<double> feats(p.size(), 1) ;
		x = rasterize(p, "", feats, NAN, false, false, false, false, false, ops);
		if (gtype == "polygons") {
			std::string etype = "inner";
			x = x.edges(false, etype, 8, 0, ops);
		}
		p = x.as_points(false, true, opt);
	}


	if (p.size() == 0) {
		out.setError("no cells to compute distance from");
		return(out);
	}
	
	bool lonlat = is_lonlat(); // m == 0
	unsigned nc = ncol();
	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

 	if (!out.writeStart(opt)) {
		readStop();
		return out;
	}
	std::vector<std::vector<double>> pxy = p.coordinates();
	std::vector<double> v, cells;
	
	for (size_t i = 0; i < out.bs.n; i++) {
		double s = out.bs.row[i] * nc;
		cells.resize(out.bs.nrows[i] * nc) ;
		std::iota(cells.begin(), cells.end(), s);
		
		if (gtype != "points") {
			v = x.readBlock(out.bs, i);
		} else if (hasValues()) {
			v = readBlock(out.bs, i);
		}
		for (size_t j=0; j<v.size(); j++) {
			if (!std::isnan(v[j])) {
				cells[j] = -1;
			}
		}
		std::vector<std::vector<double>> xy = xyFromCell(cells);
		std::vector<double> d(cells.size(), 0); 
		shortDistPoints(d, xy[0], xy[1], pxy[0], pxy[1], lonlat, m);
		if (!out.writeValues(d, out.bs.row[i], out.bs.nrows[i], 0, nc)) return out;
	}
	out.writeStop();
	readStop();
	return(out);
}


SpatRaster SpatRaster::distance(SpatOptions &opt) {
	SpatRaster out = geometry(1);
	if (!hasValues()) {
		out.setError("SpatRaster has no values");
		return out;
	}

	SpatOptions ops(opt);
	bool warn = false;
	std::string msg;
	if (nlyr() > 1) {
		warn = true;
		msg = "distance computations are only done for the first input layer";
		std::vector<unsigned> lyr = {0};
		*this = subset(lyr, ops);
	}

	out = edges(false, "inner", 8, 0, ops);
	SpatVector p = out.as_points(false, true, opt);
	out = out.distance(p, opt);
	if (warn) {
		out.addWarning(msg);
	}
	return out;
}




std::vector<double> SpatVector::distance(bool sequential) {
	std::vector<double> d;
	if (srs.is_empty()) {
		setError("crs not defined");
		return(d);
	}
	double m = srs.to_meter();
	m = std::isnan(m) ? 1 : m;
	bool lonlat = is_lonlat(); // m == 0
	
//	if ((!lonlat) || (gtype != "points")) {
	std::string gtype = type();
	if (gtype != "points") {
		d = geos_distance(sequential);
		if ((!lonlat) && (m != 1)) {
			for (double &i : d) i *= m;
		}
		return d;
	} else {
		if (sequential) {
			std::vector<std::vector<double>> p = coordinates();
			size_t n = p[0].size();
			d.reserve(n);
			d.push_back(0);
			n -= 1;
			if (lonlat) {
				for (size_t i=0; i<n; i++) {
					d.push_back(
						distance_lonlat(p[0][i], p[1][i], p[0][i+1], p[1][i+1])
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
							distance_lonlat(p[0][i], p[1][i], p[0][j], p[1][j])
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
	
	return d;
}


std::vector<double>  SpatVector::distance(SpatVector x, bool pairwise) {

	std::vector<double> d;

	if (srs.is_empty() || x.srs.is_empty()) {
		setError("SRS not defined");
		return(d);
	}
	if (! srs.is_same(x.srs, false) ) {
		setError("SRS do not match");
		return(d);
	}
	double m = srs.to_meter();
	m = std::isnan(m) ? 1 : m;
	bool lonlat = is_lonlat();

	std::string gtype = type();
	std::string xtype = x.type();
	if ((!lonlat) || (gtype != "points") || (xtype != "points")) {
		d = geos_distance(x, pairwise);
		if ((!lonlat) && (m != 1)) {
			for (double &i : d) i *= m;
		}
		return d;
	}

	size_t s = size();
	size_t sx = x.size();
	if (s != sx) {
		pairwise = false;
	}
	size_t n = pairwise ? s : s*sx;
	d.resize(n);	
	std::vector<std::vector<double>> p = coordinates();
	std::vector<std::vector<double>> px = x.coordinates();

	if (pairwise) {
		if (lonlat) {
			for (size_t i = 0; i < s; i++) {
				d[i] = distance_lonlat(p[0][i], p[1][i], px[0][i], px[1][i]);
			}
		} else { // not reached
			for (size_t i = 0; i < s; i++) {
				d[i] = distance_plane(p[0][i], p[1][i], px[0][i], px[1][i]) * m;
			}		
		} 
	} else {	
		if (lonlat) {
			for (size_t i=0; i<s; i++) {
				size_t k = i * sx;
				for (size_t j=0; j<sx; j++) {
					d[k+j] = distance_lonlat(p[0][i], p[1][i], px[0][j], px[1][j]);
				}
			}
		} else { // not reached
			for (size_t i=0; i<s; i++) {
				size_t k = i * sx;
				for (size_t j=0; j<sx; j++) {
					d[k+j] = distance_plane(p[0][i], p[1][i], px[0][j], px[1][j]) * m;
				}
			}
		} 
	}

	return d;
}




std::vector<double> broom_dist_planar(std::vector<double> &v, std::vector<double> &above, std::vector<double> res, size_t nr, size_t nc, double lindist) {

	double dx = res[0] * lindist;
	double dy = res[1] * lindist;
	double dxy = sqrt(dx * dx + dy *dy);


	std::vector<double> dist(v.size(), 0);

	//top to bottom
    //left to right

	if ( std::isnan(v[0]) ) { //first cell, no cell left of it
		dist[0] = above[0] + dy;
	}
	for (size_t i=1; i<nc; i++) { //first row, no row above it, use "above"
		if (std::isnan(v[i])) {
			dist[i] = std::min(std::min(above[i] + dy, above[i-1] + dxy), dist[i-1] + dx);
		}
	}

	for (size_t r=1; r<nr; r++) { //other rows
		size_t start=r*nc;

		if (std::isnan(v[start])) {
			dist[start] = dist[start-nc] + dy;
		}
	
		size_t end = start+nc;
		for (size_t i=(start+1); i<end; i++) {
			if (std::isnan(v[i])) {
				dist[i] = std::min(std::min(dist[i-1] + dx, dist[i-nc] + dy), dist[i-nc-1] + dxy);
			}
		}
	}
	
	//right to left
	if ( std::isnan(v[nc-1])) { //first cell
		dist[nc-1] = std::min(dist[nc-1], above[nc-1] + dy);
	}
	
	for (int i=(nc-2); i > -1; i--) { // other cells on first row
		if (std::isnan(v[i])) {
			dist[i] = std::min(std::min(std::min(dist[i+1] + dx, above[i+1] + dxy), above[i] + dy), dist[i]);
		}
	}

	for (size_t r=1; r<nr; r++) { // other rows
		size_t start=(r+1)*nc-1;
		if (std::isnan(v[start])) {
			dist[start] = std::min(dist[start], dist[start-nc] + dy);
		}
		for (size_t i=start-1; i>=(r*nc); i--) {
			if (std::isnan(v[i])) {
				dist[i] = std::min(std::min(std::min(dist[i], dist[i+1] + dx), dist[i-nc] + dy), dist[i-nc+1] + dxy);
			}
		}
	}

	size_t off = (nr-1) * nc;
	above = std::vector<double>(dist.begin()+off, dist.end());

	return dist;
}


void DxDxy(const double &lat, const int &row, double xres, double yres, const int &dir, double &dx, double &dxy) {
	double thislat = lat + row * yres * dir;
	xres /= 2;
	yres /= 2;
	dx  = distance_lonlat(-xres, thislat     , xres, thislat);
	dxy = distance_lonlat(-xres, thislat-yres, xres, thislat+yres);
	//double dy = distance_lonlat(0, 0, -yres, yres, a, f);
//	Rcpp::Rcout << thislat << " " << row << " " << xres << " " << yres << std::endl;
//	Rcpp::Rcout << dy << " " << dx << " " << dxy << std::endl;
}


std::vector<double> broom_dist_geo(std::vector<double> &v, std::vector<double> &above, std::vector<double> res, size_t nr, size_t nc, double lat, double latdir) {

	double dy = distance_lonlat(0, 0, 0, res[0]);
	double dx, dxy;
	
	std::vector<double> dist(v.size(), 0);

	//top to bottom
    //left to right
	DxDxy(lat, 0, res[0], res[1], latdir, dx, dxy);
	if ( std::isnan(v[0]) ) { //first cell, no cell left of it
		dist[0] = above[0] + dy;
	}
	for (size_t i=1; i<nc; i++) { //first row, no row above it, use "above"
		if (std::isnan(v[i])) {
			dist[i] = std::min(std::min(above[i] + dy, above[i-1] + dxy), dist[i-1] + dx);				
		}	
	}
	

	for (size_t r=1; r<nr; r++) { //other rows
		DxDxy(lat, r, res[0], res[1], latdir, dx, dxy);
		size_t start=r*nc;

		if (std::isnan(v[start])) {
			dist[start] = dist[start-nc] + dy;
		}
	
		size_t end = start+nc;
		for (size_t i=(start+1); i<end; i++) {
			if (std::isnan(v[i])) {
				dist[i] = std::min(std::min(dist[i-1] + dx, dist[i-nc] + dy), dist[i-nc-1] + dxy);
			}
		}
	}
	
	//right to left
	DxDxy(lat, 0, res[0], res[1], latdir, dx, dxy);
	if ( std::isnan(v[nc-1])) { //first cell
		dist[nc-1] = std::min(dist[nc-1], above[nc-1] + dy);
	}
	
	for (int i=(nc-2); i > -1; i--) { // other cells on first row
		if (std::isnan(v[i])) {
			dist[i] = std::min(std::min(std::min(dist[i+1] + dx, above[i+1] + dxy), above[i] + dy), dist[i]);
		}
	}

	for (size_t r=1; r<nr; r++) { // other rows
		DxDxy(lat, r, res[0], res[1], latdir, dx, dxy);
	
		size_t start=(r+1)*nc-1;
		if (std::isnan(v[start])) {
			dist[start] = std::min(dist[start], dist[start-nc] + dy);
		}
		for (size_t i=start-1; i>=(r*nc); i--) {
			if (std::isnan(v[i])) {
				dist[i] = std::min(std::min(std::min(dist[i], dist[i+1] + dx), dist[i-nc] + dy), dist[i-nc+1] + dxy);
			}
		}
	}

	size_t off = (nr-1) * nc;
	above = std::vector<double>(dist.begin()+off, dist.end());

	return dist;
}

	

SpatRaster SpatRaster::gridDistance(SpatOptions &opt) {

	SpatRaster out = geometry(1);
	SpatOptions ops(opt);
	bool warn = false;
	std::string msg;
	if (nlyr() > 1) {
		warn = true;
		msg = "distance computations are only done for the first input layer";
		std::vector<unsigned> lyr = {0};
		*this = subset(lyr, ops);
	}

	if (!hasValues()) {
		out.setError("cannot compute distance for a raster with no values");
		return out;
	}

	//bool isgeo = out.islonlat

	double m = source[0].srs.to_meter();
	m = std::isnan(m) ? 1 : m;

	std::vector<double> res = resolution();

	SpatRaster first = out.geometry();

	std::string tempfile = "";
	std::vector<double> above(ncol(), std::numeric_limits<double>::infinity());
    std::vector<double> d, v, vv;

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
	std::string filename = opt.get_filename();
	opt.set_filenames({""});
 	if (!first.writeStart(opt)) { return first; }

//	bool lonlat = is_lonlat(); 
	size_t nc = ncol();
	for (size_t i = 0; i < first.bs.n; i++) {
        v = readBlock(first.bs, i);
//		if (lonlat) {			
//			double lat = yFromRow(first.bs.row[i]);
//			d = broom_dist_geo(v, above, res, first.bs.nrows[i], nc, lat, -1);
//		} else {
			d = broom_dist_planar(v, above, res, first.bs.nrows[i], nc, m);
//		}
		if (!first.writeValues(d, first.bs.row[i], first.bs.nrows[i], 0, nc)) return first;
	}
	first.writeStop();

	if (!first.readStart()) {
		out.setError(first.getError());
		return(out);
	}

	opt.set_filenames({filename});
	above = std::vector<double>(ncol(), std::numeric_limits<double>::infinity());

  	if (!out.writeStart(opt)) {
		readStop();
		return out;
	}
	for (int i = out.bs.n; i>0; i--) {
        v = readBlock(out.bs, i-1);
		std::reverse(v.begin(), v.end());
//		if (lonlat) {			
//			double lat = yFromRow(out.bs.row[i-1] + out.bs.nrows[i-1] - 1);
//			d = broom_dist_geo(v, above, res, out.bs.nrows[i-1], nc, lat, 1);
//		} else {
			d = broom_dist_planar(v, above, res, out.bs.nrows[i-1], nc, m);
//		}
		vv = first.readBlock(out.bs, i-1);
	    std::transform (d.rbegin(), d.rend(), vv.begin(), vv.begin(), [](double a, double b) {return std::min(a,b);});
		if (!out.writeValues(vv, out.bs.row[i-1], out.bs.nrows[i-1], 0, nc)) return out;
	}
	out.writeStop();
	readStop();
	first.readStop();
	
	if (warn) out.addWarning(msg);
	return(out);
	
}


/*
std::vector<double> do_edge(std::vector<double> &d, size_t nrow, size_t ncol, bool before, bool after, bool classes, bool inner, unsigned dirs) {

	bool falseval = 0;

	size_t n = nrow * ncol;
	std::vector<double> val(n, NAN);

	// main
	int r[8] = { -1,0,0,1 , -1,-1,1,1};
	int c[8] = { 0,-1,1,0 , -1,1,-1,1};
		// first col
	int fr[5] = {-1,0,1,-1,1};
	int fc[5] = { 0,1,0, 1,1};
		// last col
	int lr[5] = { -1, 0,1, -1, 1};
	int lc[5] = {  0,-1,0, -1,-1};


	// first row
	int br[5] = {  0,0,1 , 1,1};
	int bc[5] = { -1,1,0 ,-1,1};
		// first col
	int bfr[3] = { 0,1 ,1};
	int bfc[3] = { 1,0 ,1};
		// last col
	int blr[3] = {  0,1 , 1};
	int blc[3] = { -1,0 ,-1};


	// last row
	int ar[5] = { -1,0,0, -1,-1};
	int ac[5] = { 0,-1,1, -1, 1};
		// first col
	int afr[3] = { -1,0,-1};
	int afc[3] = { 0 ,1, 1};
		// last col
	int alr[3] = { -1,0, -1};
	int alc[3] = { 0,-1, -1};
	

	size_t rowoff = 0;
	size_t nrows = nrow;
	if (before) {
		rowoff = 1;
	}
	if (after) {
		nrows++;
	}
	size_t hrdirs = dirs == 4 ? 3 : 5;
	size_t hcdirs = dirs == 4 ? 2 : 3;

	if (classes) {  // by class

		for (size_t i = 1; i < (nrows); i++) {
			for (size_t j = 1; j < (ncol-1); j++) {
				size_t cell = i * ncol+j ;
				double test = d[cell + r[0] * ncol + c[0]];
				val[cell] = std::isnan(test) ? NAN : falseval;
				for (size_t k=1; k<dirs; k++) {
					double v = d[cell+r[k]*ncol +c[k]];
					if (std::isnan(test)) {
						if (!std::isnan(v)) {
							val[cell] = 1;
							break;
						}
					} else if (test != v) {
						val[cell] = 1;
						break;
					}
				}
			}
		}		
	} else { // not by class
		if (inner) {  ////// inner //// 	
			if (!before) { // no row above
				for (size_t j = 1; j < (ncol-1); j++) {
					// cell = j
					if (!std::isnan(d[j])) {
						val[j] = 0;
						for (size_t k=0; k < hrdirs; k++) {		
							if ( std::isnan(d[j + br[k] * ncol + bc[k] ])) {
								val[j] = 1;
								break;
							}
						}
					}
				} // first column of first row
				// cell = j = 0
				if (!std::isnan(d[0])) {
					val[0] = 0;
					for (size_t k=0; k < hcdirs; k++) {		
						if ( std::isnan(d[bfr[k] * ncol + bfc[k] ])) {
							val[0] = 1;
							break;
						}
					}
				} // last column of first row
				size_t cell = ncol-1;
				if (!std::isnan(d[cell])) {
					val[cell] = 0;
					for (size_t k=0; k < hcdirs; k++) {		
						if ( std::isnan(d[cell + blr[k] * ncol + blc[k] ])) {
							val[cell] = 1;
							break;
						}
					}
				}
			}
			if (!after) { // no row below
				size_t i = nrows-1;
				for (size_t j = 1; j < (ncol-1); j++) {
					size_t cell = i * ncol + j;
					size_t outcell = (i-rowoff) * ncol + j;
					if (!std::isnan(d[cell])) {
						val[outcell] = 0;
						for (size_t k=0; k < hrdirs; k++) {		
							if ( std::isnan(d[cell+ ar[k] * ncol + ac[k] ])) {
								val[outcell] = 1;
								break;
							}
						}
					}
				} // first cell for last row
				size_t cell = (nrows-1) * ncol;
				size_t outcell = (nrows-1-rowoff) * ncol;
				if (!std::isnasn(d[cell])) {
					val[outcell] = 0;
					for (size_t k=0; k < hcdirs; k++) {		
						if ( std::isnan(d[cell + afr[k] * ncol + afc[k] ])) {
							val[outcell] = 1;
							break;
						}
					}
				} // last cell for last row
				cell += ncol-1;
				outcell += ncol-1;
				if (!std::isnan(d[cell])) {
					val[outcell] = 0;
					for (size_t k=0; k < hcdirs; k++) {		
						if ( std::isnan(d[cell+ alr[k] * ncol + alc[k] ])) {
							val[outcell] = 1;
							break;
						}
					}
				}	
			} // other rows 

			
			for (size_t i = 1; i < nrows; i++) {
				for (size_t j = 1; j < (ncol-1); j++) {
					size_t cell = i * ncol + j;
					if (!std::isnan(d[cell])) {
						size_t outcell = (i-rowoff) * ncol + j;
						val[outcell] = 0;
						for (size_t k=0; k < dirs; k++) {		
							if ( std::isnan(d[cell+ r[k] * ncol + c[k] ])) {
								val[outcell] = 1;
								break;
							}
						}
					}
				}

				// first column
				size_t cell = i * ncol;
				size_t outcell = (i-rowoff) * ncol;
				if (!std::isnan(d[cell])) {
					val[outcell] = 0;
					for (size_t k=0; k < hrdirs; k++) {		
						if ( std::isnan(d[cell + fr[k] * ncol + fc[k] ])) {
							val[outcell] = 1;
							break;
						}
					}
				}
				// last column
				cell += ncol - 1;
				outcell += ncol - 1;
				if (!std::isnan(d[cell])) {
					val[outcell] = 0;
					for (size_t k=0; k < hrdirs; k++) {		
						if ( std::isnan( d[cell + lr[k] * ncol + lc[k] ])) {
							val[outcell] = 1;
							break;
						}
					}
				}
			}

		} else { ////// outer //// 	

	
			if (!before) {
				for (size_t j = 1; j < (ncol-1); j++) {
					size_t cell = j;
					if (std::isnan(d[cell])) {
						val[cell] = NAN;
						for (size_t k=0; k < hcdirs; k++) {		
							if ( !std::isnan(d[j + br[k] * ncol + bc[k] ])) {
								val[cell] = 1;
								break;
							}
						}
					} else {
						val[cell] = 0;
					}
				}
			}
			if (!after) {
				size_t i = (nrow - 1) * ncol;
				for (size_t j = 1; j < (ncol-1); j++) {
					size_t cell = (i-rowoff) * ncol + j;
					if (std::isnan(d[cell])) {
						val[cell] = NAN;
						for (size_t k=0; k < hcdirs; k++) {		
							if (!std::isnan(d[cell+ ar[k] * ncol + ac[k] ])) {
								val[cell] = 1;
								break;
							}
						}
					} else {
						val[cell] = 0;
					}
				}
			}

			for (size_t i = 1; i < nrows; i++) {
				for (size_t j = 1; j < (ncol-1); j++) {
					size_t cell = (i-rowoff) * ncol + j;
					if (std::isnan(d[cell])) {
						val[cell] = NAN;
						for (size_t k=0; k<dirs; k++) {
							if (!std::isnan(d[cell + r[k] * ncol + c[k]])) {
								val[cell] = 1;
								break;
							}
						}
					} else {
						val[cell] = 0;
					}
				}
			}
		}
	
	}
//	val.erase(val.begin(), val.begin()+ncol);
//	val.erase(val.end()-ncol, val.end());	
	return(val);
}

*/

/*
std::vector<double> get_border(std::vector<double> xd, size_t nrows, size_t ncols, bool classes, std::string edgetype, unsigned dirs) {

	size_t n = nrows * ncols;

	std::vector<double> xval(n, 0);
	//Rcpp::Rcout << "hello" << std::endl;

	int r[8] = {-1,0,0,1, -1,-1,1,1};
	int c[8] = {0,-1,1,0, -1,1,-1,1};
	int falseval = 0;

	if (!classes) {
		if (edgetype == "inner") { 
			for (size_t i = 1; i < (nrows-1); i++) {
				for (size_t j = 1; j < (ncols-1); j++) {
					size_t cell = i*ncols+j;
					if (std::isnan(xd[cell])) {
						xval[cell] = NAN;
					} else {
						xval[cell] = falseval;
						for (size_t k=0; k< dirs; k++) {
							if (std::isnan (xd[cell + r[k] * ncols + c[k]])) {
								xval[cell] = 1;
								break;
							}
						}
					}
				}
			}
	
		} else { // if (edgetype == "outer"
			for (size_t i = 1; i < (nrows-1); i++) {
				for (size_t j = 1; j < (ncols-1); j++) {
					size_t cell = i*ncols+j;
					xval[cell] = falseval;
					if (std::isnan(xd[cell])) {
						xval[cell] = NAN;
					} else {
						for (size_t k=0; k < dirs; k++) {		
							if (std::isnan(xd[cell+ r[k] * ncols + c[k] ])) {
								xval[cell] = 1;
								break;
							}
						}
					}
				}
			}
		} 
	} else { // by class
		for (size_t i = 1; i < (nrows-1); i++) {
			for (size_t j = 1; j < (ncols-1); j++) {
				size_t cell = i*ncols+j;
				double test = xd[ cell + r[0]*ncols + c[0] ];
				if (std::isnan(test)) {
					xval[cell] = NAN;
				} else {
					xval[cell] = falseval;
					for (size_t k=1; k < dirs; k++) {
						if (test != xd[ cell+r[k]*ncols +c[k] ]) {
							xval[cell] = 1;
							break;
						}
					}
				}
			}
		}

	}
	return(xval);
}
*/


/*
SpatRaster SpatRaster::edges(bool classes, std::string type, unsigned directions, SpatOptions &opt) {

	SpatRaster out = geometry();
	if (nlyr() > 1) {
		out.setError("boundary detection can only be done for one layer at a time --- to be improved");
		return(out);
	}
	if (!hasValues()) {
		out.setError("SpatRaster has no values");
		return out;
	}

	
	if ((directions != 4) && (directions != 8)) {
		out.setError("directions should be 4 or 8");
		return(out);	
	}
	if ((type != "inner") && (type != "outer")) {
		out.setError("directions should be 'inner' or 'outer'");
		return(out);	
	}
	bool inner = type == "inner";

	size_t nc = ncol();

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

 	if (!out.writeStart(opt)) {
		readStop();
		return out;
	}
	
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> v;
		bool before = false;
		bool after = false;
		if (i == 0) {
			if (out.bs.n == 1) {
				v = readValues(out.bs.row[i], out.bs.nrows[i], 0, nc);
			} else {
				v = readValues(out.bs.row[i], out.bs.nrows[i]+1, 0, nc);
				after = true;
			}	
		} else {
			before = true;
			if (i == out.bs.n) {
				v = readValues(out.bs.row[i]-1, out.bs.nrows[i]+1, 0, nc);
			} else {
				v = readValues(out.bs.row[i]-1, out.bs.nrows[i]+2, 0, nc);
				after = true;
			}
		}
		std::vector<double> vv = do_edge(v, out.bs.nrows[i], nc, before, after, classes, inner, directions);
		if (!out.writeValues(vv, out.bs.row[i], out.bs.nrows[i], 0, nc)) return out;
	}
	out.writeStop();
	readStop();

	return(out);
}
*/

std::vector<double> do_edge(const std::vector<double> &d, const size_t nrow, const size_t ncol, const bool classes, const bool inner, const unsigned dirs, double falseval) {

	size_t n = nrow * ncol;
	std::vector<double> val(n, falseval);
	
	int r[8] = { -1,0,0,1 , -1,-1,1,1};
	int c[8] = { 0,-1,1,0 , -1,1,-1,1};	
	
	if (!classes) {
		if (inner) { // inner
			for (size_t i = 1; i < (nrow-1); i++) {
				for (size_t j = 1; j < (ncol-1); j++) {
					size_t cell = i*ncol+j;
					val[cell] = NAN;
					if ( !std::isnan(d[cell])) {
						val[cell] = falseval;
						for (size_t k=0; k< dirs; k++) {
							if ( std::isnan(d[cell + r[k] * ncol + c[k]])) {
								val[cell] = 1;
								break;
							}
						}
					}
				}
			}
		
		} else { //outer
			for (size_t i = 1; i < (nrow-1); i++) {
				for (size_t j = 1; j < (ncol-1); j++) {
					size_t cell = i*ncol+j;
					val[cell] = falseval;
					if (std::isnan(d[cell])) {
						val[cell] = NAN;
						for (size_t k=0; k < dirs; k++) {			
							if ( !std::isnan(d[cell+ r[k] * ncol + c[k] ])) {
								val[cell] = 1;
								break;
							}
						}
					}
				}
			}
		} 
	} else { // by class
		for (size_t i = 1; i < (nrow-1); i++) {
			for (size_t j = 1; j < (ncol-1); j++) {
				size_t cell = i*ncol+j;
				double test = d[cell+r[0]*ncol+c[0]];
				val[cell] = std::isnan(test) ? NAN : falseval;
				for (size_t k=1; k<dirs; k++) {
					double v = d[cell+r[k]*ncol +c[k]];
					if (std::isnan(test)) {
						if (!std::isnan(v)) {
							val[cell] = 1;
							break;
						}
					} else if (test != v) {
						val[cell] = 1;
						break;
					}
				}
			}
		}

	}
	return(val);
}



void addrowcol(std::vector<double> &v, size_t nr, size_t nc, bool rowbefore, bool rowafter, bool cols) {
	
	if (rowbefore) {
		v.insert(v.begin(), v.begin(), v.begin()+nc);
		nr++;
	}
	if (rowafter) {
		v.insert(v.end(), v.end()-nc, v.end());
		nr++;
	}
	if (cols) {
		for (size_t i=0; i<nr; i++) {
			size_t j = i*(nc+2);
			v.insert(v.begin()+j+nc, v[j+nc-1]);
			v.insert(v.begin()+j, v[j]);
		}
	}
}


void striprowcol(std::vector<double> &v, size_t nr, size_t nc, bool rows, bool cols) {
	if (rows) {
		v.erase(v.begin(), v.begin()+nc);
		v.erase(v.end()-nc, v.end());
		nr -= 2;
	}
	if (cols) {
		nc -= 2;
		for (size_t i=0; i<nr; i++) {
			size_t j = i*nc;
			v.erase(v.begin()+j);
			v.erase(v.begin()+j+nc);
		}
	}
}


SpatRaster SpatRaster::edges(bool classes, std::string type, unsigned directions, double falseval, SpatOptions &opt) {

	SpatRaster out = geometry();
	if (nlyr() > 1) {
		out.addWarning("boundary detection is only done for the first layer");
		std::vector<unsigned> lyr = {0};
		SpatOptions ops(opt);
		*this = subset(lyr, ops);
	}
	if (!hasValues()) {
		out.setError("SpatRaster has no values");
		return out;
	}

	
	if ((directions != 4) && (directions != 8)) {
		out.setError("directions should be 4 or 8");
		return(out);	
	}
	if ((type != "inner") && (type != "outer")) {
		out.setError("directions should be 'inner' or 'outer'");
		return(out);	
	}
	bool inner = type == "inner";

	size_t nc = ncol();
	size_t nr = nrow();
	

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

 	if (!out.writeStart(opt)) {
		readStop();
		return out;
	}
	
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> v;
		//bool before = false;
		//bool after = false;
		if (i == 0) {
			if (out.bs.n == 1) {
				v = readValues(out.bs.row[i], out.bs.nrows[i], 0, nc);
				addrowcol(v, nr, nc, true, true, true);			
			} else {
				v = readValues(out.bs.row[i], out.bs.nrows[i]+1, 0, nc);
				addrowcol(v, nr, nc, true, false, true);			
				//after = true;
			}	
		} else {
			//before = true;
			if (i == out.bs.n) {
				v = readValues(out.bs.row[i]-1, out.bs.nrows[i]+1, 0, nc);
				addrowcol(v, nr, nc, false, true, true);			
			} else {
				v = readValues(out.bs.row[i]-1, out.bs.nrows[i]+2, 0, nc);
				addrowcol(v, nr, nc, false, false, true);			
				//after = true;
			}
		}
		//before, after, 
		std::vector<double> vv = do_edge(v, out.bs.nrows[i]+2, nc+2, classes, inner, directions, falseval);
		striprowcol(vv, out.bs.nrows[i]+2, nc+2, true, true);
		if (!out.writeValues(vv, out.bs.row[i], out.bs.nrows[i], 0, nc)) return out;
	}
	out.writeStop();
	readStop();

	return(out);
}



SpatRaster SpatRaster::buffer(double d, SpatOptions &opt) {
	SpatRaster out = geometry(1);
	if (d <= 0) {
		out.setError("buffer size <= 0; nothing to compute");
		return out;
	}
	
	SpatOptions ops(opt);
	bool warn = false;
	std::string msg;
	if (nlyr() > 1) {
		warn = true;
		msg = "distance computations are only done for the first input layer";
		std::vector<unsigned> lyr = {0};
		*this = subset(lyr, ops);
	}
	
	std::string etype = "inner";
	SpatRaster e = edges(false, etype, 8, 0, ops);
	SpatVector p = e.as_points(false, true, opt);
	out = out.distance(p, ops);
	out = out.arith(d, "<=", false, opt);
	if (warn) addWarning(msg);
	return out;
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
					SpatVector add = v.crop(east);
					add = add.shift(360, 0);
					v = v.crop(world);
					v.geoms[i].addPart(add.geoms[0].parts[0]);
				}
				replaceGeom(v.geoms[0], i);
			}

			if (geoms[i].extent.xmax > 180) {
				SpatVector v(geoms[i]);
				if (geoms[i].extent.xmin >= 180) {
					v = v.shift(-360, 0);
				} else {
					SpatVector add = v.crop(west);
					add = add.shift(-360, 0);
					v = v.crop(world);
					v.geoms[i].addPart(add.geoms[0].parts[0]);
				}
				replaceGeom(v.geoms[0], i);
			}
		}
	}

	if ((extent.ymax > 90) || (extent.ymin < -90)) {	
		SpatVector out = crop(world);
		geoms = out.geoms;
		extent = out.extent;
		df = out.df;
		srs = out.srs;
	}
	return;
}




SpatVector SpatVector::point_buffer(std::vector<double> d, unsigned quadsegs) { 

	SpatVector out;
	std::string vt = type();
	if (vt != "points") {
		out.setError("geometry must be points");
		return out;
	}

	size_t npts = size();

/*
# taken care of by `buffer`	
	for (size_t i=0; i<d.size(); i++) {
		if (d[i] <= 0) {
			d[i] = -d[i];
		}
	}
	recycle(d, npts);
*/
	size_t n = quadsegs * 4;
	double step = 360.0 / n;
	SpatGeom g(polygons);
	g.addPart(SpatPart(0, 0));

	std::vector<std::vector<double>> xy = coordinates();

	if (is_lonlat()) {
		std::vector<double> brng(n);
		for (size_t i=0; i<n; i++) {
			brng[i] = i * step;
		}
		for (size_t i=0; i<npts; i++) {
			if (std::isnan(xy[0][i]) || std::isnan(xy[1][i])) {
				out.addGeom(SpatGeom(polygons));
			} else {
				std::vector<std::vector<double>> dp = destpoint_lonlat(xy[0][i], xy[1][i], brng, d[i], false);
				//close polygons
				dp[0].push_back(dp[0][0]);
				dp[1].push_back(dp[1][0]);
				g.setPart(SpatPart(dp[0], dp[1]), 0);
				out.addGeom(g);
			}
		}
		out.fix_lonlat_overflow();
	} else {
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





double area_polygon_lonlat(geod_geodesic &g, const std::vector<double> &lon, const std::vector<double> &lat) {
	struct geod_polygon p;
	geod_polygon_init(&p, 0);
	size_t n = lat.size();
	for (size_t i=0; i < n; i++) {
		geod_polygon_addpoint(&g, &p, lat[i], lon[i]);
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

	size_t s = size();
	size_t m = mask.size();
	bool domask = false;
	if (m > 0) {
		if (s != mask.size()) {
			setError("mask size is not correct");
			return {NAN};
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
	

	if (srs.wkt == "") {
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
				SpatVector v = project("EPSG:4326");
				if (v.hasError()) {
					setError(v.getError());
					return {NAN};
				}
				return v.area(unit, false, {});
			} else {
				double m = srs.to_meter();
				adj *= std::isnan(m) ? 1 : m * m;	
				for (size_t i=0; i<s; i++) {
					ar.push_back(area_plane(geoms[i]));
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

SpatRaster SpatRaster::rst_area(bool mask, std::string unit, bool transform, SpatOptions &opt) {

	SpatRaster out = geometry(1);
	if (out.source[0].srs.wkt == "") {
		out.setError("empty CRS");
		return out;
	}

	std::vector<std::string> f {"m", "km", "ha"};
	if (std::find(f.begin(), f.end(), unit) == f.end()) {
		out.setError("invalid unit");	
		return out;
	}

	if (opt.names.size() == 0) {
		opt.names = {"area"};
	}
	bool lonlat = is_lonlat();

	SpatOptions mopt = opt;
	if (mask) {
		if (!hasValues()) {
			mask = false;
		} else {
			if (lonlat) {
				opt = SpatOptions(opt);			
			} else {
				if (!readStart()) {
					out.setError(getError());
					return(out);
				}
			}
		}
	}
	
  	if (!out.writeStart(opt)) { return out; }


	if (lonlat) { 
		SpatExtent extent = getExtent();
		SpatExtent e = {extent.xmin, extent.xmin+xres(), extent.ymin, extent.ymax};
		SpatOptions optint(opt);
		SpatRaster onecol = out.crop(e, "near", optint);
		SpatOptions popt(opt);
		SpatVector p = onecol.as_polygons(false, false, false, false, popt);
		if (p.hasError()) {
			out.setError(p.getError());
			return out;
		}
		std::vector<double> a = p.area(unit, true, {});
		size_t nc = ncol();
		for (size_t i = 0; i < out.bs.n; i++) {
			std::vector<double> v;
			for (size_t j=0; j<out.bs.nrows[i]; j++) {
				size_t r = out.bs.row[i] + j;
				v.insert(v.end(), nc, a[r]);
			}
			if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
		}

	} else {
		if (transform) {
			SpatExtent extent = getExtent();
			double dy = yres() / 2;
			SpatOptions popt(opt);
			for (size_t i = 0; i < out.bs.n; i++) {
				double ymax = yFromRow(out.bs.row[i]) + dy;
				double ymin = yFromRow(out.bs.row[i] + out.bs.nrows[i]-1) - dy;
				SpatExtent e = {extent.xmin, extent.xmax, ymin, ymax};
				SpatRaster onechunk = out.crop(e, "near", popt);
				SpatVector p = onechunk.as_polygons(false, false, false, false, popt);
				//std::vector<double> cells(onechunk.ncell());
				//std::iota (cells.begin(), cells.end(), 0);
				//onechunk.setValues(cells);
				//SpatVector p = onechunk.as_polygons(false, true, false, false, popt);
				std::vector<double> v;
				if (mask) {
					std::vector<double> m = readValues(bs.row[i], bs.nrows[i], 0, ncol());
					v = p.area(unit, true, m);
				} else {
					v = p.area(unit, true, {});
				}
				if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
			}
		} else {
			double u = unit == "m" ? 1 : unit == "km" ? 1000000 : 10000;
			double m = out.source[0].srs.to_meter();
			double a = std::isnan(m) ? 1 : m;
			a *= xres() * yres() / u;
			for (size_t i = 0; i < out.bs.n; i++) {
				std::vector<double> v(out.bs.nrows[i]*ncol(), a);
				if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
			}
		}
	} 
	out.writeStop();
	if (mask) {
		if (lonlat) {
			out = out.mask(*this, false, NAN, NAN, opt);
		} else {
			readStop();
		}
	}
	return(out);
}


std::vector<double> SpatRaster::sum_area(std::string unit, bool transform, SpatOptions &opt) {

	if (source[0].srs.wkt == "") {
		setError("empty CRS");
		return {NAN};
	}

	std::vector<std::string> f {"m", "km", "ha"};
	if (std::find(f.begin(), f.end(), unit) == f.end()) {
		setError("invalid unit");	
		return {NAN};
	}

	std::vector<double> out(nlyr(), 0);

	if (transform) { //avoid very large polygon objects
		opt.set_memfrac(std::max(0.1, opt.get_memfrac()/2));
	}
	BlockSize bs = getBlockSize(opt);
	if (!readStart()) {
		std::vector<double> err(nlyr(), -1);
		return(err);
	}

	if (is_lonlat()) {
		SpatRaster x = geometry(1);
		SpatExtent extent = x.getExtent();
		SpatExtent e = {extent.xmin, extent.xmin+xres(), extent.ymin, extent.ymax};
		SpatOptions opt;
		SpatRaster onecol = x.crop(e, "near", opt);
		SpatVector p = onecol.as_polygons(false, false, false, false, opt);
		std::vector<double> ar = p.area(unit, true, {});
		size_t nc = ncol();
		if (!hasValues()) {
			out.resize(1);
			for (size_t i=0; i<ar.size(); i++) {
				out[0] += ar[i] * nc;
			}
		} else {
			for (size_t i=0; i<bs.n; i++) {
				std::vector<double> v = readValues(bs.row[i], bs.nrows[i], 0, ncol());
				size_t blockoff = bs.nrows[i] * nc;
				for (size_t lyr=0; lyr<nlyr(); lyr++) {
					size_t lyroff = lyr * blockoff;
					for (size_t j=0; j<bs.nrows[i]; j++) {
						size_t row = bs.row[i] + j;
						size_t offset = lyroff + row * nc;
						size_t n = offset + nc;
						for (size_t k=offset; k<n; k++) {
							if (!std::isnan(v[k])) out[lyr] += ar[row];
						}
					}
				}
			}
		}
	} else if (transform) {
		SpatRaster x = geometry(1);
		SpatExtent extent = x.getExtent();
		double dy = x.yres() / 2;
		SpatOptions popt(opt);
		if (!hasValues()) {
			out.resize(1);
			for (size_t i=0; i<bs.n; i++) {
				double ymax = x.yFromRow(bs.row[i]) + dy;
				double ymin = x.yFromRow(bs.row[i] + bs.nrows[i]-1) - dy;
				SpatExtent e = {extent.xmin, extent.xmax, ymin, ymax};
				SpatRaster onechunk = x.crop(e, "near", popt);
				SpatVector p = onechunk.as_polygons(false, false, false, false, popt);
				p = p.project("EPSG:4326");
				std::vector<double> v = p.area(unit, true, 	{});
				out[0] += accumulate(v.begin(), v.end(), 0);
			}
		} else {
			for (size_t i=0; i<bs.n; i++) {
				double ymax = x.yFromRow(bs.row[i]) + dy;
				double ymin = x.yFromRow(bs.row[i] + bs.nrows[i]-1) - dy;
				SpatExtent e = {extent.xmin, extent.xmax, ymin, ymax};
				SpatRaster onechunk = x.crop(e, "near", popt);
				SpatVector p = onechunk.as_polygons(false, false, false, false, popt);
				p = p.project("EPSG:4326");
				std::vector<double> par = p.area(unit, true, {});
			
				std::vector<double> v = readValues(bs.row[i], bs.nrows[i], 0, ncol());
				unsigned off = bs.nrows[i] * ncol() ;
				for (size_t lyr=0; lyr<nlyr(); lyr++) {
					unsigned offset = lyr * off;
					unsigned n = offset + off;
					for (size_t j=offset; j<n; j++) {
						if (!std::isnan(v[j])) {
							out[lyr] += par[j - offset];
						}
					}
				}
			}
		
		}
	} else {
		double adj = unit == "m" ? 1 : unit == "km" ? 1000000 : 10000;
		double m = source[0].srs.to_meter();
		m = std::isnan(m) ? 1 : m;
		double ar = xres() * yres() * m * m / adj;
		if (!hasValues()) {
			out.resize(1);
			out[0] = ncell() * ar;
		} else {
			for (size_t i=0; i<bs.n; i++) {
				std::vector<double> v = readValues(bs.row[i], bs.nrows[i], 0, ncol());
				unsigned off = bs.nrows[i] * ncol() ;
				for (size_t lyr=0; lyr<nlyr(); lyr++) {
					unsigned offset = lyr * off;
					unsigned n = offset + off;
					for (size_t j=offset; j<n; j++) {
						if (!std::isnan(v[j])) out[lyr]++;
					}
				}
			}
			for (size_t lyr=0; lyr<nlyr(); lyr++) {
				out[lyr] = out[lyr] * ar;
			}
		}
	}
	readStop();
	return(out);
}



//layer<value-area
std::vector<std::vector<double>> SpatRaster::area_by_value(SpatOptions &opt) {

	double m = source[0].srs.to_meter();
	m = std::isnan(m) ? 1 : m;

	if (m != 0) {
		double ar = xres() * yres() * m * m;
		std::vector<std::vector<double>> f = freq(true, false, 0, opt);
		for (size_t i=0; i<f.size(); i++) {
			size_t fs = f[i].size();
			for (size_t j=fs/2; j<fs; j++) {
				f[i][j] *= ar;
			}
		}
		return f;
	} else {
		// to do
		// combine freq and area to get area by latitudes
	
		std::vector<std::vector<double>> out(nlyr());
		return out;
	}
}



void do_flowdir(std::vector<double> &val, std::vector<double> const &d, size_t nrow, size_t ncol, double dx, double dy, unsigned seed) {

	size_t n = nrow * ncol;
	size_t add = val.size();
	val.resize(add+n, NAN);

	std::vector<double> r = {0, 0, 0, 0, 0, 0, 0, 0};
	std::vector<double> p = {1, 2, 4, 8, 16, 32, 64, 128}; // pow(2, j)
	double dxy = sqrt(dx * dx + dy * dy);

	std::default_random_engine generator(seed);
	std::uniform_int_distribution<> U(0, 1);

	for (size_t i = ncol+1; i < (ncol * (nrow-1)-1); i++) {
		if (!std::isnan(d[i])) {
			r[0] = (d[i] - d[i+1]) / dx;
			r[1] = (d[i] - d[i+1+ncol]) / dxy;
			r[2] = (d[i] - d[i+ncol]) / dy;
			r[3] = (d[i] - d[i-1+ncol]) / dxy;
			r[4] = (d[i] - d[i-1]) / dx;
			r[5] = (d[i] - d[i-1-ncol]) / dxy;
			r[6] = (d[i] - d[i-ncol]) / dy;
			r[7] = (d[i] - d[i+1-ncol]) / dxy;
			// using the lowest neighbor, even if it is higher than the focal cell.
			double dmin = r[0];
			int k = 0;
			for (size_t j=1; j<8; j++) {
				if (r[j] > dmin) {
					dmin = r[j];
					k = j;
				} else if (r[j] == dmin) {
					if (U(generator)) {
						dmin = r[j];
						k = j;
					}
				}
			}
			val[i+add] = p[k];
		}
	}
}


void do_TRI(std::vector<double> &val, std::vector<double> const &d, size_t nrow, size_t ncol) {
	size_t n = nrow * ncol;
	size_t add = val.size();
	val.resize(add+n, NAN);
	for (size_t i = ncol+1; i < (ncol * (nrow-1)-1); i++) {
		val[i+add] = (fabs(d[i-1-ncol]-d[i]) + fabs(d[i-1]-d[i]) + fabs(d[i-1+ncol]-d[i]) +  fabs(d[i-ncol]-d[i]) + fabs(d[i+ncol]-d[i]) +  fabs(d[i+1-ncol]-d[i]) + fabs(d[i+1]-d[i]) +  fabs(d[i+1+ncol]-d[i])) / 8;
	}
}

void do_TPI(std::vector<double> &val, const std::vector<double> &d, size_t nrow, size_t ncol) {
	size_t n = nrow * ncol;
	size_t add = val.size();
	val.resize(add+n, NAN);
	for (size_t i = ncol+1; i < (ncol * (nrow-1)-1); i++) {
		val[i+add] = d[i] - (d[i-1-ncol] + d[i-1] + d[i-1+ncol] + d[i-ncol]
		+ d[i+ncol] + d[i+1-ncol] + d[i+1] + d[i+1+ncol]) / 8;
	}
}



void do_roughness(std::vector<double> &val, const std::vector<double> &d, size_t nrow, size_t ncol) {
	size_t n = nrow * ncol;
	size_t add = val.size();
	val.resize(add+n, NAN);
	int incol = ncol;
	int a[9] = { -1-incol, -1, -1+incol, -incol, 0, incol, 1-incol, 1, 1+incol };
	double min, max, v;
	for (size_t i = ncol+1; i < (ncol * (nrow-1)-1); i++) {
		min = d[i + a[0]];
		max = d[i + a[0]];
		for (size_t j = 1; j < 9; j++) {
			v = d[i + a[j]]; 
			if (v > max) {
				max = v;
			} else if (v < min) {
				min = v;
			}
		}
		val[i+add] = max - min;
	}
}


#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif




void to_degrees(std::vector<double>& x, size_t start) { 
	double adj = 180 / M_PI;
	for (size_t i=start; i<x.size(); i++) {
		x[i] *= adj;
	}
}


void do_slope(std::vector<double> &val, const std::vector<double> &d, unsigned ngb, unsigned nrow, unsigned ncol, double dx, double dy, bool geo, std::vector<double> &gy, bool degrees) {

	size_t n = nrow * ncol;
	size_t add = val.size();
	val.resize(add+n, NAN);

	std::vector<double> ddx;
	if (geo) {
		ddx.resize(nrow);
		for (size_t i=0; i<nrow; i++) {
			ddx[i] = distHaversine(-dx, gy[i], dx, gy[i]) / 2 ;
		}
	} 
	double zy, zx; 
	
	if (ngb == 4) {
		if (geo) {
			int q;
			double xwi[2] = {-1,1};
			double xw[2] = {0,0};
			double yw[2] = {-1,1};

			for (size_t i=0; i<2; i++) {
				yw[i] = yw[i] / (2 * dy);
			}		
			for (size_t i = ncol; i < (ncol * (nrow-1)-1); i++) {
				if (i % ncol == 0) {
					q = i / ncol;
					for (size_t k=0; k<2; k++) {
						xw[k] = xwi[k] / (-2 * ddx[q]);
					}
				}
				zx = d[i-1] * xw[0] + d[i+1] * xw[1];
				zy = d[i-ncol] * yw[0] + d[i+ncol] * yw[1];
				val[i+add] = atan(sqrt( pow(zy, 2) + pow(zx, 2) ));
			}
		} else {
		
			double xw[2] = {-1,1};
			double yw[2] = {-1,1};
			for (size_t i=0; i<2; i++) {
				xw[i] /= -2 * dx;
				yw[i] /=  2 * dy;
			}
			for (size_t i = ncol+1; i < (ncol * (nrow-1)-1); i++) {
				zx = d[i-1] * xw[0] + d[i+1] * xw[1];
				zy = d[i-ncol] * yw[0] + d[i+ncol] * yw[1];
				val[i+add] =  atan( sqrt( pow(zy, 2) + pow(zx, 2)  ));
			}
		}	

		
	} else {
	
		if (geo) {
			int q;
			double xwi[6] = {-1,-2,-1,1,2,1};
			double xw[6] = {0,0,0,0,0,0};
			double yw[6] = {-1,1,-2,2,-1,1};
		
			for (size_t i=0; i<6; i++) {
				yw[i] = yw[i] / (8 * dy);
			}
					
			for (size_t i = ncol; i < (ncol * (nrow-1)-1); i++) {
				if (i % ncol == 0) {
					q = i / ncol;
					for (size_t k=0; k<6; k++) {
						xw[k] = xwi[k] / (8 * ddx[q]);
					}
				}
				zx = d[i-1-ncol] * xw[0] + d[i-1] * xw[1] + d[i-1+ncol] * xw[2]
						+ d[i+1-ncol] * xw[3] + d[i+1] * xw[4] + d[i+1+ncol] * xw[5];
				zy = d[i-1-ncol] * yw[0] + d[i-1+ncol] * yw[1] + d[i-ncol] * yw[2] 
						+ d[i+ncol] * yw[3] + d[i+1-ncol] * yw[4] + d[i+1+ncol] * yw[5];
				val[i+add] = atan(sqrt( pow(zy, 2) + pow(zx, 2)  ));
							
			}
		
		} else {
	
			double xw[6] = {-1,-2,-1,1,2,1};
			double yw[6] = {-1,1,-2,2,-1,1};
			for (size_t i=0; i<6; i++) {
				xw[i] /= -8 * dx;
				yw[i] /= 8 * dy;
			}
			for (size_t i = ncol+1; i < (ncol * (nrow-1)-1); i++) {
				zx = d[i-1-ncol] * xw[0] + d[i-1] * xw[1] + d[i-1+ncol] * xw[2]
						+ d[i+1-ncol] * xw[3] + d[i+1] * xw[4] + d[i+1+ncol] * xw[5];
				zy = d[i-1-ncol] * yw[0] + d[i-1+ncol] * yw[1] + d[i-ncol] * yw[2] 
						+ d[i+ncol] * yw[3] + d[i+1-ncol] * yw[4] + d[i+1+ncol] * yw[5];
				val[i+add] = atan(sqrt( pow(zy, 2) + pow(zx, 2) ));
			}
		}
	} 

	if (degrees) {
		to_degrees(val, add);
	} 
}


double dmod(double x, double n) {
	return(x - n * std::floor(x/n));
}



void do_aspect(std::vector<double> &val, const std::vector<double> &d, unsigned ngb, unsigned nrow, unsigned ncol, double dx, double dy, bool geo, std::vector<double> &gy, bool degrees) {

	size_t n = nrow * ncol;
	std::vector<double> ddx;
	if (geo) {
		ddx.resize(nrow);
		for (size_t i=0; i<nrow; i++) {
			ddx[i] = distHaversine(-dx, gy[i], dx, gy[i]) / 2 ;
		}
	} 
	double zy, zx; 
	size_t add = val.size();
	val.resize(add+n, NAN);
	
	//double const pi2 = M_PI / 2;
	double const twoPI = 2 * M_PI;
	double const halfPI = M_PI / 2;
	
	if (ngb == 4) {
		if (geo) {
			int q;
			double xwi[2] = {-1,1};
			double xw[2] = {0,0};
			double yw[2] = {-1,1};
		
			for (size_t i=0; i<2; i++) {
				yw[i] = yw[i] / (2 * dy);
			}			
			for (size_t i = ncol; i < (ncol * (nrow-1)-1); i++) {
				if (i % ncol == 0) {
					q = i / ncol;
					for (size_t k=0; k<2; k++) {
						xw[k] = xwi[k] / (-2 * ddx[q]);
					}
				}
				zx = d[i-1] * xw[0] + d[i+1] * xw[1];
				zy = d[i-ncol] * yw[0] + d[i+ncol] * yw[1];
				zx = atan2(zy, zx);
				val[i+add] = dmod( halfPI - zx, twoPI);
			}
		} else {
		
			double xw[2] = {-1,1};
			double yw[2] = {-1,1};
			for (size_t i=0; i<2; i++) {
				xw[i] = xw[i] / (-2 * dx);
				yw[i] = yw[i] / (2 * dy);
			}
			for (size_t i = ncol+1; i < (ncol * (nrow-1)-1); i++) {
				zx = d[i-1] * xw[0] + d[i+1] * xw[1];
				zy = d[i-ncol] * yw[0] + d[i+ncol] * yw[1];
				zx = atan2(zy, zx);
				val[i+add] = dmod( halfPI - zx, twoPI);
			}
		} 
	} else {
		
		if (geo) {
			int q;
			double xwi[6] = {-1,-2,-1,1,2,1};
			double xw[6] = {0,0,0,0,0,0};
			double yw[6] = {-1,1,-2,2,-1,1};
			
			for (size_t i=0; i<6; i++) {
				yw[i] = yw[i] / (8 * dy);
			}
						
			for (size_t i = ncol+1; i < (ncol * (nrow-1)-1); i++) {
				if (i % ncol == 0) {
					q = i / ncol;
					for (size_t k=0; k<6; k++) {
						xw[k] = xwi[k] / (-8 * ddx[q]);
					}
				}
				zx = d[i-1-ncol] * xw[0] + d[i-1] * xw[1] + d[i-1+ncol] * xw[2]
						+ d[i+1-ncol] * xw[3] + d[i+1] * xw[4] + d[i+1+ncol] * xw[5];
				zy = d[i-1-ncol] * yw[0] + d[i-1+ncol] * yw[1] + d[i-ncol] * yw[2] 
						+ d[i+ncol] * yw[3] + d[i+1-ncol] * yw[4] + d[i+1+ncol] * yw[5];
				zx = atan2(zy, zx);
				val[i+add] = dmod( halfPI - zx, twoPI);
			}
		
		} else {
	
			double xw[6] = {-1,-2,-1,1,2,1};
			double yw[6] = {-1,1,-2,2,-1,1};
			for (size_t i=0; i<6; i++) {
				xw[i] /= -8 * dx;
				yw[i] /= 8 * dy;
			}
			for (size_t i = ncol+1; i < (ncol * (nrow-1)-1); i++) {
				zx = d[i-1-ncol] * xw[0] + d[i-1] * xw[1] + d[i-1+ncol] * xw[2]
						+ d[i+1-ncol] * xw[3] + d[i+1] * xw[4] + d[i+1+ncol] * xw[5];
				zy = d[i-1-ncol] * yw[0] + d[i-1+ncol] * yw[1] + d[i-ncol] * yw[2] 
						+ d[i+ncol] * yw[3] + d[i+1-ncol] * yw[4] + d[i+1+ncol] * yw[5];
				zx = atan2(zy, zx);
				val[i+add] = dmod( halfPI - zx, twoPI);
			}
			
		}
	}

	if (degrees) {
		to_degrees(val, add);
	}
}


SpatRaster SpatRaster::terrain(std::vector<std::string> v, unsigned neighbors, bool degrees, unsigned seed, SpatOptions &opt) {

//aspect, flowdir, TPI, TRI, slope, roughness
	std::sort(v.begin(), v.end());
	v.erase(std::unique(v.begin(), v.end()), v.end());

	SpatRaster out = geometry(v.size());
	bool slope = false, aspect = false, rough=false, tpi=false, tri=false, flow=false;
	for (size_t i=0; i<v.size(); i++) {
		if (v[i] == "slope") {
			slope=true;
		} else if (v[i] == "aspect") {
			aspect=true;
		} else if (v[i] == "roughness") {
			rough=true;
		} else if (v[i] == "TPI") {
			tpi=true;
		} else if (v[i] == "TRI") {
			tri=true;
		} else if (v[i] == "flowdir") {
			flow=true;
		} else {
			out.setError("unknown option: " + v[i]);
			return out;
		}
	}
	out.setNames(v);
	if (nlyr() > 1) {
		out.setError("terrain needs a single layer object");
		return out;
	}
	if ((neighbors != 4) && (neighbors != 8)) {
		out.setError("neighbors should be 4 or 8");
		return out;	
	}
	bool lonlat = is_lonlat();
	double yr = yres();

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}
  	if (!out.writeStart(opt)) {
		readStop();
		return out;
	}
	std::vector<double> y;
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> d = readValues(out.bs.row[i], out.bs.nrows[i], 0, ncol());
		if (lonlat && (aspect || slope)) {
			std::vector<int_64> rows(out.bs.nrows[i]);
			std::iota(rows.begin(), rows.end(), out.bs.row[i]);
			y = yFromRow(rows);
			yr = distHaversine(0, 0, 0, yr);
		}
		std::vector<double> val;
		if (aspect) {
			do_aspect(val, d, neighbors, out.bs.nrows[i], ncol(), xres(), yr, lonlat, y, degrees);
		}
		if (slope) {
			do_slope(val, d, neighbors, out.bs.nrows[i], ncol(), xres(), yr, lonlat, y, degrees);
		}
		if (rough) {
			do_roughness(val, d, out.bs.nrows[i], ncol());
		}
		if (tpi) {
			do_TPI(val, d, out.bs.nrows[i], ncol());
		}
		if (tri) {
			do_TRI(val, d, out.bs.nrows[i], ncol());
		}
		if (flow) {
			double dx = xres();
			double dy = yres();
			if (lonlat) {
				double yhalf = yFromRow((size_t) nrow()/2);
				dx = distHaversine(0, yhalf, dx, yhalf);
				dy = distHaversine(0, 0, 0, dy);
			}
			do_flowdir(val, d, out.bs.nrows[i], ncol(), dx, dy, seed);
		}
		
		if (!out.writeValues(val, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;
	}
	out.writeStop();
	readStop();
	return out; 
}
