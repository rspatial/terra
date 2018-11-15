// Copyright (c) 2018  Robert J. Hijmans
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



#include <vector>
#include "spatraster.h"
#include <algorithm>

// todo: three dimensional focal

std::vector<double> focal_get(std::vector<double> d, std::vector<unsigned> dim, std::vector<unsigned> ngb, double fillvalue) {
  
  // object
	int nrow = dim[0];
	int ncol = dim[1];
  
  // window
	unsigned wrows = ngb[0];
	unsigned wcols = ngb[1];
	int wr = std::floor(wrows / 2);
	int wc = std::floor(wcols / 2);
  
  //	if ((wrows % 2 == 0) | (wcols % 2 == 0))
  //		error("weights matrix must have uneven sides");
  
	unsigned n = nrow * ncol * wrows * wcols;
  	std::vector<double> val(n, fillvalue);
    
	int f = 0;
	int row, col;
	for (int i = 0; i < nrow; i++) {
		for (int j = 0; j < ncol; j++) {
			for (int r=-wr; r <= wr ; r++) {
			row = i+r;
				for (int c=-wc; c <= wc ; c++) {
					col = j+c;
					if (col >= 0 && col < ncol && row >= 0 && row < nrow) {
						val[f] = d[row*ncol+col];
					}
					f++;
				}
			}
		}
    }	
	return(val);
}



std::vector<double> SpatRaster::focal_values(std::vector<unsigned> w, double fillvalue, unsigned row, unsigned nrows) {

	std::vector<unsigned> dim = {nrow, ncol}; 

	int wr = std::floor(w[0]/2);
	unsigned row2 = std::max(unsigned(0), row-wr);
	unsigned nrows2 = std::min(nrows, nrows+wr);

	readStart();
	std::vector<double> d = readValues(row2, nrows2, 0, ncol);	
	readStop();

	std::vector<double> f = focal_get(d, dim, w, fillvalue);	
	if ((row2 < row) | (nrows2 > nrows)) {
		int start = (row2-row) * ncol;
		int end = f.size() - (nrows2-nrows) * ncol;
		return std::vector<double> (&f[start], &f[end]);
//		f = f[start:end];
	}
	return(f);
}



SpatRaster SpatRaster::focal(std::vector<double> w, double fillvalue, bool narm, unsigned fun, SpatOptions opt) {
    
	bool wmat = false;
	int ww;
	std::vector<unsigned> window;
	unsigned size = w.size();
	if (size > 3) {
		wmat = true;
		ww = w.size();
		unsigned wsize = sqrt(ww);
		window.push_back(wsize);
		window.push_back(wsize);
	} else if (size == 2) {
		ww = w[0] * w[1];
		window.push_back(w[0]);
		window.push_back(w[1]);
	} else {
		w.push_back(w[0]);
		ww = w[0] * w[1];
		window.push_back(w[0]);
		window.push_back(w[1]);
	}
	
	SpatRaster out = geometry();
	if (!source[0].hasValues) { return(out); }
	std::vector<unsigned> dim = {0, ncol};
	
 	out.writeStart(opt);
	readStart();
	std::vector<double> v, f, d;
	
	
	std::vector<double> fv;
	
	for (size_t i = 0; i < out.bs.n; i++) {
		d = readValues(out.bs.row[i], out.bs.nrows[i], 0, ncol);	
		dim[0] = out.bs.nrows[i];
		f = focal_get(d, dim, window, fillvalue);
		v.resize(out.bs.nrows[i] * ncol);
		for (size_t j = 0; j < v.size(); j++) {
			double z = 0;
			int n = 0;
			fv.resize(0);
			for (int k = 0; k < ww; k++) {
				int m = j * ww + k; 
				if (std::isnan(f[m])) {
					if (!narm) {
						z = NAN;
						n = 0;
						break;
					}
				} else {
					if (wmat) {
						z = z + f[m] * w[n];
					} else {
						fv.push_back(f[m]);
					}
					n++;
				}
			}
			if (n > 0) {
				if (!wmat) {
					v[j] = z / n;
				} else {
					if (fv.size() == 0) {
						v[j] = NAN;
					} else  if (fun == 0) { //mean
						v[j] = std::accumulate(fv.begin(), fv.end(), 0.0) / fv.size();
					} else if (fun == 1) { //min
						v[j] = *std::min_element(fv.begin(), fv.end());
					} else if (fun == 2) { //max
						v[j] = *std::max_element(fv.begin(), fv.end());
					} else { // sum
						v[j] = std::accumulate(fv.begin(), fv.end(), 0.0);
					}
				}
			} else {
				v[j] = NAN;	
			}
		}
		out.writeValues(v, out.bs.row[i]);
	} 
	readStop();
	out.writeStop();
	return(out);
}

