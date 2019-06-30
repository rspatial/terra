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

#include <vector>
#include "spatRaster.h"
#include <algorithm>


// todo: three dimensional focal

std::vector<double> focal_get(std::vector<double> d, std::vector<unsigned> dim, std::vector<unsigned> ngb, double fillvalue, unsigned offset) {

  // object
	unsigned ncols = dim[1];
	unsigned nrows = dim[0];


  // window
	unsigned wrows = ngb[0];
	unsigned wcols = ngb[1];
	unsigned wr = std::floor(wrows / 2);
	unsigned wc = std::floor(wcols / 2);

  //	if ((wrows % 2 == 0) | (wcols % 2 == 0))
  //		error("weights matrix must have uneven sides");

	unsigned n = nrows * ncols * wrows * wcols;
  	std::vector<double> val(n, fillvalue);

	unsigned f = 0;
	for (size_t i = offset; i < nrows; i++) {
		for (size_t j = 0; j < ncols; j++) {
			for (size_t r=-wr; r <= wr ; r++) {
				unsigned row = i+r;
				for (size_t c=-wc; c <= wc ; c++) {
					unsigned col = j+c;
					if (col >= 0 && col < ncols && row >= 0 && row < nrows) {
						val[f] = d[row*ncols+col];
					}
					f++;
				}
			}
		}
    }
	return(val);
}



std::vector<double> SpatRaster::focal_values(std::vector<unsigned> w, double fillvalue, unsigned row, unsigned nrows) {

	unsigned wr = std::floor(w[0]/2);
	unsigned row2 = std::max(unsigned(0), row-wr);
	unsigned nrows2 = std::min(nrow(), nrows+wr);
	unsigned offset = row - row2;
	readStart();
	std::vector<double> d = readValues(row2, nrows2, 0, ncol());
	readStop();

	std::vector<unsigned> dim = {nrows, ncol()};
	std::vector<double> f = focal_get(d, dim, w, fillvalue, offset);

//	if ((row2 < row) | (nrows2 > nrows)) {
//		unsigned start = (row-row2) * dim[1] * w[0] * w[1];
//		unsigned end = start + nrows * dim[1] * w[0] * w[1];
//		return std::vector<double> (f.begin()+start, f.begin()+end);
//	}
	return(f);
}



SpatRaster SpatRaster::focal(std::vector<double> w, double fillvalue, bool narm, std::string fun, SpatOptions &opt) {

	bool wmat = false;
	int ww;
	std::vector<unsigned> window;
	unsigned size = w.size();
	if (size > 3) {
		// this does not look correct (only for square matrices)
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
	std::vector<unsigned> dim = {0, ncol()};

 	if (!out.writeStart(opt)) { return out; }
	readStart();
	std::vector<double> v, f, d;
	std::vector<double> fv;
	
	unsigned wr = std::floor(w[0]/2);

	for (size_t i = 0; i < out.bs.n; i++) {

		unsigned startrow = std::max(unsigned(0), out.bs.row[i]-wr);
		unsigned nrows = out.bs.nrows[i] + wr;
		nrows = (nrows + startrow) > nrow() ? nrow() - startrow : nrows;
		unsigned offset = out.bs.row[i]-startrow;

		d = readValues(startrow, nrows, 0, ncol());
		dim[0] = out.bs.nrows[i];
		f = focal_get(d, dim, window, fillvalue, offset);
		v.resize(out.bs.nrows[i] * ncol());
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
					} else if (fun == "mean") { //mean
						v[j] = std::accumulate(fv.begin(), fv.end(), 0.0) / fv.size();
					} else if (fun == "min") { //min
						v[j] = *std::min_element(fv.begin(), fv.end());
					} else if (fun == "max") { //max
						v[j] = *std::max_element(fv.begin(), fv.end());
					} else { // sum
						v[j] = std::accumulate(fv.begin(), fv.end(), 0.0);
					}
				}
			} else {
				v[j] = NAN;
			}
		}
		if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;

	}
	readStop();
	out.writeStop();
	return(out);
}

