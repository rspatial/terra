// Robert Hijmans, November 2011 


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

std::vector<double> get_border(std::vector<double> xd, std::vector<unsigned> dim, bool classes, std::string edgetype, unsigned dirs) {

	//R_len_t i, j;
	//SEXP val;
	//int nrows, ncols, n, cell, k;
	//int *xd, *xval;

	unsigned nrows = dim[0];
	unsigned ncols = dim[1];
	unsigned n = nrows * ncols;
	
	std::vector<double> xval(n, NAN);

	int r[8] = {-1,0,0,1, -1,-1,1,1};
	int c[8] = {0,-1,1,0, -1,1,-1,1};	
	int falseval = 0;
	
	if (!classes) {
		if (edgetype == "inner") { 
			for (size_t i = 1; i < (nrows-1); i++) {
				for (size_t j = 1; j < (ncols-1); j++) {
					size_t cell = i*ncols+j;
					if (!std::isnan(xd[cell])) {
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
				double test = xd[ cell+r[0]*ncols+c[0] ];
				if (!std::isnan(test)) {
					xval[cell] = falseval;
				}
				for (size_t k=1; k < dirs; k++) {
					if (test != xd[ cell+r[k]*ncols +c[k] ]) {
						xval[cell] = 1;
						break;
					}
				}
			}
		}

	}
	return(xval);
}




SpatRaster SpatRaster::edges(bool classes, std::string type, unsigned directions, SpatOptions &opt) {

	SpatRaster out = geometry();
	std::vector<unsigned> dim = {nrow(), ncol()}; 

 	if (!out.writeStart(opt)) { return out; }
	readStart();
	std::vector<double> v;
	unsigned nc = ncol();
	for (size_t i = 0; i < out.bs.n; i++) {
		v = readValues(out.bs.row[i], out.bs.nrows[i], 0, nc);
		std::vector<double> vv = get_border(v, dim, classes, type, directions);
		if (!out.writeValues(vv, out.bs.row[i])) return out;
	}
	out.writeStop();
	readStop();

	return(out);
}


