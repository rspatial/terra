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
#include "read_rst.h"


bool SpatRaster::readStart() {
// for now assuming a single source
// will have to become a loop over sources
	for (size_t i=0; i<nsrc(); i++) {
		if (!source[i].memory) {
		// open filestream
		}
		source[i].open_read = true;
	}
	return true;
}

bool SpatRaster::readStop() {
	for (size_t i=0; i<nsrc(); i++) {
		if (!source[0].memory) {
		// close filestream
		}
		source[i].open_read = false;
	}
	return true;
}


std::vector<double> SpatRaster::readBlock(BlockSize bs, unsigned i){
	std::vector<double> x = readValues(bs.row[i], bs.nrows[i], 0, ncol());
	return(x);
}

void bil_to_bsq(std::vector<double> &v, unsigned nrows, unsigned ncols, unsigned nlyrs) {
	std::vector<std::vector<double>> x(nlyrs);
	for (size_t r=0; r<nrows; r++) {
		unsigned off = r*nlyrs;
		for (size_t i=0; i<nlyrs; i++) {
			unsigned start = (off+i)*ncols;
			x[i].insert(x[i].end(), v.begin()+start, v.begin()+start+ncols);
		}	
	}
	v.resize(0);
	for (size_t i=0; i<nlyrs; i++) {
		v.insert(v.end(), x[i].begin(), x[i].end());
	}
}


std::vector<double> SpatRaster::readValues(unsigned row, unsigned nrows, unsigned col, unsigned ncols){

	unsigned nl = nlyr();
	row = std::min(std::max(unsigned(0), row), nrow()-1);
	col = std::min(std::max(unsigned(0), col), ncol()-1);
	nrows = std::max(unsigned(1), std::min(nrows, nrow()-row));
	ncols = std::max(unsigned(1), std::min(ncols, ncol()-col));
	std::vector<double> out;
	unsigned n = nsrc();

	for (size_t src=0; src<n; src++) {
		if (source[src].memory) {
			if (row==0 && nrows==nrow() && col==0 && ncols==ncol()) {
				out.insert(out.end(), source[src].values.begin(), source[src].values.end());
			} else {
				unsigned a, b;
				unsigned ncells = ncell();
				if (col==0 && ncols==ncol()) {
					for (size_t lyr=0; lyr < nl; lyr++) {
						unsigned add = ncells * lyr;
						a = add + row * ncol();
						b = a + nrows * ncol();
						out.insert(out.end(), source[src].values.begin()+a, source[src].values.begin()+b);
					}
				} else {
					unsigned endrow = row + nrows;
					unsigned endcol = col + ncols;
					for (size_t lyr=0; lyr < nl; lyr++) {
						unsigned add = ncells * lyr;
						for (size_t r = row; r < endrow; r++) {
							a = add + r * ncol();
							out.insert(out.end(), source[src].values.begin()+a+col, source[src].values.begin()+a+endcol);
						}
					}
				}
			}
		} else {
			// read from file
			if (source[src].driver == "raster") {
				std::string file = source[src].filename;
				if (source[src].datatype == "FLT8S") {
					std::vector<double> fvals = readFLT8(file, "BIL", 0, ncell());
					bil_to_bsq(fvals, ncols, nrows, nl);
					out.insert(out.end(), fvals.begin(), fvals.end());
				} else {
					std::vector<double> fvals = readFLT4(file, "BIL", 0, ncell());
					bil_to_bsq(fvals, ncols, nrows, nl);
					out.insert(out.end(), fvals.begin(), fvals.end());
				}

			} else {
				#ifdef useGDAL
				std::vector<double> fvals = readValuesGDAL(src, row, nrows, col, ncols);
				out.insert(out.end(), fvals.begin(), fvals.end());
				#endif // useGDAL
			}
		}
	}
	return out;
}




std::vector<double>  SpatRaster::getValues() {
	std::vector<double> out;
	unsigned n = nsrc();
	for (size_t src=0; src<n; src++) {
		if (source[src].memory) {
			out.insert(out.end(), source[src].values.begin(), source[src].values.end());
		} else if (source[0].driver == "raster") {
			std::vector<double> fvals = readValues(0, nrow(), 0, ncol());
			out.insert(out.end(), fvals.begin(), fvals.end());
		} else {
			#ifdef useGDAL
			std::vector<double> fvals = readValuesGDAL(src, 0, nrow(), 0, ncol());
			out.insert(out.end(), fvals.begin(), fvals.end());
			#endif // useGDAL
		}
	}
	return out;
}

