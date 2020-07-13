// Copyright (c) 2018-2020  Robert J. Hijmans
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

bool SpatRaster::readStart() {
	for (size_t i=0; i<nsrc(); i++) {
		if (source[i].open_read) {
			addWarning("source already open for reading");
			continue;
		}
		if (source[i].memory) {
			source[i].open_read = true;
		} else {
			if (!readStartGDAL(i)) {
				return false;
			}
		}	
	}
	return true;
}

bool SpatRaster::readStop() {
	for (size_t i=0; i<nsrc(); i++) {
		if (source[i].open_read) {
			if (source[i].memory) {
				source[i].open_read = false;
			} else {
				readStopGDAL(i); 
			}	
		}
	}
	return true;
}


// BSQ
std::vector<double> SpatRaster::readBlock(BlockSize bs, unsigned i){
	return readValues(bs.row[i], bs.nrows[i], 0, ncol());
}


// 2D BSQ
std::vector<std::vector<double>> SpatRaster::readBlock2(BlockSize bs, unsigned i) {
	std::vector<double> x = readValues(bs.row[i], bs.nrows[i], 0, ncol());
	std::vector<std::vector<double>> v(nlyr());
	size_t off = bs.nrows[i] * ncol();
	for (size_t i=0; i<nlyr(); i++) {
		v[i] = std::vector<double>(x.begin()+(i*off), x.begin()+((i+1)*off));
	}	
	return(v);
}

// BIP
std::vector<double> SpatRaster::readBlockIP(BlockSize bs, unsigned i) {
	std::vector<double> x = readValues(bs.row[i], bs.nrows[i], 0, ncol());
	std::vector<double> v(x.size());
	size_t off = bs.nrows[i] * ncol();
	size_t nl = nlyr();
	for (size_t i=0; i<nl; i++) {
		std::vector<double> lyr = std::vector<double>(x.begin()+(i*off), x.begin()+((i+1)*off));
		for (size_t j=0; j<off; j++){
			size_t jj = j * nl + i;
			v[jj] = lyr[j];
		}
	}	
	return(v);
}



std::vector<double> SpatRaster::readValues(uint_64 row, uint_64 nrows, uint_64 col, uint_64 ncols){

	std::vector<double> out;
	if (!hasValues()) return out; // or NAs?
	
	row = std::min(std::max(uint_64(0), row), nrow()-1);
	col = std::min(std::max(uint_64(0), col), ncol()-1);
	nrows = std::max(uint_64(1), std::min(nrows, nrow()-row));
	ncols = std::max(uint_64(1), std::min(ncols, ncol()-col));
	if ((nrows==0) | (ncols==0)) {
		return out;
	}
	unsigned n = nsrc();
	
	for (size_t src=0; src<n; src++) {
		unsigned nl = source[src].nlyr;
		if (source[src].memory) {
			if (row==0 && nrows==nrow() && col==0 && ncols==ncol()) {
				out.insert(out.end(), source[src].values.begin(), source[src].values.end());
			} else {
				unsigned ncells = ncell();
				if (col==0 && ncols==ncol()) {
					for (size_t lyr=0; lyr < nl; lyr++) {
						unsigned add = ncells * lyr;
						unsigned a = add + row * ncol();
						unsigned b = a + nrows * ncol();
						out.insert(out.end(), source[src].values.begin()+a, source[src].values.begin()+b);
					}
				} else {
					unsigned endrow = row + nrows;
					unsigned endcol = col + ncols;
					for (size_t lyr=0; lyr < nl; lyr++) {
						unsigned add = ncells * lyr;
						for (size_t r = row; r < endrow; r++) {
							unsigned a = add + r * ncol();
							out.insert(out.end(), source[src].values.begin()+a+col, source[src].values.begin()+a+endcol);
						}
					}
				}
			}
		} else {
			// read from file
			#ifdef useGDAL
			std::vector<double> fvals = readChunkGDAL(src, row, nrows, col, ncols);
			out.insert(out.end(), fvals.begin(), fvals.end());			
			#endif // useGDAL
		}
	}
	return out;
}


std::vector<double> SpatRaster::getValues(long lyr) {
	std::vector<double> out;
	if (lyr < 0) { // default; read all
		unsigned n = nsrc();
		for (size_t src=0; src<n; src++) {
			if (source[src].memory) {
				out.insert(out.end(), source[src].values.begin(), source[src].values.end());
			} else {
				#ifdef useGDAL
				std::vector<double> fvals = readValuesGDAL(src, 0, nrow(), 0, ncol());
				out.insert(out.end(), fvals.begin(), fvals.end());
				#endif // useGDAL
			}
		}
	} else { // read one lyr
		std::vector<unsigned> sl = findLyr(lyr);
		unsigned src=sl[0];
		if (source[src].memory) {
			size_t start = sl[1] * ncell();
			out = std::vector<double>(source[src].values.begin()+start, source[src].values.begin()+start+ncell());
		} else {
			#ifdef useGDAL
			out = readValuesGDAL(src, 0, nrow(), 0, ncol(), sl[1]);
			#endif // useGDAL
		}
	}
	return out;
}


bool SpatRaster::getValuesSource(size_t src, std::vector<double> &out) {
	
	unsigned n = nsrc();
	if (src > n) {
		return false;
	}
	if (source[src].memory) {
		out = std::vector<double>(source[src].values.begin(), source[src].values.end());
	} else {
		#ifdef useGDAL
		out = readValuesGDAL(src, 0, nrow(), 0, ncol());
		#endif // useGDAL
	}	
	return true;
}

