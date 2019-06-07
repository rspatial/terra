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

SpatRaster SpatRaster::init(std::string value, bool plusone, SpatOptions &opt) {

	SpatRaster out = geometry();
	
	std::vector<std::string> f {"row", "col", "cell", "x", "y", "chess"};
	bool test = std::find(f.begin(), f.end(), value) == f.end();
	if (test) {
		out.setError("not a valid init option");
		return out;
	}
 	if (!out.writeStart(opt)) { return out; }
	
	if (value == "row") {
		std::vector<double> v(ncol());
		size_t nr = nrow();
		for (size_t i = 0; i < nr; i++) {
			std::fill(v.begin(), v.end(), i+plusone);				
			if (!out.writeValues(v, i)) return out;
		}
	} else if (value == "col") {
		std::vector<double> col(ncol());
		double start = plusone ? 1 : 0;
		std::iota(col.begin(), col.end(), start);
		unsigned nr = nrow();
		for (unsigned i = 0; i < nr; i++) {
			if (!out.writeValues(col, i)) return out;
		}
	} else if (value == "cell") {
		std::vector<unsigned> col(ncol());
		std::iota(col.begin(), col.end(), 0);
		std::vector<unsigned> row(1);
		unsigned nr = nrow();
		for (unsigned i = 0; i < nr; i++) {
			row[0] = i;
			std::vector<double> v = cellFromRowCol(row, col);
			if (plusone) for(double& d : v) d=d+1;
			if (!out.writeValues(v, i)) return out;
		}
	} else if (value == "x") {
		std::vector<unsigned> col(ncol());
		std::iota(col.begin(), col.end(), 0);
		std::vector<double> x = xFromCol(col);
		unsigned nr = nrow();
		for (unsigned i = 0; i < nr; i++) {
			if (!out.writeValues(x, i)) return out;
		}
	} else if (value == "y") {
		std::vector<double> v(ncol());
		unsigned nr = nrow();
		for (unsigned i = 0; i < nr; i++) {
			double y = yFromRow(i);
			std::fill(v.begin(), v.end(), y);				
			if (!out.writeValues(v, i)) return out;
		}
	} else if (value == "chess") {
		std::vector<double> a(ncol());
		std::vector<double> b(ncol());
		size_t nr = nrow();
		size_t nc = ncol();
		a[0] = 1;
		b[0] = 0;
		for (size_t i=1; i<nc; i++) {
			bool test = i%2 == 0;
			a[i] = test;
			b[i] = !test;
		}		
		for (unsigned i=0; i<(nr-1); i=i+2) {
			if (!out.writeValues(a, i)) return out;
			if (!out.writeValues(b, i+1)) return out;
		}
		if (nr%2 == 0) {
			if (!out.writeValues(a, nr-2)) return out;			
			if (!out.writeValues(b, nr-1)) return out;			
		} else {
			if (!out.writeValues(a, nr-1)) return out;						
		}
	}

	out.writeStop();
	return(out);
}



SpatRaster SpatRaster::init(double value, SpatOptions &opt) {
	SpatRaster out = geometry();
 	if (!out.writeStart(opt)) { return out; }
	unsigned nc = ncol();
	std::vector<double> v(out.bs.nrows[0]*nc, value);	
	for (size_t i = 0; i < out.bs.n; i++) {
		if (i > 0 && i == (out.bs.n-1)) {
			v.resize(bs.nrows[i]*nc);
		}
		if (!out.writeValues(v, i)) return out;
	}
	out.writeStop();
	return(out);
}

