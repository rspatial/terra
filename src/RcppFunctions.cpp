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

#include <Rcpp.h>
#include "spatRaster.h"


/*
NumericMatrix getValuesM(SpatRaster* r) {
	NumericMatrix x(r->ncell(), r->nlyr() );
	std::vector<double> v;
	v = r->getValues();
	std::copy(v.begin(), v.end(), x.begin());
	return(x);
}
*/

Rcpp::List getBlockSizeR(SpatRaster* r, unsigned n) { 
    BlockSize bs = r->getBlockSize(n);
	Rcpp::List L = Rcpp::List::create(Rcpp::Named("row") = bs.row, Rcpp::Named("nrows") = bs.nrows, Rcpp::Named("n") = bs.n);
	return(L);
}



Rcpp::List getAttributes(SpatVector* v) {
	unsigned n = v->ncol();
	Rcpp::List out(n);	
	std::vector<unsigned> itype = v->getItype();
	for (size_t i=0; i < n; i++) {
		if (itype[i] == 0) {
			out[i] = v->getDv(i);
		} else if (itype[i] == 1) {
			out[i] = v->getIv(i);
		} else {
			out[i] = v->getSv(i);
		}
	}	
	// todo: deal with NAs in int and str
	return out;
// df is nice, but no no of variables is <= 20, and no "stringsAsFactors"=false
//	Rcpp::DataFrame result(out);
//	result.attr("names") = v->names();
//	return result;
}


Rcpp::DataFrame getGeometry(SpatVector* v) {
	SpatDataFrame df = v->getGeometryDF();

	Rcpp::DataFrame out = Rcpp::DataFrame::create(
			Rcpp::Named("id") = df.iv[0], 
			Rcpp::Named("part") = df.iv[1], 
			Rcpp::Named("x") = df.dv[0],
			Rcpp::Named("y") = df.dv[1],
			Rcpp::Named("hole") = df.iv[2]
	);
	return out;
}


SpatRaster rcppReclassify(SpatRaster* x, Rcpp::NumericMatrix rcl, unsigned right, bool lowest, SpatOptions opt) {
	
	unsigned nc = rcl.ncol();
	unsigned nr = rcl.nrow();
	printf("nc %u nr %u \n", nr, nr);
	std::vector< std::vector<double>> rc(nc);
	for (size_t c=0; c<nc; c++) {
		for (size_t r=0; r<nr; r++) {
			rc[c].push_back(rcl(r,c));
		}
	}

	SpatRaster out = x->reclassify(rc, right, lowest, opt);
	return out;
}


