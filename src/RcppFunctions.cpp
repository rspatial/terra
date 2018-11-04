#include <Rcpp.h>
#include "spatraster.h"


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



Rcpp::List getAttributes(SpatLayer* v) {
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


Rcpp::DataFrame getGeometry(SpatLayer* v) {
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


