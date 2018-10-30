#include <Rcpp.h>
#include "spat.h"

using namespace Rcpp;

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
	List L = List::create(Named("row") = bs.row, Named("nrows") = bs.nrows, Named("n") = bs.n);
	return(L);
}


//Rcpp::DataFrame getDataFrame(SpatVector* v) {
Rcpp::List getAttributes(SpatVector* v) {
	unsigned n = v->ncol();
	List out(n);	
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

Rcpp::NumericVector getGeometry(SpatVector* v) {
	// for points only
	std::vector<double> xy = v->pts.x;
	xy.insert( xy.end(), v->pts.y.begin(), v->pts.y.end() );
	NumericVector out = wrap(xy);	
	out.attr("dim") = Dimension(v->pts.x.size(), 2);
	return out;
}

