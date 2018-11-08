#include "spatraster.h"

std::vector<double> SpatRaster::extractCell(std::vector<double> &cell) {
	std::vector<double> out;
	//for (size_t i=0; i<nsrc(); i++) {
		size_t i = 0;
		if (source[i].driver == "memory") {
			size_t n = cell.size();
			size_t nc = ncell();
			unsigned nlyrs = source[i].nlyr;
			out.resize(n * nlyrs);
			for (size_t j=0; j<nlyrs; j++) {
				size_t jj = j * nc;
				size_t nn = j * n;
				for (size_t k=0; k<n; k++) {
					size_t kk = nn + k;
					size_t c = jj + cell[k];
					out[kk] = source[i].values[c];
				}
			}
		} else {
			std::vector<std::vector<unsigned> > rc = rowColFromCell(cell);
			std::vector<unsigned> rows = rc[0];
			std::vector<unsigned> cols = rc[1];
			out = readRowColGDAL(rows, cols);
		}
	//}
	return out;
}



std::vector<double> SpatRaster::extractLayer(SpatLayer v, std::string fun) {

	std::vector<double> out;
	std::string gtype = v.type();
	if (gtype == "points") {
		//for (size_t i=0; i<nsrc(); i++) {
		size_t i = 0;
		SpatDataFrame vd = v.getGeometryDF();
		std::vector<double> x = vd.getD(0);
		std::vector<double> y = vd.getD(1);
		if (source[i].driver == "memory") {
			std::vector<double> cell = cellFromXY(x, y);
			out = extractCell(cell);
		} else {
			std::vector<double> x = vd.getD(2);
			std::vector<double> y = vd.getD(3);
			std::vector<unsigned> rows = rowFromY(y);
			std::vector<unsigned> cols = colFromX(x);
			out = readRowColGDAL(rows, cols);
		}
	} else if (gtype == "lines") {
		
		
	} else {
		
		
	}
	return out;
}


