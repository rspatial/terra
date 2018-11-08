#include "spatraster.h"
#include "NA.h"

std::vector<double> SpatRaster::extractCell(std::vector<double> &cell) {

	unsigned n = cell.size();
	unsigned nc = ncell();
	std::vector<double> out(n * nlyr(), NAN);
	unsigned ns = nsrc();
	unsigned offset = 0;

	std::vector<std::vector<unsigned> > rc = rowColFromCell(cell);
	std::vector<unsigned> rows = rc[0];
	std::vector<unsigned> cols = rc[1];

	for (size_t src=0; src<ns; src++) {

		unsigned slyrs = source[src].nlyr;
		if (source[src].driver == "memory") {
			for (size_t i=0; i<slyrs; i++) {
				size_t j = i * nc;
				for (size_t k=0; k<n; k++) {
					if (!is_NA(cell[k]) && cell[k] >= 0 && cell[k] < nc) {
						out[offset + k] = source[src].values[j + cell[k]];
					}
				}
				offset += n;
			}
		} else {
		#ifdef useGDAL
			std::vector<double> srcout = readRowColGDAL(src, rows, cols);
			std::copy(srcout.begin(), srcout.end(), out.begin()+offset);
			offset += n;
        #endif
		}
	}
	return out;
}



std::vector<double> SpatRaster::extractLayer(SpatLayer v, std::string fun) {

	std::vector<double> out;
	std::vector<double> srcout;

	std::string gtype = v.type();
	if (gtype == "points") {
		SpatDataFrame vd = v.getGeometryDF();
		std::vector<double> x = vd.getD(0);
		std::vector<double> y = vd.getD(1);
		unsigned n = x.size();
		for (size_t src=0; src<nsrc(); src++) {
			if (source[src].driver == "memory") {
				std::vector<double> cell = cellFromXY(x, y);
				srcout = extractCell(cell);
			} else {
				std::vector<double> x = vd.getD(2);
				std::vector<double> y = vd.getD(3);
				std::vector<unsigned> rows = rowFromY(y);
				std::vector<unsigned> cols = colFromX(x);
				#ifdef useGDAL
				srcout = readRowColGDAL(src, rows, cols);
				#endif
			}
		}
		out.insert(out.end(), srcout.begin(), srcout.end());
		
	} else if (gtype == "lines") {


	} else {


	}
	return out;
}


