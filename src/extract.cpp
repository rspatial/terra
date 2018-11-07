#include "spatraster.h"


std::vector<double> SpatRaster::extract(SpatLayer v, std::string fun) {

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
			return cell; // for now
		} else {
			std::vector<double> x = vd.getD(0);
			std::vector<double> y = vd.getD(1);
			std::vector<unsigned> cols = colFromX(x);
			std::vector<unsigned> rows = rowFromY(y);
			out = readRowColGDAL(rows, cols);
		}
//		vd.cbind(df);
//		return vd;
	} else if (gtype == "lines") {
		
		
	} else {
		
		
	}
	return out;
}


