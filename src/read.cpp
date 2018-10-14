using namespace std;
#include "spat.h"
#include "rst_read.h"


bool SpatRaster::readStart() {
// for now assuming a single source
// will have to become a loop over sources
	if (!source[0].memory) {
		// open filestream
	}
	return true;
}

bool SpatRaster::readStop() {
	if (!source[0].memory) {
		// close filestream
	}
	return true;
}


std::vector<double> SpatRaster::readValues(unsigned row, unsigned nrows, unsigned col, unsigned ncols){

	row = std::min(std::max(unsigned(0), row), nrow-1);
	col = std::min(std::max(unsigned(0), col), ncol-1);
	nrows = std::max(unsigned(1), std::min(nrows, nrow-row));
	ncols = std::max(unsigned(1), std::min(ncols, ncol-col));
	std::vector<double> out;

	if (source[0].memory) {

		if (row==0 && nrows==nrow && col==0 && ncols==ncol) {
			out.insert(out.end(), source[0].values.begin(), source[0].values.end());
		} else {
			unsigned i, j;
			unsigned ncells = ncell();
			if (col==0 && ncols==ncol) {
				for (size_t lyr=0; lyr < nlyr(); lyr++) {
					unsigned add = ncells * lyr;
					i = add + row * ncol;
					j = i + nrows * ncol;
					out.insert(out.end(), source[0].values.begin()+i, source[0].values.begin()+j);
				}
			} else {
				unsigned endrow = row + nrows;
				unsigned endcol = col + ncols;
				for (size_t lyr=0; lyr < nlyr(); lyr++) {
					unsigned add = ncells * lyr;
					for (size_t r = row; r < endrow; r++) {
						i = add + r * ncol;
						out.insert(out.end(), source[0].values.begin()+i+col, source[0].values.begin()+i+endcol);
					}
				}
			}
		}
	} else {
		// read from file
		if (source[0].driver == "raster") {
			string file = source[0].filename;
			std::vector<unsigned> lyrs = {1};
			if (source[0].datatype == "FLT8S") {
				return readFLT8(file, 0, ncell(), nlyr(), order="BIL");
			} else {
				return readFLT4(file, "BIL", 0, ncell(), lyrs);
			}

		} else {
			return readValuesGDAL(row, nrows, col, ncols);
		}
	}
	return out;
}




std::vector<double>  SpatRaster::getValues() {
	if (source[0].memory) {
		return source[0].values;
	} else if (source[0].driver == "raster") {
		return readValues(0, nrow, 0, ncol);
	} else {
//		return readValuesGDAL(0, nrow, 0, ncol);
	}
}


