using namespace std;
#include "spat.h"
#include "rst_read.h"


bool SpatRaster::readStart() {
// for now assuming a single source
// will have to become a loop over sources
	if (!source.memory[0]) {
		// open filestream
	}
	return true;
}

bool SpatRaster::readStop() {
	if (!source.memory[0]) {
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
	
	if (source.memory[0]) {

		if (row==0 && nrows==nrow && col==0 && ncols==ncol) {
			out.insert(out.end(), values.begin(), values.end());
		} else { 
			unsigned i, j;
			unsigned ncells = ncell();
			if (col==0 && ncols==ncol) {
				for (size_t lyr=0; lyr < nlyr; lyr++) {
					unsigned add = ncells * lyr;
					i = add + row * ncol;
					j = i + nrows * ncol;
					out.insert(out.end(), values.begin()+i, values.begin()+j);
				}				
			} else {
				unsigned endrow = row + nrows;
				unsigned endcol = col + ncols;
				for (size_t lyr=0; lyr < nlyr; lyr++) {
					unsigned add = ncells * lyr;
					for (size_t r = row; r < endrow; r++) {
						i = add + r * ncol;
						out.insert(out.end(), values.begin()+i+col, values.begin()+i+endcol);
					}
				}
			}
		} 
	} else {
		// read from file
		if (source.driver[0] == "raster") {
			string file = source.filename[0];
			if (source.datatype[0] == "FLT8S") {
				return readFLT8(file, 0, nlyr * ncell());
			} else {
				return readFLT4(file, 0, nlyr * ncell());
			}

		} else {
			return readValuesGDAL(row, nrows, col, ncols);			
		}
	}
	return out;
}




std::vector<double>  SpatRaster::getValues() { 
	if (source.memory[0]) {
		return values; 
	} else {
		if (source.driver[0] == "raster") {
			return readValues(0, nrow, 0, ncol);
		} else {
			return readValuesGDAL(0, nrow, 0, ncol);			
		}
	}
}

