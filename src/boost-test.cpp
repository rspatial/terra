/*
#include "boost/multi_array.hpp"
#include <cassert>


class BoostAr {

typedef boost::multi_array<double, 3> array3D;
typedef array3D::index index;

	public:
	  // Create a 3D array that is 3 x 4 x 2
	  
		BoostAr() {
			array3D A(boost::extents[3][4][2]);

	  // Assign values to the elements
		  int values = 0;
		  for(index i = 0; i != 3; ++i) 
			for(index j = 0; j != 4; ++j)
			  for(index k = 0; k != 2; ++k)
				A[i][j][k] = values++;
	}
	
	int verifyValues() {
	  // Verify values
	  int verify = 0;
	  for(index i = 0; i != 3; ++i) 
		for(index j = 0; j != 4; ++j)
		  for(index k = 0; k != 2; ++k)
			assert(A[i][j][k] == verify++);

		return verify;
	}
};



typedef boost::multi_array<double, 3> array3D;
typedef array3D::index index3D;

array3D SpatRaster::readValues3D(unsigned row, unsigned nrows, unsigned col, unsigned ncols){
	// to use zero based indexing at the C level
	// have to make sure that this happens at the right places
	// and not more than once ..
	
	// for now
	unsigned nlyrs = nlyr();
	
	unsigned nr = std::min(nrows, nrow-row);
	unsigned nc = std::min(ncols, ncol-col);
	if ((nr != nrows) || (nc != ncols)) {
		// message
		nrows = nr;
		ncols = nc;
	}
	unsigned endrow = row+nrows-1;
	unsigned endcol = col+ncols-1;
	
	array3D out(boost::extents[nlyrs][nrows][ncols]);

	if (source.memory[0]) {
		size_t ii = 0;
		size_t jj = 0; 
		size_t kk = 0;
		for(size_t i = 0; i != nlyrs; ++i) {
			for(size_t j = row; j != endrow; ++j) {
				for(size_t k = col; k != endcol; ++k) {
					out[ii][jj][kk] = values3D[i][j][k];
					kk++;
				}
				jj++;
			}
			ii++;
		}
	} else {
	// read from file
		string file = source.filename[0];
		
		// for now, read all values 
		std::vector<double> v = readFLT4(file, 0, size());

//		std::copy(m.data(), m.data()+m.num_elements(), v.begin());
		
		boost::multi_array<double, 3> out(boost::extents[nlyrs][nrows][ncols]);

		for(size_t i = 0; i < nlyrs; ++i) {
			for(size_t j = row; j < endrow; ++j) {
				for(size_t k = col; k < endcol; ++k) {
					out[i][j][k] = v[i*nrow*ncol + j*ncol + k];
				}
			}
		}
	}
	
	return(out);	
}

std::vector<double> SpatRaster::readValues3D2R(unsigned row, unsigned nrows, unsigned col, unsigned ncols){
	
	boost::multi_array<double, 3> m = SpatRaster::readValues3D(row, nrows, col, ncols);
	std::vector<double> v(m.num_elements());
	std::copy(m.data(), m.data()+m.num_elements(), v.begin());
	return v;
}



*/