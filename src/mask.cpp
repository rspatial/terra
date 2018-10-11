using namespace std;
#include <vector>
#include "spat.h"






SpatRaster SpatRaster::mask(SpatRaster x, string filename, bool overwrite) {

// check for size; need for recycling
	//SpatRaster out = *this;
	SpatRaster out = geometry();

//	out.source.resize(1);
 //   our.source[0].nlyr = 1;
//	out.values.resize(0);
  	out.writeStart(filename, overwrite);
	readStart();
	x.readStart();
	std::vector<double> v, m;
	for (size_t i = 0; i < out.bs.n; i++) {
		v = readValues(out.bs.row[i], out.bs.nrows[i], 0, ncol);
		m = x.readValues(out.bs.row[i], out.bs.nrows[i], 0, ncol);
		for (size_t i=0; i < v.size(); i++) {
			if (std::isnan(m[i])) {
				v[i] = NAN;
			}
		}
		out.writeValues(v, out.bs.row[i]);
	}
	out.writeStop();
	readStop();
	x.readStop();
	return(out);
}


/*

std::vector<std::vector<double> > matrix(int nrow, int ncol) {
	std::vector<std::vector<double> > m (nrow, std::vector<double>(ncol));
	return(m);
}


int main() {
	std::vector<vector<double> > d = matrix(10, 2);
	std::vector<double> m (10);
	m[1] = 1;
	m[5] = 1;
	d = mask(d, m, 1, 9, false);
	for (int i=0; i < d.size(); i++) {
		for (int j=0; j < d[0].size(); j++) {
			std::cout << ' ' << d[i][j];
		}
		std::cout << '\n';
	}
	std::cout << '\n';
}
*/
