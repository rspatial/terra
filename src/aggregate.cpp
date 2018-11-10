/* Robert Hijmans, October 2014 */

#include <vector>
#include <limits>
#include <cmath>
#include "spatraster.h"


template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T>>& v) {
    std::size_t total_size = 0;
    for (const auto& sub : v)
        total_size += sub.size(); // I wish there was a transform_accumulate
    std::vector<T> result;
    result.reserve(total_size);
    for (const auto& sub : v)
        result.insert(result.end(), sub.begin(), sub.end());
    return result;
}

std::vector<double> flat(std::vector<std::vector<double>> v) {
    unsigned s1 = v.size();
    unsigned s2 = v[0].size();

	std::size_t s = s1 * s2;
    std::vector<double> result(s);
    for (size_t i=0; i<s1; i++) {
		for (size_t j=0; j<s2; j++) {
			result[i*s2+j] = v[i][j];
		}
	}
	return result;
}


std::vector<unsigned> SpatRaster::get_aggregate_dims( std::vector<unsigned> fact ) {
	// fact has aggregation factors in the three dimensions
	unsigned fs = fact.size();
	fact.resize(6);
	if (fs == 1) {
		fact[1] = fact[0];
		fact[2] = 1;
	} else if (fs == 2) {
		fact[2] = 1;
	}
	// int dy = dim[0], dx = dim[1], dz = dim[2];
	fact[0] = std::max(unsigned(1), std::min(fact[0], nrow));
	fact[1] = std::max(unsigned(1), std::min(fact[1], ncol));
	fact[2] = std::max(unsigned(1), std::min(fact[2], nlyr()));
	// new dimensions: rows, cols, lays
	fact[3] = std::ceil(double(nrow) / fact[0]);
	fact[4] = std::ceil(double(ncol) / fact[1]);
	fact[5] = std::ceil(double(nlyr()) / fact[2]);
	return fact;
}


std::vector<std::vector<double> > SpatRaster::get_aggregates(std::vector<unsigned> dim) {
	// blocks per row (=ncol), col (=nrow)

//	dim = get_aggregate_dims(dim);

	unsigned dy = dim[0], dx = dim[1], dz = dim[2];
	unsigned bpC = dim[3];
	unsigned bpR = dim[4];
	unsigned bpL = bpR * bpC;

	// new number of layers
	unsigned newNL = dim[5];

	// new number of rows, adjusted for additional (expansion) rows
	unsigned adjnr = bpC * dy;

	// number of aggregates
	unsigned nblocks = (bpR * bpC * newNL);
	// cells per aggregate
	unsigned blockcells = dx * dy * dz;

	// output: each row is a block
	std::vector< std::vector<double> > a(nblocks, std::vector<double>(blockcells, std::numeric_limits<double>::quiet_NaN()));

	for (unsigned b = 0; b < nblocks; b++) {
		unsigned lstart = dz * (b / bpL);
		unsigned rstart = (dy * (b / bpR)) % adjnr;
		unsigned cstart = dx * (b % bpR);

		unsigned lmax   = std::min(nlyr(), (lstart + dz));
		unsigned rmax   = std::min(nrow, (rstart + dy));
		unsigned cmax   = std::min(ncol, (cstart + dx));

		unsigned f = 0;
		unsigned nc = ncell();
		for (unsigned j = lstart; j < lmax; j++) {
			unsigned lj = j * nc;
			for (unsigned r = rstart; r < rmax; r++) {
				unsigned cell = lj + r * ncol;
				for (unsigned c = cstart; c < cmax; c++) {
					a[b][f] = source[0].values[cell + c];
					f++;
				}
			}
		}
	}
	return(a);
}


SpatRaster SpatRaster::aggregate(std::vector<unsigned> fact, std::string fun, bool narm, std::string filename, std::string format, std::string datatype, bool overwrite) {

//std::vector<double> SpatRaster::aggregate(std::vector<unsigned> fact, bool narm, string fun, string filename) {

//	fact = get_aggregate_dims(fact);

	unsigned f = 1, mean = 0; // sum
	if (fun == "mean") {
		f = 1;
		mean = 1;
	} else if (fun == "min") {
		f = 2;
	} else if (fun == "max") {
		f = 3;
	}


	double xmax = extent.xmin + fact[4] * yres();
	double ymin = extent.ymax - fact[5] * xres();
	SpatExtent e = SpatExtent(extent.xmin, xmax, ymin, extent.ymax);
	SpatRaster r = SpatRaster(fact[3], fact[4], fact[5], e, crs);

	if (!source[0].hasValues) { return r; }

	// output: each row is a new cell
	std::vector< std::vector<double> > v(fact[5], std::vector<double>(fact[3]*fact[4]));

	// get the aggregates
	std::vector<std::vector< double > > a = get_aggregates(fact);

	int nblocks = a.size();
	int naggs = a[0].size();

	for (int i = 0; i < nblocks; i++) {
		unsigned row = (i / ncol) % nrow;
		unsigned col = i % ncol;
		unsigned cell = row * ncol + col;
		unsigned lyr = std::floor(i / (nrow * ncol));

		double x = 0;
		if (f==2) { // min
			x = std::numeric_limits<double>::infinity();
		} else if (f==3) { // max
			x = - std::numeric_limits<double>::infinity() ;
		}

		double cnt = 0;
		for (int j = 0; j < naggs; j++) {
			if (std::isnan(a[i][j])) {
				if (!narm) {
					x = NAN;
					goto breakout;
				}
			} else {
				if (f==2) { // min
					x = std::min(x, a[i][j]);
				} else if (f==3) { // max
					x = std::max(x, a[i][j]);
				} else { // sum or mean
					x += a[i][j];
				}
				cnt++;
			}
		}
		if (cnt > 0) {
			if (mean) {
				x = x / cnt;
			}
		} else {
			x = NAN;
		}
		breakout:
		v[lyr][cell] = x;
	}


	r.setValues( flat(v) );

	return(r);

}


