
#include "spatRaster.h"

/* Compute surface area, using method from:
Jeff S. Jenness, 2004. Calculating Landscape Surface Area from Digital Elevation Models. Wildlife Society Bulletin 32(3):829-839. http://www.jstor.org/stable/3784807

With edge adjustments. Adapted from code in R package "sp" by Barry Rowlingson 2010 <b.rowlingson@lancaster.ac.uk>
*/

inline double height(const std::vector<double> &heights, const size_t &ncols, const size_t &row, const size_t &col) {
	return heights[col + ncols*row];
}

inline double triarea(const double &a, const double &b, const double &c) {
	// triangle area given side lengths
	double s = (a+b+c) / 2.0;
	return sqrt(s * (s-a) * (s-b) * (s-c));
}


void sarea(std::vector<double> &heights, const size_t &nrow, const size_t &ncol, const double &w, const double &h, 
	std::vector<double> &sa, const bool &startadd, const bool &endadd) {

// given an nx by ny matrix of heights with single-cell edge border, compute the surface area.

// point values
	double z1, z2, z3;
// side lengths
	double l1, l2, l3; 
// diagonal length
	double s2 = sqrt((w*w)+(h*h));

// offsets to neighbours
	int dxv[] = {-1, 0, 1, 1, 1, 0, -1, -1, -1};
	int dyv[] = {-1, -1, -1, 0, 1, 1, 1, 0, -1};

// triangle side lengths
// first the radial sides
	double side[] = {s2, h, s2, w, s2, h, s2, w, s2};
// outer edges
	double l3v[] = {w, w, h, h, w, w, h, h};

	size_t outsize = heights.size() - ncol * (startadd + endadd);
	sa = std::vector<double>(outsize, NAN);
	size_t cellI = startadd ? 0 : ncol;

	for (size_t i=1; i<(nrow-1); i++){
		cellI++;
		for (size_t j=1; j<(ncol-1); j++){
			z1 = height(heights, ncol, i, j);
			if (!std::isnan(z1)) {
				double cellArea = 0;
				for (size_t tri=0; tri<8; tri++){
					z2 = height(heights, ncol, i+dxv[tri], j+dyv[tri]);
					// replace missing adjacent values with the current cell value
					if (std::isnan(z2)) z2=z1;
					z3 = height(heights, ncol, i+dxv[tri+1], j+dyv[tri+1]);
					if (std::isnan(z3)) z3=z1;
					l1 = 0.5 * sqrt(side[tri] * side[tri] + (z1-z2)*(z1-z2));
					l2 = 0.5 * sqrt(side[tri+1] * side[tri+1] + (z1-z3)*(z1-z3));
					l3 = 0.5 * sqrt(l3v[tri] * l3v[tri] + (z2-z3)*(z2-z3));
					cellArea += triarea(l1, l2, l3);
				}
//				sa[cellI] = cellArea;
				sa[i*ncol+j] = cellArea;
			}
			cellI++;
		}	
		cellI++;
	}
}



SpatRaster SpatRaster::surfaceArea(SpatOptions &opt) {

	SpatRaster out = geometry(1, false);
	if (is_lonlat()) {
		out.setError("not yet implemented for lonlat data");
		return out;
	}
	
	if (!hasValues()) {
		out.setError("cannot compute surfaceArea for a raster with no values");
		return out;
	}

	size_t nl = nlyr();
	if (nl != 1) {
		out.setError("can only compute surfaceArea for a single raster layer");
		return out;		
	}

	if (!readStart()) {
		out.setError(getError());
		return(out);
	}

  	if (!out.writeStart(opt, filenames())) {
		readStop();
		return out;
	}

	BlockSize cbs = out.bs;
	for (size_t i = 0; i < (cbs.n-1); i++) {
		cbs.nrows[i] += 1;
	}	
	for (size_t i = 1; i < cbs.n; i++) {
		cbs.row[i] -= 1;
		cbs.nrows[i] += 1;
	}
	
	std::vector<double> wh = resolution();
	for (size_t i = 0; i < cbs.n; i++) {
		std::vector<double> v;
		readBlock(v, cbs, i);
		std::vector<double> sa;
		sarea(v, cbs.nrows[i], ncol(), wh[0], wh[1], sa, i>0, i<(cbs.n-1));
		if (!out.writeBlock(sa, i)) return out;
	}
	readStop();
	out.writeStop();
	return(out);
}

