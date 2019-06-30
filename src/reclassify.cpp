// Copyright (c) 2018-2019  Robert J. Hijmans
//
// This file is part of the "spat" library.
//
// spat is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// spat is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with spat. If not, see <http://www.gnu.org/licenses/>.

#include "spatRaster.h"
#include <cmath>

void reclass_vector(std::vector<double> &v, std::vector<std::vector<double>> rcl, unsigned doright, bool lowest, bool othNA) {

	size_t nc = rcl.size(); // should be 2 or 3
	if (nc == 2) {
		doright = 3;
	}

	bool right;
	bool leftright = false;
	if ((doright != 0) & (doright != 1)) {
		leftright = true;
	} else if (doright == 0) {
		right = true;
	} else {
		right = false;
	}

//	bool hasNA = false;
	double NAval = NAN;

	size_t n = v.size();
	unsigned nr = rcl[0].size();

	// "is - becomes"
	if (nc == 2) {
		for (size_t i=0; i<n; i++) {
			if (std::isnan(v[i])) {
				v[i] = NAval;
			} else {
				bool found = false;
				for (size_t j=0; j<nr; j++) {
					if (v[i] == rcl[0][j]) {
						v[i] = rcl[1][j];
						found = true;
						break;
					}
				}
				if ((othNA) & (!found)) {
					v[i] = NAval;
				}
			}
		}

	// "from - to - becomes"
	} else if (leftright) {   // interval closed at left and right

		for (size_t i=0; i<n; i++) {
			if (std::isnan(v[i])) {
				v[i] = NAval;
			} else {
				bool found = false;
				for (size_t j=0; j<nr; j++) {
					if ((v[i] >= rcl[0][j]) & (v[i] <= rcl[1][j])) {
						v[i] = rcl[1][j];
						found = true;
						break;
					}
				}
				if ((othNA) & (!found))  {
					v[i] = NAval;
				}				
			}
		}
		} else if (right) {  // interval closed at right
			if (lowest) {  // include lowest value (left) of interval

			double lowval = rcl[0][0];
			double lowres = rcl[2][0];
			for (size_t i=1; i<nr; i++) {
				if (rcl[0][i] < lowval) {
					lowval = rcl[0][i];
					lowres = rcl[2][i];
				}
			}

			for (size_t i=0; i<n; i++) {
				if (std::isnan(v[i])) {
					v[i] = NAval;
				} else if (v[i] == lowval) {
					v[i] = lowres;
				} else {
					bool found = false;
					for (size_t j=0; j<nr; j++) {
						if ((v[i] > rcl[0][j]) & (v[i] <= rcl[1][j])) {
							v[i] = rcl[2][j];
							found = true;
							break;
						}
					}
					if  ((othNA) & (!found))  {
						v[i] = NAval;
					}
				}
			}

		} else { // !dolowest
				for (size_t i=0; i<n; i++) {
				if (std::isnan(v[i])) {
					v[i] = NAval;
				} else {
					bool found = false;
					for (size_t j=0; j<nr; j++) {
						if ((v[i] > rcl[0][j]) & (v[i] <= rcl[1][j])) {
							v[i] = rcl[2][j];
							found = true;
							break;
						}
					}
					if  ((othNA) & (!found))  {
						v[i] = NAval;
					}
				}
			}
		}

	} else { // !doright

		if (lowest) { // which here means highest because right=FALSE

			double lowval = rcl[1][0];
			double lowres = rcl[2][0];
			for (size_t i=0; i<nr; i++) {
				if (rcl[1][i] > lowval) {
					lowval = rcl[1][i];
					lowres = rcl[2][i];
				}
			}

			for (size_t i=0; i<n; i++) {
				if (std::isnan(v[i])) {
					v[i] = NAval;
				} else if (v[i] == lowval) {
					v[i] = lowres;
				} else {
					bool found = false;
					for (size_t j=0; j<nr; j++) {
						if ((v[i] >= rcl[0][j]) & (v[i] < rcl[1][j])) {
							v[i] = rcl[2][j];
							found = true;
							break;
						}
					}
					if  ((othNA) & (!found))  {
						v[i] = NAval;
					}
				}
			}

		} else { //!dolowest

			for (size_t i=0; i<n; i++) {
				if (std::isnan(v[i])) {
					v[i] = NAval;
				} else {
					bool found = false;
					for (size_t j=0; j<nr; j++) {
						if ((v[i] >= rcl[0][j]) & (v[i] < rcl[1][j])) {
							v[i] = rcl[2][j];
							found = true;
							break;
						}
					}
					if  ((othNA) & (!found))  {
						v[i] = NAval;
					}
				}
			}
		}
	}
}


SpatRaster SpatRaster::reclassify(std::vector<std::vector<double>> rcl, unsigned right, bool lowest, bool othersNA, SpatOptions &opt) {

	SpatRaster out = geometry();
	size_t nc = rcl.size();
	size_t nr = rcl[0].size();
	if (nc < 2 || nc > 3 || nr < 1) {
		out.setError("reclassification matrix must have 2 or 3 columns, and at least one row");
		return out;
	}
	for (size_t i=0; i<nc; i++) {
		if (rcl[i].size() != nr) {
			out.setError("reclassification matrix is not rectangular");
			return out;
		}
	}
	if (nc == 3) {
		for (size_t i=0; i<nr; i++) {
			if (rcl[0][i] > rcl[1][i]) {
				out.setError("'from' smaller than 'to' in reclassification matrix");
				return out;
			}
		}
	}

  	if (!out.writeStart(opt)) { return out; }
	readStart();
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> v = readBlock(out.bs, i);
		reclass_vector(v, rcl, right, lowest, othersNA);
		if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i], 0, ncol())) return out;		
	}
	readStop();
	out.writeStop();
	return(out);

}


SpatRaster SpatRaster::reclassify(std::vector<double> rcl, unsigned nc, unsigned right, bool lowest, bool othersNA, SpatOptions &opt) {
	
	SpatRaster out;
	std::vector< std::vector<double>> rc(nc);
	if ((rcl.size() % nc) != 0) {
		out.setError("incorrect length of reclassify matrix");
		return(out);
	}
	unsigned i  = rcl.size() / nc;
	
	if (nc == 1) {
		rc[0] = rcl;
	} else if (nc==2) {
		rc[0] = std::vector<double>(rcl.begin(), rcl.begin()+i);
		rc[1] = std::vector<double>(rcl.begin()+i, rcl.end());
	} else if (nc==3) {
		rc[0] = std::vector<double>(rcl.begin(), rcl.begin()+i);
		rc[1] = std::vector<double>(rcl.begin()+i, rcl.begin()+2*i);
		rc[2] = std::vector<double>(rcl.begin()+2*1, rcl.end());
	} else {
		out.setError("incorrect number of columns in reclassify matrix");
		return(out);
	}
	out = reclassify(rc, right, lowest, othersNA, opt);
	return out;
}