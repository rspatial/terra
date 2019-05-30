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

std::vector<std::vector<double>> SpatRaster::adjacent(std::vector<double> cells, std::string directions, bool include) {

	unsigned n = cells.size();
	std::vector<std::vector<double>> out(n);

	std::vector<std::string> f {"rook", "queen", "bishop", "16"};
	if (std::find(f.begin(), f.end(), directions) == f.end()) {
        setError("argument directions is not valid");
        return(out);
	}

	std::vector<std::vector<unsigned>> rc = rowColFromCell(cells);
	std::vector<unsigned> r = rc[0];
	std::vector<unsigned> c = rc[1];
	bool globlatlon = is_global_lonlat();
    unsigned nc = ncol();
    unsigned lc = nc-1;
    std::vector<unsigned> cols, rows;
	if (directions == "rook") {
		for (size_t i=0; i<n; i++) {
			rows = {r[i]-1, r[i]   , r[i]  , r[i]+1};
            cols = {c[i]  , c[i]-1 , c[i]+1, c[i]};
            if (globlatlon) {
                if (c[i]==0) {
                    cols[1] = lc;
                } else if (c[i]==lc) {
                    cols[2] = 0;
                }
            }
            if (include) {
                rows.push_back(r[i]);
                cols.push_back(c[i]);
            }
			out[i] = cellFromRowCol(rows, cols);
			//std::sort(out[i].begin(), out[i].end());
		}
	} else if (directions == "queen") {
		for (size_t i=0; i<n; i++) {
            rows = {r[i]-1, r[i]-1, r[i]-1, r[i], r[i], r[i]+1, r[i]+1, r[i]+1};
            cols = {c[i]-1, c[i], c[i]+1, c[i]-1, c[i]+1, c[i]-1, c[i], c[i]+1};
            if (globlatlon) {
                if (c[i]==0) {
                    cols = {lc, c[i], c[i]+1, lc, c[i]+1, lc, c[i], c[i]+1};
                } else if (c[i]==lc) {
                    cols = {c[i]-1, c[i], 0, c[i]-1, 0, c[i]-1, c[i], 0};
                }
            }
            if (include) {
                rows.push_back(r[i]);
                cols.push_back(c[i]);
            }
			out[i] = cellFromRowCol(rows, cols);
			//std::sort(out[i].begin(), out[i].end());
		}
	} else if (directions == "bishop") {
		for (size_t i=0; i<n; i++) {
            rows = {r[i]-1, r[i]-1, r[i]+1, r[i]+1};
            cols = {c[i]-1, c[i]+1, c[i]-1, c[i]+1};
            if (globlatlon) {
                if (c[i]==0) {
                    cols = {lc, c[i]+1, lc, c[i]+1};
                } else if (c[i]==lc) {
                    cols = {c[i]-1, 0, c[i]-1, 0};
                }
            }
            if (include) {
                rows.push_back(r[i]);
                cols.push_back(c[i]);
            }
			out[i] = cellFromRowCol(rows, cols);
			//std::sort(out[i].begin(), out[i].end());
		}
	} else if (directions == "16") {
		for (size_t i=0; i<n; i++) {
            rows = {r[i]-2, r[i]-2, r[i]-1, r[i]-1, r[i]-1, r[i]-1, r[i]-1, r[i]  , r[i]  , r[i]+1, r[i]+1, r[i]+1, r[i]+1, r[i]+1, r[i]+2, r[i]+2};
            cols = {c[i]-1, c[i]+1, c[i]-2, c[i]-1, c[i],   c[i]+1, c[i]+2, c[i]-1, c[i]+1, c[i]-2, c[i]-1, c[i]  , c[i]+1, c[i]+2, c[i]-1, c[i]+1};
            if (globlatlon) {
                if ((c[i]==0) | (c[i]==1)) {
                    for (size_t j=0; j<16; j++) {
                        cols[j] = (cols[j] < 0) ? nc-cols[j] : cols[j];
                    }
                } else if (c[i]==nc) {
                    for (size_t j=0; j<16; j++) {
                        cols[j] = (cols[j] > lc) ? cols[j]-nc : cols[j];
                    }
                }
            }
            if (include) {
                rows.push_back(r[i]);
                cols.push_back(c[i]);
            }
			out[i] = cellFromRowCol(rows, cols);
			//std::sort(out[i].begin(), out[i].end());
		}
	}
	return(out);
}

/*
std::vector<std::vector<double>> SpatRaster::adjacent(std::vector<double> cells, std::string directions, bool include) {

	double xr = xres();
	double yr = yres();
	unsigned n = cells.size();
	std::vector<std::vector<double>> xy = xyFromCell(cells);
	std::vector<double> x = xy[0];
	std::vector<double> y = xy[1];
	std::vector<std::vector<double>> out(n);

	if (directions == "4") {
		if (include) {
			for (size_t i=0; i<n; i++) {
				std::vector<double> adjx = {x[i], x[i]-xr, x[i]+xr, x[i]   , x[i]    };
				std::vector<double> adjy = {y[i], y[i]   , y[i],    y[i]-yr, y[i]+yr };
				out[i] = cellFromXY(adjx, adjy);
				std::sort(out[i].begin(), out[i].end());
			}
		} else {
			for (size_t i=0; i<n; i++) {
				std::vector<double> adjx = {x[i]-xr, x[i]+xr, x[i]   , x[i]    };
				std::vector<double> adjy = {y[i]   , y[i],    y[i]-yr, y[i]+yr };
				out[i] = cellFromXY(adjx, adjy);
				std::sort(out[i].begin(), out[i].end());
			}
		}
	}
	return(out);
}
*/

