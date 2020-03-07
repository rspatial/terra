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
#include "recycle.h"
#include "NA.h"

std::vector<double> SpatRaster::cellFromXY (std::vector<double> x, std::vector<double> y) {
// size of x and y should be the same

	size_t size = x.size();
	std::vector<double> cells(size);

	double yr_inv = nrow() / (extent.ymax - extent.ymin);
	double xr_inv = ncol() / (extent.xmax - extent.xmin);

	for (size_t i = 0; i < size; i++) {
		// cannot use trunc here because trunc(-0.1) == 0
		long row = std::floor((extent.ymax - y[i]) * yr_inv);
		// points in between rows go to the row below
		// except for the last row, when they must go up
		if (y[i] == extent.ymin) {
			row = nrow()-1 ;
		}

		long col = floor((x[i] - extent.xmin) * xr_inv);
		// as for rows above. Go right, except for last column
		if (x[i] == extent.xmax) {
			col = ncol() - 1 ;
		}

		if (row < 0 || (unsigned)row >= nrow() || col < 0 || (unsigned)col >= ncol()) {
			cells[i] = NAN;
		} else {
			// result[i] = static_cast<int>(row) * ncols + static_cast<int>(col) + 1;
			cells[i] = row * ncol() + col;
		}
	}

	return cells;
}


double SpatRaster::cellFromXY (double x, double y) {
	std::vector<double> X = {x};
	std::vector<double> Y = {y};
	std::vector<double> cell = cellFromXY(X, Y);
	return  cell[0];
}


std::vector<double> SpatRaster::cellFromRowCol(std::vector<unsigned> row, std::vector<unsigned> col) {
	recycle(row, col);
	size_t n = row.size();
	std::vector<double> result(n);
	for (size_t i=0; i<n; i++) {
		result[i] = (row[i]<0 || row[i] >= nrow() || col[i]<0 || col[i] >= ncol()) ? NAN : row[i] * ncol() + col[i];
	}
	return result;
}


double SpatRaster::cellFromRowCol (unsigned row, unsigned col) {
	std::vector<unsigned> rows = {row};
	std::vector<unsigned> cols = {col};
	std::vector<double> cell = cellFromRowCol(rows, cols);
	return  cell[0];
}

std::vector<double> SpatRaster::cellFromRowColCombine(std::vector<unsigned> row, std::vector<unsigned> col) {
	recycle(row, col);
	size_t n = row.size();
	size_t nc = ncol();
	size_t nr = nrow();

	std::vector<double> result(n * n);
	for (size_t i=0; i<n; i++) {
		for (size_t j=0; j<n; j++) {
			result[i*n+j] = (row[i]<0 || row[i] >= nr || col[j]<0 || col[j] >= nc) ? NAN : row[i] * nc + col[j];
		}
	}
	return result;
}


double SpatRaster::cellFromRowColCombine(unsigned row, unsigned col) {
	return cellFromRowCol(row, col);
}


std::vector<double> SpatRaster::yFromRow(std::vector<unsigned> &row) {
	size_t size = row.size();
	std::vector<double> result( size );
	double ymax = extent.ymax;
	double yr = yres();
	for (size_t i = 0; i < size; i++) {
		result[i] = (row[i] < 0 || row[i] >= nrow() ) ? NAN : ymax - ((row[i]+0.5) * yr);
	}
	return result;
}

double SpatRaster::yFromRow (unsigned row) {
	std::vector<unsigned> rows = {row};
	std::vector<double> y = yFromRow(rows);
	return y[0];
}



std::vector<double> SpatRaster::xFromCol(std::vector<unsigned> &col) {
	size_t size = col.size();
	std::vector<double> result( size );
	double xmin = extent.xmin;
	double xr = xres();
	for (size_t i = 0; i < size; i++) {
		result[i] = (col[i] < 0 || col[i] >= ncol() ) ? NAN : xmin + ((col[i]+0.5) * xr);
	}
	return result;
}

double SpatRaster::xFromCol(unsigned col) {
	std::vector<unsigned> cols = {col};
	std::vector<double> x = xFromCol(cols);
	return x[0];
}

std::vector<unsigned> SpatRaster::colFromX(std::vector<double> &x) {
	size_t size = x.size();
	unsigned navalue = NA<unsigned>::value;

	std::vector<unsigned> result(size, navalue);
	double xmin = extent.xmin;
	double xmax = extent.xmax;
	double xr = xres();

	for (size_t i = 0; i < size; i++) {
		if (x[i] == xmax) {
			result[i] = ncol()-1;
		} else if (x[i] >= xmin || x[i] < xmax ) {
			result[i] =  trunc((x[i] - xmin) / xr);
		}
	}
	return result;
}


unsigned SpatRaster::colFromX(double x) {
	std::vector<double> X = {x};
	return colFromX(X)[0];
}


std::vector<unsigned> SpatRaster::rowFromY(std::vector<double> &y) {
	size_t ysize = y.size();
	unsigned navalue = NA<unsigned>::value;

	std::vector<unsigned> result(ysize, navalue);
	double ymin = extent.ymin;
	double ymax = extent.ymax;
	double yr = yres();

	for (size_t i = 0; i < ysize; i++) {
		if (y[i] == ymin) {
			result[i] = nrow()-1;
		} else if (y[i] > ymin || y[i] <= ymax ) {
			result[i] = trunc((ymax - y[i]) / yr);
		}
	}
	return result;
}


unsigned SpatRaster::rowFromY(double y) {
	std::vector<double> Y = {y};
	return rowFromY(Y)[0];
}


std::vector< std::vector<double> > SpatRaster::xyFromCell( std::vector<double> &cell) {
	size_t n = cell.size();
	double xmin = extent.xmin;
	double ymax = extent.ymax;
	double yr = yres();
	double xr = xres();
    double ncells = ncell();
    unsigned nc = ncol();
	std::vector< std::vector<double> > out(2, std::vector<double> (n, NAN) );
	for (size_t i = 0; i<n; i++) {
		if ((cell[i] < 0) || (cell[i] >= ncells)) continue;
        unsigned row = cell[i] / nc;
        unsigned col = cell[i] - (row * nc);
        out[0][i] = xmin + (col + 0.5) * xr;
        out[1][i] = ymax - (row + 0.5) * yr;
	}
	return out;
}


std::vector< std::vector<double> > SpatRaster::xyFromCell( double cell) {
	std::vector<double> Cell = {cell};
	return xyFromCell(Cell);
}


std::vector< std::vector<unsigned> > SpatRaster::rowColFromCell(std::vector<double> &cell) {
	size_t size = cell.size();
	unsigned navalue = NA<unsigned>::value;
	std::vector< std::vector<unsigned> > result(2, std::vector<unsigned> (size, navalue) );
	double nc = ncell();
	for (size_t i = 0; i < size; i++) {
		if ((cell[i] >= 0) || (cell[i] < nc )) {
			result[0][i] = trunc(cell[i]/ ncol());
			result[1][i] = (cell[i] - ((result[0][i]) * ncol()));
		}
	}
	return result;
}




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

