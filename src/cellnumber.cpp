// Copyright (c) 2018-2020  Robert J. Hijmans
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
//#include "NA.h"

std::vector<double> SpatRaster::cellFromXY (std::vector<double> x, std::vector<double> y) {
// size of x and y should be the same

	size_t size = x.size();
	std::vector<double> cells(size);

	SpatExtent extent = getExtent();
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

		long col = std::floor((x[i] - extent.xmin) * xr_inv);
		// as for rows above. Go right, except for last column
		if (x[i] == extent.xmax) {
			col = ncol() - 1 ;
		}
		long nr = nrow();
		long nc = ncol();
		if (row < 0 || row >= nr || col < 0 || col >= nc) {
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


std::vector<double> SpatRaster::cellFromRowCol(std::vector<int_64> row, std::vector<int_64> col) {
	recycle(row, col);
	size_t n = row.size();
	std::vector<double> result(n);
	int_64 nr = nrow();
	int_64 nc = ncol();
	for (size_t i=0; i<n; i++) {
		result[i] = (row[i]<0 || row[i] >= nr || col[i]<0 || col[i] >= nc) ? NAN : row[i] * nc + col[i];
	}
	return result;
}


double SpatRaster::cellFromRowCol (int_64 row, int_64 col) {
	std::vector<int_64> rows = {row};
	std::vector<int_64> cols = {col};
	std::vector<double> cell = cellFromRowCol(rows, cols);
	return  cell[0];
}

std::vector<double> SpatRaster::cellFromRowColCombine(std::vector<int_64> row, std::vector<int_64> col) {
	recycle(row, col);
	size_t n = row.size();
	int_64 nc = ncol();
	int_64 nr = nrow();

	std::vector<double> x(n * n);
	for (size_t i=0; i<n; i++) {
		for (size_t j=0; j<n; j++) {
			x[i*n+j] = (row[i]<0 || row[i] >= nr || col[j]<0 || col[j] >= nc) ? NAN : row[i] * nc + col[j];
		}
	}
	// duplicates occur if recycling occurs
	// could be avoided by smarter combination
	x.erase(std::remove_if(x.begin(), x.end(),
            [](const double& value) { return std::isnan(value); }), x.end());
	
	std::sort(x.begin(), x.end());
	x.erase(std::unique(x.begin(), x.end()), x.end());
	return x;
}


double SpatRaster::cellFromRowColCombine(int_64 row, int_64 col) {
	return cellFromRowCol(row, col);
}


std::vector<double> SpatRaster::yFromRow(std::vector<int_64> &row) {
	size_t size = row.size();
	std::vector<double> result( size );
	SpatExtent extent = getExtent();
	double ymax = extent.ymax;
	double yr = yres();
	int_64 nr = nrow();
	
	for (size_t i = 0; i < size; i++) {
		result[i] = (row[i] < 0 || row[i] >= nr ) ? NAN : ymax - ((row[i]+0.5) * yr);
	}
	return result;
}

double SpatRaster::yFromRow (int_64 row) {
	std::vector<int_64> rows = {row};
	std::vector<double> y = yFromRow(rows);
	return y[0];
}



std::vector<double> SpatRaster::xFromCol(std::vector<int_64> &col) {
	size_t size = col.size();
	std::vector<double> result( size );
	SpatExtent extent = getExtent();	
	double xmin = extent.xmin;
	double xr = xres();
	int_64 nc = ncol();
	for (size_t i = 0; i < size; i++) {
		result[i] = (col[i] < 0 || col[i] >= nc ) ? NAN : xmin + ((col[i]+0.5) * xr);
	}
	return result;
}

double SpatRaster::xFromCol(int_64 col) {
	std::vector<int_64> cols = {col};
	std::vector<double> x = xFromCol(cols);
	return x[0];
}

std::vector<int_64> SpatRaster::colFromX(std::vector<double> &x) {

	SpatExtent extent = getExtent();

	double xmin = extent.xmin;
	double xmax = extent.xmax;
	double xr = xres();
	size_t xs = x.size();
	std::vector<int_64> result(xs, -1);

	for (size_t i = 0; i < xs; i++) {
		if (x[i] >= xmin && x[i] < xmax ) {
			result[i] =  trunc((x[i] - xmin) / xr);
		} else if (x[i] == xmax) {
			result[i] = ncol()-1;
		}
	}
	return result;
}


int_64 SpatRaster::colFromX(double x) {
	std::vector<double> xv = {x};
	return colFromX(xv)[0];
}


std::vector<int_64> SpatRaster::rowFromY(std::vector<double> &y) {

	SpatExtent extent = getExtent();
	double ymin = extent.ymin;
	double ymax = extent.ymax;
	double yr = yres();
	size_t ys = y.size();
	std::vector<int_64> result(ys, -1);

	for (size_t i = 0; i < ys; i++) {
		if (y[i] > ymin && y[i] <= ymax) {
			result[i] = trunc((ymax - y[i]) / yr);
		} else if (y[i] == ymin) {
			result[i] = nrow() - 1;
		}		
	}
	return result;
}


int_64 SpatRaster::rowFromY(double y) {
	std::vector<double> Y = {y};
	return rowFromY(Y)[0];
}


std::vector<std::vector<double>> SpatRaster::xyFromCell( std::vector<double> &cell) {
	size_t n = cell.size();
	SpatExtent extent = getExtent();
	
	double xmin = extent.xmin;
	double ymax = extent.ymax;
	double yr = yres();
	double xr = xres();
    double ncells = ncell();
    long nc = ncol();
	std::vector< std::vector<double> > out(2, std::vector<double> (n, NAN) );
	for (size_t i = 0; i<n; i++) {
		if (std::isnan(cell[i]) || (cell[i] < 0) || (cell[i] >= ncells)) continue;
        long row = cell[i] / nc;
        long col = cell[i] - (row * nc);
        out[0][i] = xmin + (col + 0.5) * xr;
        out[1][i] = ymax - (row + 0.5) * yr;
	}
	return out;
}


std::vector< std::vector<double>> SpatRaster::xyFromCell( double cell) {
	std::vector<double> vcell = {cell};
	return xyFromCell(vcell);
}


std::vector<std::vector<int_64>> SpatRaster::rowColFromCell(std::vector<double> &cell) {
	size_t cs = cell.size();
	std::vector<std::vector<int_64>> result(2, std::vector<int_64> (cs, -1) );
	double nc = ncell();
	for (size_t i = 0; i < cs; i++) {
		if ((cell[i] >= 0) && (cell[i] < nc )) {
			result[0][i] = trunc(cell[i]/ ncol());
			result[1][i] = (cell[i] - ((result[0][i]) * ncol()));
		}
	}
	return result;
}


std::vector<std::vector<int_64>>  SpatRaster::rowColFromExtent(SpatExtent e) {
	std::vector<std::vector<double>> xy = e.asPoints();
	std::vector<int_64> col = colFromX(xy[0]); 
	std::vector<int_64> row = rowFromY(xy[1]); 
	std::vector<std::vector<int_64>> out = { row, col };
	return out;
}



std::vector<std::vector<double>> SpatRaster::adjacent(std::vector<double> cells, std::string directions, bool include) {

	unsigned n = cells.size();
	std::vector<std::vector<double>> out(n);

	std::vector<std::string> f {"rook", "queen", "bishop", "16"};
	if (std::find(f.begin(), f.end(), directions) == f.end()) {
        setError("argument directions is not valid");
        return(out);
	}

	std::vector<std::vector<int_64>> rc = rowColFromCell(cells);
	std::vector<int_64> r = rc[0];
	std::vector<int_64> c = rc[1];
	bool globlatlon = is_global_lonlat();
    int_64 nc = ncol();
    int_64 lc = nc-1;
    std::vector<int_64> cols, rows;
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

