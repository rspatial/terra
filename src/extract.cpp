// Copyright (c) 2018  Robert J. Hijmans
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

#include <functional>

#include "spatRaster.h"
#include "NA.h"
#include "distance.h"

std::vector<double> SpatRaster::extractCell(std::vector<double> &cell) {

	unsigned n = cell.size();
	unsigned nc = ncell();
	std::vector<double> out(n * nlyr(), NAN);
	unsigned ns = nsrc();
	unsigned offset = 0;

	std::vector<std::vector<unsigned> > rc = rowColFromCell(cell);
	std::vector<unsigned> rows = rc[0];
	std::vector<unsigned> cols = rc[1];

	for (size_t src=0; src<ns; src++) {

		unsigned slyrs = source[src].nlyr;
		if (source[src].driver == "memory") {
			for (size_t i=0; i<slyrs; i++) {
				size_t j = i * nc;
				for (size_t k=0; k<n; k++) {
					if (!is_NA(cell[k]) && cell[k] >= 0 && cell[k] < nc) {
						out[offset + k] = source[src].values[j + cell[k]];
					}
				}
				offset += n;
			}
		} else {
		#ifdef useGDAL
			std::vector<double> srcout = readRowColGDAL(src, rows, cols); //
			std::copy(srcout.begin(), srcout.end(), out.begin()+offset);
			offset += n;
        #endif
		}
	}
	return out;
}


std::vector<double> SpatRaster::extractXY(std::vector<double> &x, std::vector<double> &y, std::string method) {
	std::vector<double> out, srcout;
	if (method == "bilinear") {

        size_t nr = nrow();
        size_t nc = ncol();
        bool lonlat = could_be_lonlat();

        std::vector<double> distance_plane(std::vector<double> &x1, std::vector<double> &y1, std::vector<double> &x2, std::vector<double> &y2);

        std::function<std::vector<double>(std::vector<double>&,std::vector<double>&,double,double)> distFun;
        if (lonlat) {
            distFun = distance_lonlat_vd;
        } else {
            distFun = distance_plane_vd;
        }

        bool isGlobalLonLat = extent.is_global_lonlat(getCRS());
        size_t n = x.size();
        double ymax = extent.ymax;
        double xmin = extent.xmin;

        double yres_inv = nr / (extent.ymax - extent.ymin);
        double xres_inv = nc / (extent.xmax - extent.xmin);

        out.resize(n, NAN);
        std::vector<double> v, d, cells(4);
        std::vector<std::vector<double> > cxy;

        for (size_t i=0; i<n; i++) {
            double row = (ymax - y[i]) * yres_inv - 0.5;
            double col = (x[i] - xmin) * xres_inv - 0.5;
            double roundRow = round(row);
            double roundCol = round(col);
            // Check for out-of-bounds.
            if (roundRow < 0 || roundRow > (nr-1) || roundCol < 0 || roundCol > (nc-1)) {
                continue;
            }
// roundRow and roundCol are now the nearest row/col to x/y. That gives us one corner. We will find the other corner by starting
// at roundRow/roundCol and moving in the direction of row/col, stopping at the next integral values.
// >0 if row is greater than the nearest round row, 0 if equal
            double vertDir = row - roundRow;
// >0 if col is greater than the nearest round col, 0 if equal
            double horizDir = col - roundCol;
// roundRow and roundCol will be one corner; posRow and posCol will be the other corner. Start out by moving left/right or up/down relative to roundRow/roundCol.
            double posRow = roundRow + (vertDir > 0 ? 1 : vertDir < 0 ? -1 : 0);
            double posCol = roundCol + (horizDir > 0 ? 1 : horizDir < 0 ? -1 : 0);
// Now, some fixups in case posCol/posRow go off the edge of the raster.
            if (isGlobalLonLat) {
                if (posCol < 0) {
                    posCol = nc-1;
                } else if (posCol > (nc-1)) {
                    posCol = 0;
                }
            } else {
                if (posCol < 0) {
                    posCol = 1;
                } else if (posCol > (nc-1)) {
                    posCol = nc - 2;
                }
            }
            if (posRow < 0) {
                posRow = 1;
            } else if (posRow > (nr-1)) {
                posRow = nr - 2;
            }

            cells[0] = nc * roundRow + roundCol;
            cells[1] = nc * posRow + roundCol;
            cells[2] = nc * posRow + posCol;
            cells[3] = nc * roundRow + posCol;
            cxy = xyFromCell(cells);
            d = distFun(cxy[0], cxy[1], x[i], y[i]);
            v = extractCell(cells);

            double a, b;
            for (size_t j=0; j<4; j++) {
                a += v[j] * d[j];
                b += d[j];
            }
            out[i] = a / b;
        }
	} else {
        for (size_t src=0; src<nsrc(); src++) {
            if (source[src].driver == "memory") {
               std::vector<double> cell = cellFromXY(x, y);
               srcout = extractCell(cell);
            } else {
               std::vector<unsigned> rows = rowFromY(y);
               std::vector<unsigned> cols = colFromX(x);
               #ifdef useGDAL
               srcout = readRowColGDAL(src, rows, cols);
               #endif
                std::vector<double> srcout(x.size());
            }
            out.insert(out.end(), srcout.begin(), srcout.end());
        }
	}
    return out;
}



std::vector<double> SpatRaster::extractVector(SpatVector v, std::string fun) {

	std::vector<double> out;
	std::vector<double> srcout;

	std::string gtype = v.type();
	if (gtype == "points") {
		SpatDataFrame vd = v.getGeometryDF();
		std::vector<double> x = vd.getD(0);
		std::vector<double> y = vd.getD(1);
		for (size_t src=0; src<nsrc(); src++) {
			if (source[src].driver == "memory") {
				std::vector<double> cell = cellFromXY(x, y);
				srcout = extractCell(cell);
			} else {
				std::vector<double> x = vd.getD(2);
				std::vector<double> y = vd.getD(3);
				std::vector<unsigned> rows = rowFromY(y);
				std::vector<unsigned> cols = colFromX(x);
				#ifdef useGDAL
				srcout = readRowColGDAL(src, rows, cols);
				#endif
			}
		}
		out.insert(out.end(), srcout.begin(), srcout.end());

	} else if (gtype == "lines") {


	} else {


	}
	return out;
}


