// Copyright (c) 2018  Robert J. Hijmans
//
// This file is part of the "spat" library.v
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


double bilinear(std::vector<double> v, std::vector<double> e, double dxdy, double x, double y) {
    double dx1, dx2, dy1, dy2;
    dx1 = e[1] - x;
    dx2 = x - e[0];
    dy1 = y - e[2];
    dy2 = e[3] - y;
    return (v[2] * dx2 * dy2 + v[3] * dx1 * dy2 +  v[0] * dx2 * dy1 + v[1] * dx1 * dy1) / dxdy;
}


std::vector<double> SpatRaster::extractXY(std::vector<double> &x, std::vector<double> &y, std::string method) {

	std::vector<double> out;

	if ((method == "idw") | (method == "bilinear")) {

// if nrow or nocl =1 disaggregate first

        size_t nr = nrow();
        size_t nc = ncol();
        bool lonlat = could_be_lonlat();
        bool globalLonLat = is_global_lonlat();
        size_t n = x.size();
        out.resize(n, NAN);


        if (method == "idw") {

            std::function<std::vector<double>(std::vector<double>&,std::vector<double>&,double,double)> distFun;
 //           std::vector<double> distance_plane(std::vector<double> &x1, std::vector<double> &y1, std::vector<double> &x2, std::vector<double> &y2);
            if (lonlat) {
                distFun = distance_lonlat_vd;
            } else {
                distFun = distance_plane_vd;
            }

/*
                cxy = xyFromCell(cells);
                d = distFun(cxy[0], cxy[1], x[i], y[i]);
                v = extractCell(cells);

                double a=0, b=0;
                for (size_t j=0; j<4; j++) {
                    a += v[j] * d[j];
                    b += d[j];
                }
                out[i] = a / b;
*/
        } else if (method == "bilinear") {
            double ymax = extent.ymax;
            double xmin = extent.xmin;
            double yrs = yres();
            double xrs = xres();

            SpatOptions opt;
            SpatRaster g = geometry(1);
            SpatRaster gd = g.disaggregate(std::vector<unsigned>{2,2}, opt);
            double dyrs = gd.yres();
            double dxrs =  gd.xres();
            std::vector<double> v, d, cells(4);
            std::vector<std::vector<double> > cxy;

            std::vector<double> rc(4);
            unsigned mnr = nr-1;
            unsigned mnc = nc-1;

            // needs row-wise adjustment for lonlat
            double dxdy = xres() * yres();


            for (size_t i=0; i<n; i++) {
                long row_d = ((ymax - y[i]) / dyrs);
                long col_d = ((x[i] - xmin) / dxrs);
                unsigned rq = row_d % 2;
                unsigned cq = col_d % 2;


                double row1 = floor((ymax - y[i]) / yrs);
                double col1 = floor((x[i] - xmin) / xrs);
                if ((row1 < 0) || (row1 > mnr)) { continue; }

                double row2 = (rq == 0) ? row1-1 : row1+1;
                row2 = row2 < 0 ? row1+1 : row2==nr ? row1-1 : row2;
                double col2;
                  if (globalLonLat) {
                    if ((col1 < -1) | (col1 > nc)) { continue; }
                    col1 = col1 < 0 ? nc+col1 : col1 > mnc ? col1-nc : col1;
                    col2 = (cq == 0) ? col1-1 : col1 + 1;
                    col2 = col2 < 0 ? nc+col2 : col2 > mnc ? col2-nc : col2;
                } else {
                     if ((col1 < 0) | (col1 > mnc)) { continue; }
                     col2 = (cq == 0) ? col1-1 : col1 + 1;
                     col2 = col2 < 0 ? col1+1 : col2 == nc ? col1-1 : col2;
                }
                cells[0] = nc * row1 + col1;
                cells[1] = nc * row1 + col2;
                cells[2] = nc * row2 + col1;
                cells[3] = nc * row2 + col2;
                std::sort(cells.begin(), cells.end());
                std::vector<std::vector<double> > xy = xyFromCell(cells);
                v = extractCell(cells);
//                std::vector<double> e = extent.asVector();
                std::vector<double> e = {xy[0][0], xy[0][1], xy[1][2], xy[1][0]};
                out[i] = bilinear(v, e, dxdy, x[i], y[i]);
            }
        }
	} else {

        std::vector<double> srcout;

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



std::vector<std::vector<double> > SpatRaster::extractVector(SpatVector v, std::string fun) {

    std::vector<std::vector<double> > out;
	std::vector<double> srcout;

	std::string gtype = v.type();
	if (gtype == "points") {
		SpatDataFrame vd = v.getGeometryDF();
		std::vector<double> x = vd.getD(0);
		std::vector<double> y = vd.getD(1);
		out.resize(1);
		out[0] = extractXY(x, y, "simple");

	} else { //if (gtype == "lines" | polys) {
	    SpatRaster r = geometry(1);
	    SpatRaster x, y;
	    SpatVector p;
	    SpatGeom g;
	    SpatVector points;
	    size_t vs = v.size();
	    out.resize(vs);
	    SpatOptions opt;
        for (size_t i=0; i<vs; i++) {
            g = v.getGeom(i);
            x = r.crop(g.extent, "out", opt);
            p.setGeom(g);
            y = x.rasterizePolygons(p, NAN, opt);
            points = y.as_points(false, true);
            std::vector<std::vector<double> > vp = extractVector(points);
            out[i] = vp[0];
        }
	}
	return out;
}



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

