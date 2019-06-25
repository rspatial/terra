// Copyright (c) 2018-2019  Robert J. Hijmans
//
// This file is part of the "spat" library
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
#include "distance.h"
#include "vecmath.h"


double bilinear(const std::vector<double> &v, const  std::vector<double> &e, const double &dxdy, const double &x, const double &y) {
	// values
	// v[0] v[1]
	// v[2] v[3]

	// coordinates
	//           e[3] (ymax)
	// (xmin)e[0]  e[1] (xmax)
	//           e[2] (ymin)

    double dx1 = x - e[0];
    double dx2 = e[1] - x;
    double dy1 = y - e[2];
    double dy2 = e[3] - y;
    return (v[2] * dx2 * dy2 + v[3] * dx1 * dy2 + v[0] * dx2 * dy1 + v[1] * dx1 * dy1) / dxdy;
}



std::vector<double> SpatRaster::line_cells(SpatGeom& g) {

	unsigned nrows = nrow();
	unsigned ncols = ncol();
	double xmin = extent.xmin;
	double ymax = extent.ymax;
	double rx = xres();
	double ry = yres();
	std::vector<double> out;

	unsigned np = g.size();
	for (size_t prt=0; prt<np; prt++) {
        SpatPart p = g.getPart(prt);
        double miny = vmin(p.y, true);
        double maxy = vmax(p.y, true);

        double minrow = rowFromY(miny);
        double maxrow = rowFromY(maxy);
        if (minrow > nrows || maxrow < 0) {
            return(out);
        }
        size_t startrow = minrow < 0 ? 0 : minrow;
        size_t endrow = maxrow >= nrows ? (nrows-1) : maxrow;
        unsigned n = p.x.size();
        out.reserve(2*(startrow-endrow+1));

        for (size_t row=startrow; row<endrow; row++) {
            double y = ymax - (row+0.5) * ry;
            unsigned rowcell = ncols * row;
            for (size_t i=1; i<n; i++) {
                size_t j = i-1;
                if (((p.y[i] < y) && (p.y[j] >= y)) || ((p.y[j] < y) && (p.y[i] >= y))) {
                    double col = ((p.x[i] - xmin + (y-p.y[i])/(p.y[j]-p.y[i]) * (p.x[j]-p.x[i])) + 0.5 * rx ) / rx;
                    if ((col >= 0) & (col < ncols)) {
                        out.push_back(rowcell + col);
                    }
                }
            }
        }
	}
	return(out);
}




std::vector<double> SpatRaster::polygon_cells(SpatGeom& g) {

// does not deal with holes yet.

	unsigned nrows = nrow();
	unsigned ncols = ncol();
	double xmin = extent.xmin;
	double ymax = extent.ymax;
	double rx = xres();
	double ry = yres();
	std::vector<double> out;

	unsigned np = g.size();
	for (size_t prt=0; prt<np; prt++) {

        SpatPart p = g.getPart(prt);
        double miny = vmin(p.y, true);
        double maxy = vmax(p.y, true);
        double minrow = rowFromY(miny);
        double maxrow = rowFromY(maxy);
        if (minrow > nrows || maxrow < 0) {
            return(out);
        }
        size_t startrow = minrow < 0 ? 0 : minrow;
        size_t endrow = maxrow >= nrows ? (nrows-1) : maxrow;
        unsigned n = p.x.size();
        out.reserve(5*(startrow-endrow+1));

        std::vector<unsigned> nCol(n);
        for (size_t row=0; row<nrows; row++) {
            double y = ymax - (row+0.5) * ry;
            // find nodes.
            unsigned nodes = 0;
            size_t j = n-1;
            for (size_t i=0; i<n; i++) {
                if (((p.y[i] < y) && (p.y[j] >= y)) || ((p.y[j] < y) && (p.y[i] >= y))) {
                //	nCol[nodes++]=(int)  (((pX[i] - xmin + (y-pY[i])/(pY[j]-pY[i]) * (pX[j]-pX[i])) + 0.5 * rx ) / rx);
                    double nds = ((p.x[i] - xmin + (y-p.y[i])/(p.y[j]-p.y[i]) * (p.x[j]-p.x[i])) + 0.5 * rx ) / rx;
                    nds = nds < 0 ? 0 : nds;
                    nds = nds > ncols ? ncols : nds;
                    nCol[nodes] = (unsigned) nds;
                    nodes++;
                }
                j = i;
            }

            // now remove the holes?

            std::sort(nCol.begin(), nCol.begin()+nodes);
            unsigned rowcell = ncols * row;

            // fill  cells between node pairs.
            for (size_t i=0; i < nodes; i+=2) {
                if (nCol[i+1] > 0 && nCol[i] < ncols) { // surely should be >= 0?
                    for (size_t col = nCol[i]; col < nCol[i+1]; col++) {
                        out.push_back(col + rowcell);
                    }
                }
            }
        }
	}
	return(out);
}




// <layer<values>>
std::vector<std::vector<double>> SpatRaster::extractXY(std::vector<double> &x, std::vector<double> &y, std::string method) {

    unsigned nl = nlyr();
    unsigned np = x.size();
	std::vector<std::vector<double>> out(nl, std::vector<double>(np, NAN));
	if (!hasValues()) return out;

	if ((method == "idw") | (method == "bilinear")) {

// if nrow or nocl =1 disaggregate first

        bool lonlat = could_be_lonlat();
        bool globalLonLat = is_global_lonlat();
        size_t n = x.size();

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
            SpatRaster g = geometry();
			std::vector<unsigned> f = {2,2};
            SpatRaster gd = g.disaggregate(f, opt);

            double dyrs = gd.yres();
            double dxrs =  gd.xres();
            std::vector<double> d, cells(4);

            std::vector<std::vector<double> > cxy;

            std::vector<double> rc(4);
			unsigned nr = nrow();
			unsigned nc = ncol();
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
                    col1 = col1 < 0 ? mnc : col1 > mnc ? 0 : col1;
                    col2 = (cq == 0) ? col1-1 : col1 + 1;
                    col2 = col2 < 0 ? mnc : col2 > mnc ? 0 : col2;
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
                std::vector<std::vector<double>> xy = xyFromCell(cells);
                std::vector<std::vector<double>> v = extractCell(cells);
                std::vector<double> e = {xy[0][0], xy[0][1], xy[1][2], xy[1][0]};
                for (size_t j=0; j<nl; j++) {
                    out[j][i] = bilinear(v[j], e, dxdy, x[i], y[i]);
                }
            }
        }
	} else {

        std::vector<double> cell = cellFromXY(x, y);
        out = extractCell(cell);
	}

    return out;
}


// <geom<layer<values>>>
std::vector<std::vector<std::vector<double>>> SpatRaster::extractVector(SpatVector v) {

    unsigned nl = nlyr();
    unsigned ng = v.size();
    std::vector<std::vector<std::vector<double>>> out(ng, std::vector<std::vector<double>>(nl));

	if (!hasValues()) return out;

	std::vector<std::vector<double>> srcout;
	std::string gtype = v.type();
	if (gtype == "points") {
		SpatDataFrame vd = v.getGeometryDF();
		std::vector<double> x = vd.getD(0);
		std::vector<double> y = vd.getD(1);
		srcout = extractXY(x, y, "simple");
        for (size_t i=0; i<ng; i++) {
            for (size_t j=0; j<nl; j++) {
                out[i][j].push_back( srcout[j][i] );
            }
        }
	} else if (gtype == "lines") {
	    SpatRaster r = geometry(1);
	    SpatGeom g;
        for (size_t i=0; i<ng; i++) {
            g = v.getGeom(i);
            std::vector<double> cells = line_cells(g);
            srcout = extractCell(cells);
            for (size_t j=0; j<nl; j++) {
                out[i][j] = srcout[j];
            }
        }
	} else { // polys

	    SpatRaster r = geometry(1);
	    SpatRaster rc, rcr;
	    SpatVector p;
	    SpatGeom g;
	    SpatVector pts;
	    SpatOptions opt;
        for (size_t i=0; i<ng; i++) {
            g = v.getGeom(i);
            rc = r.crop(g.extent, "out", opt);
            p.setGeom(g);
            rcr = rc.rasterize(p, NAN, opt); // rather have a method that returns the cell numbers directly?
            pts = rcr.as_points(false, true);
            SpatDataFrame vd = pts.getGeometryDF();
            std::vector<double> x = vd.getD(0);
            std::vector<double> y = vd.getD(1);
            srcout = extractXY(x, y, "simple");
            //unsigned np = x.size();
            for (size_t j=0; j<nl; j++) {
               // unsigned off = j * np;
                out[i][j] = srcout[j];
//                std::copy(srcout.begin()+off, srcout.begin()+off+np-1, out[i][j].begin());
            }
        }
	}
	return out;
}



std::vector<std::vector<double>> SpatRaster::extractCell(std::vector<double> &cell) {

	unsigned n  = cell.size();
	unsigned nc = ncell();
	std::vector<std::vector<double>> out(nlyr(), std::vector<double>(n, NAN));
	if (!hasValues()) return out;
	
	unsigned ns = nsrc();

	unsigned lyr = 0;
	for (size_t src=0; src<ns; src++) {
		unsigned slyrs = source[src].layers.size();

		if (source[src].driver == "memory") {
			for (size_t i=0; i<slyrs; i++) {
				size_t j = i * nc;
				for (size_t k=0; k<n; k++) {
					if (!is_NA(cell[k]) && cell[k] >= 0 && cell[k] < nc) {
						out[lyr][k] = source[src].values[j + cell[k]];
						//out[offset + k] = source[src].values[j + cell[k]];
					}
				}
				lyr++;
			}
		} else {
			std::vector<std::vector<double>> srcout;
			if (source[0].driver == "raster") {
				srcout = readCellsBinary(src, cell); 
			} else {
			#ifdef useGDAL
				std::vector<std::vector<unsigned>> rc = rowColFromCell(cell);
				srcout = readRowColGDAL(src, rc[0], rc[1]); 
			#endif
			}
			for (size_t i=0; i<slyrs; i++) {
				out[lyr] = srcout[i];
				lyr++;
			}
		}
	}
	return out;
}


