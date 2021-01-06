// Copyright (c) 2018-2020  Robert J. Hijmans
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

#include "spatRasterMultiple.h"
#include "distance.h"
#include "vecmath.h"



//#include <iostream>

double rowColToCell(unsigned ncols, unsigned row, unsigned col) {
  return row * ncols + col;
}


double getRow(const unsigned& nrow, const double& y, const double& ymin, const double& ymax, const double& invyr) {
  double result = NAN;
  if (y == ymin) {
    result = nrow-0.0000001;
  } else if ((y > ymin) && (y <= ymax)) {
    result = (ymax - y) * invyr;
  }
  return result;
}


double getCol(const unsigned& ncol, const double& x, const double& xmin, const double& xmax, const double& invxr) {
    double result = NAN;
    if (x == xmax) {
      result = ncol-0.0000001;
    } else if ((x >= xmin) & (x < xmax)) {
      result = (x - xmin) * invxr;
    }
  return result;
}


std::vector<double> fourCellsFromXY (unsigned ncols, unsigned nrows, double xmin, double xmax, double ymin, double ymax,
                                     std::vector<double> x, std::vector<double> y, bool duplicates, bool isGlobalLonLat ) {

  size_t n = x.size();
  double yres_inv = nrows / (ymax - ymin);
  double xres_inv = ncols / (xmax - xmin);
  std::vector<double> result(4*n, NAN);
  unsigned maxrow = nrows-1;
  unsigned maxcol = ncols-1;

  for (size_t i = 0; i < n; i++) {

    double row = getRow(nrows, y[i], ymin, ymax, yres_inv);
    double col = getCol(ncols, x[i], xmin, xmax, xres_inv);
//    std::cout << row << " " << col << "\n";
    if (std::isnan(row) || std::isnan(col)) {
      continue;
    }
    //double row = (ymax - y[i]) * yres_inv;
    //double col = (x[i] - xmin) * xres_inv;

    double roundRow = std::floor(row);
    double roundCol = std::floor(col);

    if (roundRow < 0 || roundRow > maxrow || roundCol < 0 || roundCol > maxcol) {
      continue;
    }

    // roundRow and roundCol are now the nearest row/col to x/y.
    // That gives us one corner. We will find the other corner by starting
    // at roundRow/roundCol and moving in the direction of row/col, stopping
    // at the next integral values.

    // >0 if row is greater than the nearest round row, 0 if equal
    double vertDir = row - roundRow  - 0.5;
    // >0 if col is greater than the nearest round col, 0 if equal
    double horizDir = col - roundCol - 0.5;

    // If duplicates are not allowed, make sure vertDir and horizDir
    // are not 0
    if (!duplicates) {
      if (vertDir == 0)
        vertDir = 1;
      if (horizDir == 0)
        horizDir = 1;
    }

    // roundRow and roundCol will be one corner; posRow and posCol will be
    // the other corner. Start out by moving left/right or up/down relative
    // to roundRow/roundCol.
    double posRow = roundRow + (vertDir > 0 ? 1 : vertDir < 0 ? -1 : 0);
    double posCol = roundCol + (horizDir > 0 ? 1 : horizDir < 0 ? -1 : 0);

    // Now, some fixups in case posCol/posRow go off the edge of the raster.
    if (isGlobalLonLat) {
      if (posCol < 0) {
        posCol = maxcol;
      } else if (posCol > maxcol) {
        posCol = 0;
      }
    } else {
      if (posCol < 0) {
        posCol = 1;
      } else if (posCol > maxcol) {
        posCol = maxcol - 1;
      }
    }

    if (posRow < 0) {
      posRow = 1;
    } else if (posRow > maxrow) {
      posRow = maxrow - 1;
    }

    // Fixups done--just store the results.

    std::vector<double> res(4);
    res[0] = rowColToCell(ncols, roundRow, roundCol);
    res[1] = rowColToCell(ncols, posRow, roundCol);
    res[2] = rowColToCell(ncols, posRow, posCol);
    res[3] = rowColToCell(ncols, roundRow, posCol);
    std::sort(res.begin(), res.end());

    size_t j = i * 4;
    result[j]   = res[0];
    result[j+1] = res[1];
    result[j+2] = res[2];
    result[j+3] = res[3];
  }
  return result;
}


double linearInt(const double& d, const double& x, const double& x1, const double& x2, const double& v1, const double& v2) {
  double result = (v2 * (x - x1) + v1 * (x2 - x)) / d;
  return result;
}


double bilinearInt(const double& x, const double& y, const double& x1, const double& x2, const double& y1, const double& y2, const double& v11, const double& v21, const double& v12, const double& v22) {
  double d = x2-x1;
  double h1 =  linearInt(d, x, x1, x2, v11, v21);
  double h2 =  linearInt(d, x, x1, x2, v12, v22);
  d = y2-y1;
  double v =  linearInt(d, y, y1, y2, h1, h2);
  return v;
}

double distInt(double d, double pd1, double pd2, double v1, double v2) {
  double result = (v2 * pd1 + v1 * pd2) / d;
  return result;
}

double bilinear_geo(double x, double y, double x1, double x2, double y1, double y2, double halfyres, std::vector<double> vv) {

    double a = 6378137.0;
    double f = 1/298.257223563;
	double hy = y1 - halfyres;
	double d = distance_lonlat(x1, hy, x2, hy, a, f);

    std::vector<double> dist(4);
	double pd1 = distance_lonlat(x, hy, x1, hy, a, f);
	double pd2 = distance_lonlat(x, hy, x2, hy, a, f);
	double h1 = distInt(d, pd1, pd2, vv[0], vv[1]);
	double h2 = distInt(d, pd1, pd2, vv[2], vv[3]);
	d = y2 - y1;
	double v =  linearInt(d, y, y1, y2, h1, h2);

	return v;
}



std::vector<std::vector<double>> SpatRaster::bilinearValues(std::vector<double> x, std::vector<double> y) {

    bool glob = is_global_lonlat();
//    bool lonlat = could_be_lonlat();
	SpatExtent extent = getExtent();
	std::vector<double> four = fourCellsFromXY(ncol(), nrow(), extent.xmin, extent.xmax, extent.ymin, extent.ymax, x, y, false, glob);
	std::vector<std::vector<double>> xy = xyFromCell(four);
	std::vector<std::vector<double>> v = extractCell(four);
	size_t n = x.size();
//	double halfyres = yres()/2;
	std::vector<std::vector<double>> res(nlyr(), std::vector<double>(n));
/*	if (lonlat) {
        for (size_t i=0; i<n; i++) {
            size_t ii = i * 4;
            for (size_t j=0; j<nlyr(); j++) {
                std::vector<double> vv(v[j].begin()+ii, v[j].begin()+ii+4);
        if (i==0) {
            std::cout << x[i] << " "<< y[i] << "\n";
            std::cout << xy[0][ii] << " " << xy[0][ii+1] << " " << xy[1][ii] << " " << xy[1][ii+3] << "\n";
            std::cout << v[j][ii] << " " << v[j][ii+1] << " " << v[j][ii+2] << " " << v[j][ii+3]  << "\n";
        }
                res[j][i] = bilinear_geo(x[i], y[i], xy[0][ii], xy[0][ii+1], xy[1][ii], xy[1][ii+3], halfyres, vv);

            }
        }
	} else {
 */       for (size_t i=0; i<n; i++) {
            size_t ii = i * 4;
            for (size_t j=0; j<nlyr(); j++) {

  /*      if (i==0) {
            std::cout << x[i] << " "<< y[i] << "\n";
            std::cout << xy[0][ii] << " " << xy[0][ii+1] << " " << xy[1][ii] << " " << xy[1][ii+3] << "\n";
            std::cout << v[j][ii] << " " << v[j][ii+1] << " " << v[j][ii+2] << " " << v[j][ii+3]  << "\n";
        }
*/
                res[j][i] = bilinearInt(x[i], y[i], xy[0][ii], xy[0][ii+1], xy[1][ii], xy[1][ii+3], v[j][ii], v[j][ii+1], v[j][ii+2], v[j][ii+3]);
            }
        }
//	}
	return res;
}




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
	SpatExtent extent = getExtent();

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
	SpatExtent extent = getExtent();

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


/*
idw
        bool lonlat = could_be_lonlat();
        //bool globalLonLat = is_global_lonlat();
        //size_t n = x.size();

        if (method == "idw") {
            std::function<std::vector<double>(std::vector<double>&,std::vector<double>&,double,double)> distFun;
 //           std::vector<double> distance_plane(std::vector<double> &x1, std::vector<double> &y1, std::vector<double> &x2, std::vector<double> &y2);
            if (lonlat) {
                distFun = distance_lonlat_vd;
            } else {
                distFun = distance_plane_vd;
            }
*/
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


// <layer<values>>
std::vector<std::vector<double>> SpatRaster::extractXY(std::vector<double> &x, std::vector<double> &y, std::string method, bool cells) {

    unsigned nl = nlyr();
    unsigned np = x.size();
	if (!hasValues()) {
		std::vector<std::vector<double>> out(nl+cells, std::vector<double>(np, NAN));
		return out;
	}
	std::vector<std::vector<double>> out;

    if (method == "bilinear") {
			out = bilinearValues(x, y);

	} else {

        std::vector<double> cell = cellFromXY(x, y);
        out = extractCell(cell);
		if (cells) {
			out.push_back(cell);
		}
	}

    return out;
}


// <geom<layer<values>>>
std::vector<std::vector<std::vector<double>>> SpatRaster::extractVector(SpatVector v, bool touches, std::string method, bool cells, bool weights) {

	std::string gtype = v.type();
	if (gtype != "polygons") weights = false;

    unsigned nl = nlyr();
    unsigned ng = v.size();
    std::vector<std::vector<std::vector<double>>> out(ng, std::vector<std::vector<double>>(nl + cells + weights));

	if (!hasValues()) return out;
	#if GDAL_VERSION_MAJOR < 3
	if (weights) {
		setError("extract with weights not supported for your GDAL version");
		return out;
	}
	#endif

	std::vector<std::vector<double>> srcout;
	if (gtype == "points") {
		if (method != "bilinear") method = "simple";
		SpatDataFrame vd = v.getGeometryDF();
		if (vd.nrow() == ng) {  // single point geometry
			std::vector<double> x = vd.getD(0);
			std::vector<double> y = vd.getD(1);
			srcout = extractXY(x, y, method, cells);
			for (size_t i=0; i<ng; i++) {
				for (size_t j=0; j<nl; j++) {
					out[i][j].push_back( srcout[j][i] );
				}
				if (cells) {
					out[i][nl].push_back( srcout[nl][i] );			
				}
			}
		} else { // multipoint
			for (size_t i=0; i<ng; i++) {
				SpatVector vv = v.subset_rows(i);
				SpatDataFrame vd = vv.getGeometryDF();
				std::vector<double> x = vd.getD(0);
				std::vector<double> y = vd.getD(1);
				srcout = extractXY(x, y, method, cells);
				for (size_t j=0; j<nl; j++) {
					out[i][j] = srcout[j];
				}
				if (cells) {
					out[i][nl] = srcout[nl];			
				}
			}
		}
	} else {
	    SpatRaster r = geometry(1);
	    //SpatOptions opt;
		//std::vector<double> feats(1, 1) ;		
        for (size_t i=0; i<ng; i++) {
            SpatGeom g = v.getGeom(i);
            SpatVector p(g);
			p.srs = v.srs;
			std::vector<double> cell, wgt;
			if (weights) {
				std::vector<std::vector<double>> cw = rasterizeCellsWeights(p, touches);
				cell = cw[0];
				wgt = cw[1];
			} else {
				cell = rasterizeCells(p, touches);
            }
			srcout = extractCell(cell);
            for (size_t j=0; j<nl; j++) {
                out[i][j] = srcout[j];
            }
			if (cells) {
				out[i][nl] = cell;
			}
			if (weights) {
				out[i][nl+cells] = wgt;
			}
        }
	}
	return out;
}



std::vector<std::vector<double>> SpatRaster::extractCell(std::vector<double> &cell) {

	std::vector<double> wcell;
	std::vector<std::vector<int_64>> rc, wrc;
	rc = rowColFromCell(cell);

	size_t n  = cell.size();
	std::vector<std::vector<double>> out(nlyr(), std::vector<double>(n, NAN));
	if (!hasValues()) return out;

	unsigned ns = nsrc();
	unsigned lyr = 0;
	size_t nc;
	for (size_t src=0; src<ns; src++) {
		unsigned slyrs = source[src].layers.size();
		bool win = source[src].hasWindow;
		if (win) {
			nc = source[src].window.full_ncol * source[src].window.full_nrow;
			wrc = rc;
			wcell.reserve(cell.size());
			for (size_t i=0; i<cell.size(); i++) {
				if ((wrc[0][i] < 0) || (wrc[1][i] <0)) {
					wcell.push_back( NAN );				
				} else {
					wrc[0][i] += source[src].window.off_row;
					wrc[1][i] += source[src].window.off_col;
					wcell.push_back( wrc[0][i] * source[src].window.full_ncol + wrc[1][i] );
				}
			}
		} else {
			nc = ncell();
		}
		if (source[src].memory) {
			for (size_t i=0; i<slyrs; i++) {
				size_t j = i * nc;
				if (win) {
					for (size_t k=0; k<n; k++) {
						if (!is_NA(wcell[k]) && wcell[k] >= 0 && wcell[k] < nc) {
							out[lyr][k] = source[src].values[j + wcell[k]];
						}
					}				
				} else {
					for (size_t k=0; k<n; k++) {
						if (!is_NA(cell[k]) && cell[k] >= 0 && cell[k] < nc) {
							out[lyr][k] = source[src].values[j + cell[k]];
						}
					}
				}
				lyr++;
			}
		} else {
			std::vector<std::vector<double>> srcout;
			//if (source[0].driver == "raster") {
			//	srcout = readCellsBinary(src, cell);
			//} else {
			#ifdef useGDAL
			if (win) {
				srcout = readRowColGDAL(src, wrc[0], wrc[1]);
			} else {
				srcout = readRowColGDAL(src, rc[0], rc[1]);			
			}
			#endif
			if (hasError()) return out;
			//}
			for (size_t i=0; i<slyrs; i++) {
				out[lyr] = srcout[i];
				lyr++;
			}
		}
	}
	return out;
}



std::vector<double> SpatRaster::vectCells(SpatVector v, bool touches, std::string method, bool weights) {

	std::string gtype = v.type();
	if (gtype != "polygons") weights = false;
	std::vector<double> out, cells, wghts;
	if (gtype == "points") {
		if (method != "bilinear") method = "simple";
		SpatDataFrame vd = v.getGeometryDF();
		cells = cellFromXY(vd.getD(0), vd.getD(1));
		std::vector<long> id = vd.getI(0);
		out.insert(out.end(), id.begin(), id.end());
		out.insert(out.end(), cells.begin(), cells.end());
	} else {
		unsigned ng = v.size();
	    SpatOptions opt;
		SpatRaster r = geometry(1);
		std::vector<double> feats(1, 1) ;		
        for (size_t i=0; i<ng; i++) {
            SpatGeom g = v.getGeom(i);
            SpatVector p(g);
			p.srs = v.srs;
			if (weights) {
				std::vector<std::vector<double>> cw = rasterizeCellsWeights(p, touches);			
				std::vector<double> id(cw[0].size(), i);
				out.insert(out.end(), id.begin(), id.end());
				cells.insert(cells.end(), cw[0].begin(), cw[0].end());			
				wghts.insert(wghts.end(), cw[1].begin(), cw[1].end());			
			} else {
				std::vector<double> geomc = rasterizeCells(p, touches);
				std::vector<double> id(geomc.size(), i);
				out.insert(out.end(), id.begin(), id.end());
				cells.insert(cells.end(), geomc.begin(), geomc.end());
			}
        }
		if (weights) {
			out.insert(out.end(), cells.begin(), cells.end());
			out.insert(out.end(), wghts.begin(), wghts.end());
		} else {
			out.insert(out.end(), cells.begin(), cells.end());
		}
	}
	return out;
}


std::vector<double> SpatRaster::extCells(SpatExtent ext) {

	std::vector<double> out;
	ext = align(ext, "near");
	ext.intersect(getExtent());
	if (!ext.valid()) {
		return(out);
	}
	double resx = xres() / 2;
	double resy = yres() / 2;
	std::vector<double> e = ext.asVector();
	e[0] += resx;
	e[1] -= resx;
	e[2] += resy;
	e[3] -= resy;
	std::vector<double> ex = {e[0], e[1]};
	std::vector<double> ey = {e[3], e[2]};
	std::vector<int_64> r = rowFromY(ey);
	std::vector<int_64> c = colFromX(ex);
	int_64 nc = ncol();
	out.reserve((r[1]-r[0]) * (c[1]-c[0]));
	for (int_64 i=r[0]; i <= r[1]; i++) {
		for (int_64 j=c[0]; j <= c[1]; j++) {
			out.push_back(i*nc+j); 
		}
	}
	return out;
}



std::vector<std::vector<std::vector<double>>> SpatRasterStack::extractXY(std::vector<double> &x, std::vector<double> &y, std::string method) {
	unsigned ns = nsds();
	std::vector<std::vector<std::vector<double>>> out(ns);
	for (size_t i=0; i<ns; i++) {
		SpatRaster r = getsds(i);
		out[i] = r.extractXY(x, y, method);
	}
	return out;
}

std::vector<std::vector<std::vector<double>>> SpatRasterStack::extractCell(std::vector<double> &cell) {
	unsigned ns = nsds();
	std::vector<std::vector<std::vector<double>>> out(ns);
	for (size_t i=0; i<ns; i++) {
		SpatRaster r = getsds(i);
		out[i] = r.extractCell(cell);
	}
	return out;
}


// this is rather inefficient (repeated rasterization)
std::vector<std::vector<std::vector<std::vector<double>>>> SpatRasterStack::extractVector(SpatVector v, bool touches, std::string method) {
	unsigned ns = nsds();
	std::vector<std::vector<std::vector<std::vector<double>>>> out(ns);
	for (size_t i=0; i<ns; i++) {
		SpatRaster r = getsds(i);
		out[i] = r.extractVector(v, touches, method);
	}
	return out;
}



/*
        } else if (method == "oldbilinear") {

// this is much too slow
			SpatExtent extent = getExtent();

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

                double row1 = std::floor((ymax - y[i]) / yrs);
                double col1 = std::floor((x[i] - xmin) / xrs);
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
		*/


/*
std::vector<double> SpatRaster::extractCell(std::vector<double> &cell) {

	unsigned n  = cell.size();
	unsigned nc = ncell();

	std::vector<double> out;
	if (!hasValues()) {
		out = std::vector<double>(n * nlyr(), NAN)
		return out;
	}

	unsigned ns = nsrc();
	for (size_t src=0; src<ns; src++) {
		unsigned slyrs = source[src].layers.size();
		std::vector<double> srcout;
		if (source[src].memory) {
			srcout = std::vector<double>(n * slyrs, NAN)
			std::vector<size_t> off1(slyrs);
			std::vector<size_t> off2(slyrs);
			for (size_t i=0; i<slyrs; i++) {
				off1[i] = i * n;
				off2[i] = i * nc;
			}
			for (size_t i=0; i<n; i++) {
				if (!is_NA(cell[i]) && cell[i] >= 0 && cell[i] < nc) {
					for (size_t j=0; j<slyrs; j++) {
						out[off1[j]+i] = source[src].values[off2[j] + cell[i]];
					}
				}
			}
		} else {
			#ifdef useGDAL
			std::vector<std::vector<int_64>> rc = rowColFromCell(cell);
			srcout = readRowColGDAL(src, rc[0], rc[1]);
			#endif
			if (hasError()) return out;
			//}
		}
		out.insert(out.end(), srcout.begin(), srcout.end());	
	}
	return out;
}
*/
