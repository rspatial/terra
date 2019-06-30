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

// Based on  public-domain code by Darel Rex Finley, 2007
// http://alienryderflex.com/polygon_fill/

#include <vector>
#include "spatRaster.h"


std::vector<double> rasterize_polygon(std::vector<double> r, double value, const std::vector<double> &pX, const std::vector<double> &pY, const unsigned nrows, const unsigned ncols, const double xmin, const double ymax, const double rx, const double ry) {

	unsigned n = pX.size();
	std::vector<unsigned> nCol(n);
	for (size_t row=0; row<nrows; row++) {
		double y = ymax - (row+0.5) * ry;

		// find nodes.
		unsigned nodes = 0;
		size_t j = n-1;
		for (size_t i=0; i<n; i++) {
			if (((pY[i] < y) && (pY[j] >= y)) || ((pY[j] < y) && (pY[i] >= y))) {
			//	nCol[nodes++]=(int)  (((pX[i] - xmin + (y-pY[i])/(pY[j]-pY[i]) * (pX[j]-pX[i])) + 0.5 * rx ) / rx);
				double nds = ((pX[i] - xmin + (y-pY[i])/(pY[j]-pY[i]) * (pX[j]-pX[i])) + 0.5 * rx ) / rx;
				nds = nds < 0 ? 0 : nds;
		        nds = nds > ncols ? ncols : nds;
				nCol[nodes] = (unsigned) nds;
				nodes++;
			}
			j = i;
		}

		std::sort(nCol.begin(), nCol.begin()+nodes);
		unsigned ncell = ncols * row;

		//  Fill the cells between node pairs.
		for (size_t i=0; i < nodes; i+=2) {
			if (nCol[i+1] > 0 && nCol[i] < ncols) {
				for (size_t col = nCol[i]; col < nCol[i+1]; col++) {
					r[col + ncell] = value;
				}
			}
		}
	}
	return(r);
}




SpatRaster rasterizePolygons(SpatVector p, SpatRaster r, double background, SpatOptions &opt) {

	SpatRaster out = r.geometry(1);

  	if (!out.writeStart(opt)) { return out; }
	double value = 0;
	double resx = out.xres();
	double resy = out.yres();
	SpatGeom poly;
	SpatPart part;
	SpatHole hole;
	unsigned n = p.size();
	unsigned nc = out.ncol();
	SpatExtent extent = out.getExtent();
	
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> v(out.bs.nrows[i] * nc, background);

		for (size_t j = 0; j < n; j++) {
			poly = p.getGeom(j);
			//value = p.getAtt(j);
			//if (std::isnan(value)) { value = j;}
			value = j+1;
			unsigned np = poly.size();

			for (size_t k = 0; k < np; k++) {
				part = poly.getPart(k);
				if (part.hasHoles()) {
					std::vector<double> vv = rasterize_polygon(v, value, part.x, part.y, out.nrow(), out.ncol(), extent.xmin, extent.ymax, resx, resy);
					for (size_t h=0; h < part.nHoles(); h++) {
						hole = part.getHole(h);
						vv = rasterize_polygon(vv, background, hole.x, hole.y, out.nrow(), out.ncol(), extent.xmin, extent.ymax, resx, resy);
					}
					for (size_t q=0; q < vv.size(); q++) {
						if (vv[q] != background) {
							v[q] = vv[q];
						}
					}
				} else {
					v = rasterize_polygon(v, value, part.x, part.y, out.nrow(), out.ncol(), extent.xmin, extent.ymax, resx, resy);
				}
			}
		}
		if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i], 0, out.ncol())) return out;
		
	}
	out.writeStop();
	return(out);
}




std::vector<double> rasterize_line(std::vector<double> r, double value, const std::vector<double> &pX, const std::vector<double> &pY, const unsigned nrows, const unsigned ncols, const double xmin, const double ymax, const double rx, const double ry) {
	unsigned n = pX.size();
	for (size_t row=0; row<nrows; row++) {
		double y = ymax - (row+0.5) * ry;
		unsigned ncell = ncols * row;
		for (size_t i=1; i<n; i++) {
            size_t j = i-1;
			if (((pY[i] < y) && (pY[j] >= y)) || ((pY[j] < y) && (pY[i] >= y))) {
				double col = ((pX[i] - xmin + (y-pY[i])/(pY[j]-pY[i]) * (pX[j]-pX[i])) + 0.5 * rx ) / rx;
				if ((col >= 0) & (col < ncols)) {
                    r[ncell + col] = value;
				}
			}
		}
	}
	return(r);
}



SpatRaster rasterizeLines(SpatVector p, SpatRaster r, double background, SpatOptions &opt) {

	SpatRaster out = r.geometry(1);
  	if (!out.writeStart(opt)) { return out; }
	double value = 0;
	double resx = out.xres();
	double resy = out.yres();
	SpatGeom line;
	SpatPart part;
	unsigned n = p.size();
	SpatExtent extent = out.getExtent();

	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> v(out.bs.nrows[i] * out.ncol(), background);
		for (size_t j = 0; j < n; j++) {
			line = p.getGeom(j);
			//value = p.getAtt(j);
			//if (std::isnan(value)) { value = j;}
			value = j+1;
			unsigned nln = line.size();
			for (size_t k = 0; k < nln; k++) {
				part = line.getPart(k);
				v = rasterize_line(v, value, part.x, part.y, out.nrow(), out.ncol(), extent.xmin, extent.ymax, resx, resy);
			}
		}
		if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i], 0, out.ncol())) return out;
		
	}
	out.writeStop();
	return(out);
}


SpatRaster rasterizePoints(SpatVector p, SpatRaster r, double background, SpatOptions &opt) {
	r.setError("not implented yet");
	return(r);
}



SpatRaster SpatRaster::rasterize(SpatVector p, double background, SpatOptions &opt) {
	std::string gtype = p.type();
	SpatRaster out = geometry(1);
	if (gtype == "polygons") {
		out = rasterizePolygons(p, out, background, opt);
	} else if (gtype == "lines") {
		out = rasterizeLines(p, out, background, opt);		
	}  else {
		out = rasterizePoints(p, out, background, opt);		
	}
	return out;
}
