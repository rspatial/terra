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

// inspired by http://alienryderflex.com/polygon_fill/

#include <vector>
#include "spatRaster.h"


std::vector<double> rasterize_polygon(std::vector<double> r, double value, const std::vector<double> &pX, const std::vector<double> &pY, const unsigned startrow, const unsigned nrows, const unsigned ncols, const double xmin, const double ymax, const double rx, const double ry) {

	unsigned n = pX.size();
	std::vector<unsigned> nCol(n);
	for (size_t row=0; row < (nrows); row++) {
		double y = ymax - (startrow+row+0.5) * ry;

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




SpatRaster rasterizePolygons(SpatVector p, SpatRaster r, std::vector<double> value, double background, SpatOptions &opt) {

	SpatRaster out = r.geometry(1);

  	if (!out.writeStart(opt)) { return out; }
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
			unsigned np = poly.size();

			for (size_t k = 0; k < np; k++) {
				part = poly.getPart(k);
				if (part.hasHoles()) {
					std::vector<double> vv = rasterize_polygon(v, value[j], part.x, part.y, out.bs.row[i], out.bs.nrows[i], out.ncol(), extent.xmin, extent.ymax, resx, resy);
					for (size_t h=0; h < part.nHoles(); h++) {
						hole = part.getHole(h);
						vv = rasterize_polygon(vv, background, hole.x, hole.y, out.bs.row[i], out.bs.nrows[i], out.ncol(), extent.xmin, extent.ymax, resx, resy);
					}
					for (size_t q=0; q < vv.size(); q++) {
						if ((vv[q] != background) && (!std::isnan(vv[q]))) {
							v[q] = vv[q];
						}
					}
				} else {
					v = rasterize_polygon(v, value[j], part.x, part.y, out.bs.row[i], out.bs.nrows[i], out.ncol(), extent.xmin, extent.ymax, resx, resy);
				}
			}
		}
		if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i], 0, out.ncol())) return out;

	}
	out.writeStop();
	return(out);
}




std::vector<double> rasterize_line(std::vector<double> r, double value, const std::vector<double> &pX, const std::vector<double> &pY, const unsigned startrow, const unsigned nrows, const unsigned ncols, const double xmin, const double ymax, const double rx, const double ry) {
	unsigned n = pX.size();
	for (size_t row=0; row<nrows; row++) {
		double y = ymax - (startrow+row+0.5) * ry;
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



SpatRaster rasterizeLines(SpatVector p, SpatRaster r, std::vector<double> value, double background, SpatOptions &opt) {

	SpatRaster out = r.geometry(1);
  	if (!out.writeStart(opt)) { return out; }
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
			unsigned nln = line.size();
			for (size_t k = 0; k < nln; k++) {
				part = line.getPart(k);
				v = rasterize_line(v, value[j], part.x, part.y, out.bs.row[i], out.bs.nrows[i], out.ncol(), extent.xmin, extent.ymax, resx, resy);
			}
		}
		if (!out.writeValues(v, out.bs.row[i], out.bs.nrows[i], 0, out.ncol())) return out;
	}
	out.writeStop();
	return(out);
}


SpatRaster rasterizePoints(SpatVector p, SpatRaster r, std::vector<double> values, double background, SpatOptions &opt) {
	r.setError("not implemented yet");
	return(r);
}



SpatRaster SpatRaster::rasterize(SpatVector p, std::vector<double> values, double background, bool update, SpatOptions &opt) {
	std::string gtype = p.type();
	SpatRaster out = geometry(1);
	SpatOptions opts(opt);
	if (!update) {
		opts = opt;
	}
	if (gtype == "polygons") {
		out = rasterizePolygons(p, out, values, background, opts);
	} else if (gtype == "lines") {
		out = rasterizeLines(p, out, values, background, opts);
	}  else {
		out = rasterizePoints(p, out, values, background, opts);
	}
	if (update) out = cover(out, background, opt);
	return out;
}






SpatVector SpatRaster::as_points(bool values, bool narm) {

// for now assuming one layer

	BlockSize bs = getBlockSize(4);
	std::vector<double> v, vout;
	vout.reserve(v.size());
	SpatVector pv;
	SpatGeom g;
	g.gtype = points;


    std::vector<std::vector<double>> xy;
	if ((!values) && (!narm)) {
        double nc = ncell();
        for (size_t i=0; i<nc; i++) {
            xy = xyFromCell(i);
			SpatPart p(xy[0], xy[1]);
			g.addPart(p);
			pv.addGeom(g);
			g.parts.resize(0);
        }
		return pv;
	}

	if (values) {
        std::vector<std::string> nms = getNames();
        for (size_t i=0; i<nlyr(); i++) {
            pv.lyr.df.add_column(0, nms[i]);
        }
	}
	readStart();
	unsigned nc = ncol();
	unsigned nl = nlyr();
	for (size_t i = 0; i < bs.n; i++) {
		v = readValues(bs.row[i], bs.nrows[i], 0, nc);
        unsigned off1 = (bs.row[i] * nc);
 		unsigned vnc = bs.nrows[i] * nc;
		if (narm) {
            bool foundna = false;
			for (size_t j=0; j<vnc; j++) {
				for (size_t lyr=0; lyr<nl; lyr++) {
                    unsigned off2 = lyr*nc;
                    foundna = false;
                    if (std::isnan(v[off2+j])) {
                        foundna = true;
                        continue;
                    }
                }
                if (foundna) continue;
                xy = xyFromCell(off1+j);
                SpatPart p(xy[0], xy[1]);
                g.addPart(p);
                pv.addGeom(g);
                g.parts.resize(0);
                if (values) {
                    for (size_t lyr=0; lyr<nl; lyr++) {
                        unsigned off2 = lyr*nc;
                        pv.lyr.df.dv[lyr].push_back(v[off2+j]);
                    }
                }
			}
		} else { // if (values) {
			for (size_t j=0; j<vnc; j++) {
                xy = xyFromCell(off1+j);
                SpatPart p(xy[0], xy[1]);
                g.addPart(p);
                pv.addGeom(g);
                g.parts.resize(0);
                for (size_t lyr=0; lyr<nl; lyr++) {
                    unsigned off2 = lyr*nc;
                    pv.lyr.df.dv[lyr].push_back(v[off2+j]);
                }
			}
		}
	}
	readStop();
	return(pv);
}





void getCorners(std::vector<double> &x,  std::vector<double> &y, const double &X, const double &Y, const double &xr, const double &yr) {
	x[0] = X - xr;
	y[0] = Y - yr;
	x[1] = X - xr;
	y[1] = Y + yr;
	x[2] = X + xr;
	y[2] = Y + yr;
	x[3] = X + xr;
	y[3] = Y - yr;
	x[4] = x[0];
	y[4] = y[0];
}

/*
SpatVector SpatRaster::as_polygons(bool values, bool narm) {
	if (!values) narm=false;
	SpatVector v;
	SpatGeom g;
	g.gtype = polygons;
	double xr = xres()/2;
	double yr = yres()/2;
	std::vector<double> x(5);
	std::vector<double> y(5);
	if (!values) {
		std::vector<double> cells(ncell()) ;
		std::iota (std::begin(cells), std::end(cells), 0);
		std::vector< std::vector<double> > xy = xyFromCell(cells);
		for (size_t i=0; i<ncell(); i++) {
			getCorners(x, y, xy[0][i], xy[1][i], xr, yr);
			SpatPart p(x, y);
			g.addPart(p);
			v.addGeom(g);
			g.parts.resize(0);
		}
	} else {
		SpatRaster out = geometry();
		unsigned nl = nlyr();
		std::vector<std::vector<double> > att(ncell(), std::vector<double> (nl));

		BlockSize bs = getBlockSize(4);
		std::vector< std::vector<double> > xy;
		std::vector<double> atts(nl);
		for (size_t i=0; i<out.bs.n; i++) {
			std::vector<double> vals = readBlock(out.bs, i);
			unsigned nc=out.bs.nrows[i] * ncol();
			for (size_t j=0; j<nc; j++) {
				for (size_t k=0; k<nl; k++) {
					size_t kk = j + k * nl;
					att[nc+j][k] = vals[kk];
				}
				xy = xyFromCell(nc+j);
				getCorners(x, y, xy[0][0], xy[1][0], xr, yr);
				SpatPart p(x, y);
				g.addPart(p);
				v.addGeom(g);
				g.parts.resize(0);

			}
		}
		SpatDataFrame df;
		std::vector<std::string> nms = getNames();
		for (size_t i=0; i<att.size(); i++) {
			df.add_column(att[i], nms[i]);
		}
	}
	v.setCRS(getCRS());
	return(v);
}

*/

SpatVector SpatRaster::as_polygons(bool trunc, bool dissolve, bool values, bool narm) {

	if (!hasValues()) {
		values = false;
		narm = false;
		dissolve=false;
	}
	
	if (dissolve) {
		return polygonize(trunc);
	}

	SpatVector vect;
	if (!canProcessInMemory(12)) {
		vect.setError("the raster is too large");
		return vect;
	}

	bool remove_values = false;
	if (narm) {
		if (!values) remove_values = true;
		values=true;
	}

	unsigned nl = nlyr();
	unsigned nc = ncell();
	if (values) {
		std::vector<double> v = getValues();
		std::vector<std::string> nms = getNames();
		for (size_t i=0; i<nl; i++) {
			size_t offset = i * nc;
			std::vector<double> vv(v.begin()+offset, v.begin()+offset+nc);
			vect.add_column(vv, nms[i]);
		}
	}
	

	SpatGeom g;
	g.gtype = polygons;
	double xr = xres()/2;
	double yr = yres()/2;
	std::vector<double> x(5);
	std::vector<double> y(5);

	std::vector<double> cells(ncell()) ;
	std::iota (std::begin(cells), std::end(cells), 0);
	std::vector< std::vector<double> > xy = xyFromCell(cells);
	for (int i=nc-1; i>=0; i--) {
		if (narm) {
			bool erase = false;
			for (size_t j=0; j<nl; j++) {
				if (std::isnan(vect.lyr.df.dv[j][i])) {
					erase=true;
					break;
				}
			}
			if (erase) {
				for (size_t j=0; j<nl; j++) {
					vect.lyr.df.dv[j].erase (vect.lyr.df.dv[j].begin()+i);
				}
				continue; // skip the geom
			}
		}
		getCorners(x, y, xy[0][i], xy[1][i], xr, yr);
		SpatPart p(x, y);
		g.addPart(p);
		vect.addGeom(g);
		g.parts.resize(0);
	}

	std::reverse(std::begin(vect.lyr.geoms), std::end(vect.lyr.geoms));			

	if (remove_values) {
		vect.lyr.df = SpatDataFrame();		
	}
	vect.lyr.srs = srs;
	return(vect);
}




SpatVector SpatVector::as_lines() {
	SpatVector v = *this;
	if (lyr.geoms[0].gtype != polygons) {
		v.setError("this only works for polygons");
		return v;
	}
	for (size_t i=0; i<size(); i++) {
		for (size_t j=0; j < v.lyr.geoms[i].size(); j++) {
			SpatPart p = v.lyr.geoms[i].parts[j];
			if (p.hasHoles()) {
				for (size_t k=0; k < p.nHoles(); k++) {
					SpatHole h = p.getHole(k);
					SpatPart pp(h.x, h.y);
					v.lyr.geoms[i].addPart(pp);
				}
				p.holes.resize(0);
				v.lyr.geoms[i].parts[j] = p;
			}
		}
		v.lyr.geoms[i].gtype = lines;
	}
	return(v);
}


SpatVector SpatVector::as_points() {
	SpatVector v = *this;
	if (lyr.geoms[0].gtype == points) {
		v.addWarning("returning a copy");
		return v;
	}
	if (lyr.geoms[0].gtype == polygons) {
		v = v.as_lines();
	}

	for (size_t i=0; i<size(); i++) {
		v.lyr.geoms[i].gtype = points;
	}
	return(v);
}

