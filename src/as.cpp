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

#include <vector>
#include <cmath>
#include "spatRaster.h"

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

